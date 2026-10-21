// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Fibers.h"

#include "ComMod.h"
#include "Core/Exception.h"
#include "FE/Common/FEException.h"
#include "Parameters.h"
#include "utils.h"

#include "eigen3/Eigen/Eigenvalues"

#include <mpi.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <cmath>
#include <string>

namespace {

Eigen::MatrixXd orthonormalize_fiber_matrix(
    const Eigen::Ref<const Eigen::MatrixXd>& input,
    const std::string& location)
{
  const int nsd = input.rows();
  const int num_input_fibers = input.cols();

  Eigen::MatrixXd basis = Eigen::MatrixXd::Zero(nsd, nsd);
  constexpr double linear_dependence_tolerance = 1.0e-12;

  // Modified Gram-Schmidt. The first direction is copied without
  // modification; every later input direction is made orthogonal to it and
  // to the previously accepted directions.
  basis.col(0) = input.col(0);
  for (int direction = 1; direction < num_input_fibers; ++direction) {
    Eigen::VectorXd candidate = input.col(direction);
    for (int pass = 0; pass < 2; ++pass) {
      for (int previous = 0; previous < direction; ++previous) {
        candidate -= candidate.dot(basis.col(previous)) * basis.col(previous);
      }
    }

    const double norm = candidate.norm();
    if (norm <= linear_dependence_tolerance) {
      svmp::raise<svmp::ParseException>(
        "Fiber direction " + std::to_string(direction + 1) +
        " is linearly dependent " + location + ".");
    }
    basis.col(direction) = candidate / norm;
  }

  // In 3D, the normal must be the cross product of the orthonormalized fiber
  // and sheet directions so that the completed frame is right-handed.
  if (nsd == 3) {
    const Eigen::Vector3d fiber = basis.col(0);
    const Eigen::Vector3d sheet = basis.col(1);
    basis.col(2) = fiber.cross(sheet);
  } else if (nsd == 2) {
    basis(0, 1) = -basis(1, 0);
    basis(1, 1) = basis(0, 0);
  } else {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "Fiber basis completion is only supported for spatial dimensions 2 and 3.");
  }

  return basis;
}

// Interpolates one nodal fiber direction to a quadrature point. The dyads are
// interpolated before the principal direction is recovered, so opposite
// representations of the same fiber direction do not cancel.
Vector<double> interpolate_dyadic(
    const Array<double>& nodal_fibers,
    const Vector<double>& shape_functions)
{
  const int nsd = nodal_fibers.nrows();
  const int num_nodes = nodal_fibers.ncols();
  if (nsd == 0 || num_nodes == 0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "Nodal fibers must be a non-empty array.");
  }
  if (shape_functions.size() != num_nodes) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The number of shape functions must match the number of nodal fibers.");
  }

  Eigen::MatrixXd interpolated_dyad = Eigen::MatrixXd::Zero(nsd, nsd);
  Eigen::VectorXd reference_direction = Eigen::VectorXd::Zero(nsd);
  double reference_weight = -1.0;

  for (int node = 0; node < num_nodes; ++node) {
    const double shape_value = shape_functions[node];
    if (utils::is_zero(shape_value)) {
      continue;
    }

    Eigen::VectorXd direction(nsd);
    for (int component = 0; component < nsd; ++component) {
      direction[component] = nodal_fibers(component, node);
    }

    interpolated_dyad.noalias() += shape_value * direction * direction.transpose();

    if (std::abs(shape_value) > reference_weight) {
      reference_direction = direction;
      reference_weight = std::abs(shape_value);
    }
  }

  if (reference_weight < 0.0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "At least one shape function must be nonzero.");
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(interpolated_dyad);
  if (eigen_solver.info() != Eigen::Success) {
    svmp::raise<svmp::InternalErrorException>(
      "Failed to compute the principal direction of the interpolated fiber dyad.");
  }

  Eigen::VectorXd interpolated_fiber = eigen_solver.eigenvectors().col(nsd - 1);
  const double interpolated_fiber_norm = interpolated_fiber.norm();
  if (!utils::is_zero(interpolated_fiber_norm, 1.0)) {
    interpolated_fiber /= interpolated_fiber_norm;
  }
  if (interpolated_fiber.dot(reference_direction) < 0.0) {
    interpolated_fiber = -interpolated_fiber;
  }

  Vector<double> result(nsd);
  for (int component = 0; component < nsd; ++component) {
    result[component] = interpolated_fiber[component];
  }
  return result;
}

} // namespace


void Fibers::distribute(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm,
    mshType& mesh)
{
  const int nG = mesh.nG;
  const int num_proc = cm.np();
  const int nEl = mesh.nEl;
  int location = static_cast<int>(data_location);
  cm.bcast(cm_mod, &location);
  data_location = static_cast<DataLocation>(location);

  if (!has_fibers()) {
    return;
  }

  int num_components = 0;
  if (cm.mas(cm_mod)) {
    num_components = data_location == DataLocation::point
      ? quadrature_fibers.nrows()
      : input_fibers.nrows();
  }
  cm.bcast(cm_mod, &num_components);

  if (data_location == DataLocation::constant) {
    if (cm.slv(cm_mod)) {
      input_fibers.resize(num_components, 1);
    }
    MPI_Bcast(input_fibers.data(), num_components, cm_mod::mpreal,
      cm_mod.master, cm.com());
    return;
  }

  Vector<int> scount(num_proc), disp(num_proc);

  if (data_location == DataLocation::cell) {
    Array<double> reordered;
    if (cm.mas(cm_mod)) {
      reordered.resize(num_components, mesh.gnEl);
      for (int e = 0; e < mesh.gnEl; ++e) {
        reordered.set_col(mesh.otnIEN[e], input_fibers.col(e));
      }
      input_fibers.clear();
    }

    input_fibers.resize(num_components, nEl);
    for (int i = 0; i < num_proc; ++i) {
      disp[i] = mesh.eDist[i] * num_components;
      scount[i] = (mesh.eDist[i + 1] - mesh.eDist[i]) * num_components;
    }
    MPI_Scatterv(reordered.data(), scount.data(), disp.data(), cm_mod::mpreal,
      input_fibers.data(), nEl * num_components, cm_mod::mpreal,
      cm_mod.master, cm.com());
    return;
  }

  if (data_location != DataLocation::point) {
    svmp::raise<svmp::InternalErrorException>(
      "Fiber input data has no cell, point, or constant direction association.");
  }

  Array3<double> reordered;
  Array3<double> reordered_nodal;
  if (cm.mas(cm_mod)) {
    reordered.resize(num_components, nG, mesh.gnEl);
    reordered_nodal.resize(num_components, mesh.eNoN, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      reordered.set_slice(mesh.otnIEN[e], quadrature_fibers.slice(e));
      reordered_nodal.set_slice(
        mesh.otnIEN[e], element_nodal_fibers.slice(e));
    }
    quadrature_fibers.clear();
    element_nodal_fibers.clear();
  }

  quadrature_fibers.resize(num_components, nG, nEl);
  for (int i = 0; i < num_proc; ++i) {
    disp[i] = mesh.eDist[i] * nG * num_components;
    scount[i] = (mesh.eDist[i + 1] - mesh.eDist[i]) * nG * num_components;
  }
  MPI_Scatterv(reordered.data(), scount.data(), disp.data(), cm_mod::mpreal,
      quadrature_fibers.data(), nEl * nG * num_components,
      cm_mod::mpreal, cm_mod.master, cm.com());

  element_nodal_fibers.resize(num_components, mesh.eNoN, nEl);
  for (int i = 0; i < num_proc; ++i) {
    disp[i] = mesh.eDist[i] * mesh.eNoN * num_components;
    scount[i] = (mesh.eDist[i + 1] - mesh.eDist[i]) *
      mesh.eNoN * num_components;
  }
  MPI_Scatterv(reordered_nodal.data(), scount.data(), disp.data(),
      cm_mod::mpreal, element_nodal_fibers.data(),
      nEl * mesh.eNoN * num_components, cm_mod::mpreal,
      cm_mod.master, cm.com());
}


Array<double> Fibers::gather_fibers(
  const ComMod& com_mod, const CmMod& cm_mod, const mshType& mesh,
  const int nsd) const
{
  Array<double> gathered;

  if (data_location == DataLocation::none) {
    return gathered;
  }

  const int num_components = nsd * nsd;
  auto& cm = com_mod.cm;

  if (data_location == DataLocation::point) {
    if (cm.mas(cm_mod)) {
      if (input_fibers.nrows() != num_components ||
          input_fibers.ncols() != mesh.gnNo) {
        svmp::raise<svmp::InternalErrorException>(
          "Point fiber basis dimensions are inconsistent with the mesh.");
      }
      gathered = input_fibers;
    }
    return gathered;
  }

  if (data_location != DataLocation::cell &&
      data_location != DataLocation::constant) {
    return gathered;
  }

  if (mesh.nEl != 0 &&
      (input_fibers.nrows() != num_components ||
       (input_fibers.ncols() != 1 &&
        input_fibers.ncols() != mesh.nEl))) {
    svmp::raise<svmp::InternalErrorException>(
      "Cell fiber basis dimensions are inconsistent with the mesh.");
  }

  Array<double> local_basis(num_components, mesh.nEl);
  for (int e = 0; e < mesh.nEl; ++e) {
    const int input_location = input_fibers.ncols() == 1 ? 0 : e;
    local_basis.set_col(e, input_fibers.col(input_location));
  }

  Vector<int> counts(cm.np());
  Vector<int> displacements(cm.np());
  for (int rank = 0; rank < cm.np(); ++rank) {
    displacements[rank] = mesh.eDist[rank] * num_components;
    counts[rank] =
      (mesh.eDist[rank + 1] - mesh.eDist[rank]) * num_components;
  }

  Array<double> partitioned_basis;
  if (cm.mas(cm_mod)) {
    partitioned_basis.resize(num_components, mesh.gnEl);
  }
  MPI_Gatherv(local_basis.data(), mesh.nEl * num_components,
    cm_mod::mpreal, partitioned_basis.data(), counts.data(),
    displacements.data(), cm_mod::mpreal, cm_mod.master, cm.com());

  if (cm.mas(cm_mod)) {
    gathered.resize(num_components, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      gathered.set_col(e, partitioned_basis.col(mesh.otnIEN[e]));
    }
  }
  return gathered;
}


void Fibers::clear()
{
  input_fibers.clear();
  quadrature_fibers.clear();
  element_nodal_fibers.clear();
  nfibers = 0;
  data_location = DataLocation::none;
}


void Fibers::read(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  // Count the number of input modes.
  const int num_point_paths =
      mesh_parameters.fiber_direction_file_paths_points.size();
  const int num_cell_paths =
      mesh_parameters.fiber_direction_file_paths_cells.size();
  const int num_dirs = mesh_parameters.fiber_directions.size();

  // If no input modes are specified it means no fibers are needed, return early.
  const int num_input_modes =
      (num_point_paths != 0) + (num_cell_paths != 0) + (num_dirs != 0);
  if (num_input_modes == 0) {
    return;
  }

  // If multiple input modes are specified, throw an error.
  if (num_input_modes > 1) {
    svmp::raise<svmp::ParseException>(
      "Specify fiber directions with only one of Fiber_direction_constant, "
      "Fiber_direction_from_file_points, or Fiber_direction_from_file_cells.");
  }

  // Require exactly nsd - 1 input directions; the last basis vector is completed
  // during orthonormalization.
  const int num_input_fibers = num_dirs + num_point_paths + num_cell_paths;
  const int expected_fibers = nsd - 1;
  if (num_input_fibers != expected_fibers) {
    svmp::raise<svmp::ParseException>(
      "The number of input fiber directions (" +
      std::to_string(num_input_fibers) +
      ") must be exactly nsd - 1 (" + std::to_string(expected_fibers) + ").");
  }

  if (num_dirs != 0) {
    read_directions(mesh_parameters, nsd);
  } else if (num_point_paths != 0) {
    read_from_file(mesh_parameters, mesh, nsd, DataLocation::point);
  } else {
    read_from_file(mesh_parameters, mesh, nsd, DataLocation::cell);
  }

  check_input_fibers_orthonormal(nsd);
  complete_fiber_basis(input_fibers, nsd);

  if (data_location == DataLocation::point) {
    populate_quadrature_fibers(mesh, nsd);
  } else if (data_location != DataLocation::cell &&
             data_location != DataLocation::constant) {
    svmp::raise<svmp::InternalErrorException>(
      "Fiber input data has no cell, point, or constant direction association.");
  }
}


void Fibers::read_from_file(const MeshParameters& mesh_parameters, mshType& mesh,
  const int nsd, const DataLocation expected_location)
{
  if (expected_location != DataLocation::point &&
      expected_location != DataLocation::cell) {
    svmp::raise<svmp::InternalErrorException>(
      "Fiber direction files must be associated with cells or points.");
  }

  const auto& path_parameter = expected_location == DataLocation::point
      ? mesh_parameters.fiber_direction_file_paths_points
      : mesh_parameters.fiber_direction_file_paths_cells;
  const int num_paths = path_parameter.size();
  nfibers = num_paths;
  data_location = expected_location;

  const auto paths = path_parameter();
  const bool expect_cells = expected_location == DataLocation::cell;
  const std::string association = expect_cells ? "cells" : "points";
  const std::string xml_tag = expect_cells
      ? "Fiber_direction_from_file_cells"
      : "Fiber_direction_from_file_points";

  for (int idx = 0; idx < nfibers; ++idx) {
    const auto& file_name = paths[idx];

    if (FILE* file = fopen(file_name.c_str(), "r")) {
      fclose(file);
    } else {
      svmp::raise<svmp::FileNotFoundException>(file_name);
    }

    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file_name.c_str());
    reader->Update();
    auto grid = reader->GetOutput();

    if (mesh.gnEl != grid->GetNumberOfCells()) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The number of elements (" +
        std::to_string(grid->GetNumberOfCells()) +
        ") is not equal to the number of elements (" +
        std::to_string(mesh.gnEl) + ") for the mesh named '" + mesh.name + "'.");
    }

    vtkDataArray* array = nullptr;
    vtkIdType num_locations = 0;
    if (expect_cells) {
      array = grid->GetCellData()->GetArray("FIB_DIR");
      num_locations = grid->GetNumberOfCells();
    } else {
      array = grid->GetPointData()->GetArray("FIB_DIR");
      num_locations = grid->GetNumberOfPoints();
    }
    if (array == nullptr) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "No 'FIB_DIR' " + association + " data found as required by " +
        xml_tag + ".");
    }
    if (array->GetNumberOfTuples() != num_locations) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The number of 'FIB_DIR' tuples (" +
        std::to_string(array->GetNumberOfTuples()) +
        ") does not match the number of associated " + association + " (" +
        std::to_string(num_locations) + ").");
    }
    if (array->GetNumberOfComponents() < nsd) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The 'FIB_DIR' data has fewer than " + std::to_string(nsd) +
        " components.");
    }

    const vtkIdType expected_num_locations = expect_cells ? mesh.gnEl : mesh.gnNo;
    if (num_locations != expected_num_locations) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The number of " + association + " (" +
        std::to_string(num_locations) + ") does not match the mesh (" +
        std::to_string(expected_num_locations) + ").");
    }

    if (idx == 0) {
      // Allocate space for the completed nsd-dimensional basis. Input directions
      // fill the first nfibers columns; orthonormalization completes the rest.
      input_fibers.resize(nsd * nsd, num_locations);
      input_fibers = 0.0;
    }

    const int offset = idx * nsd;
    for (vtkIdType location = 0; location < num_locations; ++location) {
      const auto tuple = array->GetTuple(location);
      for (int component = 0; component < nsd; ++component) {
        input_fibers(offset + component, location) = tuple[component];
      }
    }
  }
}


void Fibers::read_directions(const MeshParameters& mesh_parameters, const int nsd)
{
  const auto& directions = mesh_parameters.fiber_directions;
  const int num_dirs = directions.size();

  for (int i = 0; i < num_dirs; ++i) {
    const int num_components = directions[i].size();
    if (num_components != nsd) {
      svmp::raise<svmp::ParseException>(
        "Fiber_direction_constant " + std::to_string(i + 1) + " has " +
        std::to_string(num_components) + " components, but exactly " +
        std::to_string(nsd) + " are required for a " +
        std::to_string(nsd) + "D mesh.");
    }
  }

  nfibers = num_dirs;
  data_location = DataLocation::constant;

  // Allocate space for the completed nsd-dimensional basis.
  input_fibers.resize(nsd * nsd, 1);
  input_fibers = 0.0;
  for (int i = 0; i < nfibers; ++i) {
    for (int j = 0; j < nsd; ++j) {
      input_fibers(i * nsd + j, 0) = directions[i][j];
    }
  }
}


void Fibers::check_input_fibers_orthonormal(const int nsd) const
{
  constexpr double tolerance = 1.0e-6;

  // Check that the fiber directions are unit vectors.
  for (int location = 0; location < input_fibers.ncols(); ++location) {
    for (int i = 0; i < nfibers; ++i) {
      double norm_sq = 0.0;
      for (int component = 0; component < nsd; ++component) {
        const double value = input_fibers(i * nsd + component, location);
        norm_sq += value * value;
      }
      const double norm = std::sqrt(norm_sq);
      if (std::abs(norm - 1.0) > tolerance) {
        svmp::raise<svmp::ParseException>(
          "Fiber direction " + std::to_string(i + 1) + " at location " +
          std::to_string(location) + " is not a unit vector (norm = " +
          std::to_string(norm) + ").");
      }
    }
  }

  // Check that the fiber directions are orthogonal.
  const int nlocations = input_fibers.ncols();
  for (int location = 0; location < nlocations; ++location) {
    for (int i = 0; i < nfibers; ++i) {
      for (int j = 0; j < i; ++j) {
        double dot = 0.0;
        for (int component = 0; component < nsd; ++component) {
          dot += input_fibers(i * nsd + component, location) *
                 input_fibers(j * nsd + component, location);
        }
        if (std::abs(dot) > tolerance) {
          svmp::raise<svmp::ParseException>(
            "Fiber directions " + std::to_string(j + 1) + " and " +
            std::to_string(i + 1) + " at location " +
            std::to_string(location) + " are not orthogonal (dot product = " +
            std::to_string(dot) + ").");
        }
      }
    }
  }
}


Array<double> Fibers::get_quadrature_fibers(
  const int nsd, const int g, const int e) const
{
  if (nsd <= 0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The spatial dimension must be positive when accessing fiber directions.");
  }

  double* data = nullptr;
  int num_components = 0;

  if (data_location == DataLocation::constant) {
    num_components = input_fibers.nrows();
    data = input_fibers.data();
  } else if (data_location == DataLocation::cell) {
    if (e < 0 || e >= input_fibers.ncols()) {
      svmp::raise<svmp::IndexOutOfRangeException>(
        "The requested fiber element index " + std::to_string(e) +
        " is out of bounds.");
    }
    num_components = input_fibers.nrows();
    data = input_fibers.data() + e * num_components;
  } else if (data_location == DataLocation::point) {
    if (g < 0 || g >= quadrature_fibers.ncols() ||
        e < 0 || e >= quadrature_fibers.nslices()) {
      svmp::raise<svmp::IndexOutOfRangeException>(
        "The requested quadrature fiber indices (" + std::to_string(g) +
        ", " + std::to_string(e) + ") are out of bounds.");
    }
    num_components = quadrature_fibers.nrows();
    const auto offset =
      (e * quadrature_fibers.ncols() + g) * num_components;
    data = quadrature_fibers.data() + offset;
  } else {
    svmp::raise<svmp::InternalErrorException>(
      "Fiber input data has not been initialized.");
  }

  if (num_components == 0 || num_components % nsd != 0) {
    svmp::raise<svmp::InternalErrorException>(
      "The fiber data have " + std::to_string(num_components) +
      " components, which is incompatible with spatial dimension " +
      std::to_string(nsd) + ".");
  }
  return {nsd, num_components / nsd, data};
}


Array<double> Fibers::get_quadrature_fibers(
  const int nsd, const Vector<double>& basis_functions, const int e) const
{
  if (data_location != DataLocation::point) {
    return get_quadrature_fibers(nsd, 0, e);
  }
  if (nsd <= 0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The spatial dimension must be positive when accessing fiber directions.");
  }
  if (e < 0 || e >= element_nodal_fibers.nslices()) {
    svmp::raise<svmp::IndexOutOfRangeException>(
      "The requested fiber element index " + std::to_string(e) +
      " is out of bounds.");
  }
  const int num_basis_functions = basis_functions.size();
  if (num_basis_functions <= 0 ||
      num_basis_functions > element_nodal_fibers.ncols()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The number of basis-function values must be positive and cannot exceed "
      "the number of element nodes.");
  }

  Array<double> interpolated_fibers(nsd, nfibers);
  Array<double> nodal_fibers(nsd, num_basis_functions);

  for (int iFn = 0; iFn < nfibers; ++iFn) {
    for (int a = 0; a < num_basis_functions; ++a) {
      for (int component = 0; component < nsd; ++component) {
        nodal_fibers(component, a) =
          element_nodal_fibers(component + iFn * nsd, a, e);
      }
    }

    const auto interpolated =
      interpolate_dyadic(nodal_fibers, basis_functions);
    for (int component = 0; component < nsd; ++component) {
      interpolated_fibers(component, iFn) = interpolated[component];
    }
  }

  const Eigen::Map<const Eigen::MatrixXd> input(
    interpolated_fibers.data(), nsd, nfibers);
  const auto basis = orthonormalize_fiber_matrix(
    input, "in element " + std::to_string(e));

  Array<double> result(nsd, nsd);
  Eigen::Map<Eigen::MatrixXd>(result.data(), nsd, nsd) = basis;
  return result;
}


Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
Fibers::get_eigen_fibers(const int nsd, const int g, const int e) const
{
  const auto fibers = get_quadrature_fibers(nsd, g, e);
  return {fibers.data(), fibers.nrows(), fibers.ncols()};
}


void Fibers::populate_quadrature_fibers(const mshType& mesh, const int nsd)
{
  if (mesh.nG <= 0) {
    svmp::raise<svmp::InternalErrorException>(
      "Point fiber directions cannot be populated before mesh quadrature points are initialized.");
  }
  quadrature_fibers.resize(nsd * nsd, mesh.nG, mesh.gnEl);
  quadrature_fibers = 0.0;
  interpolate_point_fibers_mesh(mesh, nsd);
}


void Fibers::interpolate_point_fibers_mesh(
  const mshType& mesh, const int nsd)
{
  if (mesh.N.nrows() != mesh.eNoN || mesh.N.ncols() != mesh.nG) {
    svmp::raise<svmp::InternalErrorException>(
      "Mesh shape functions are inconsistent with its nodes and quadrature points.");
  }

  element_nodal_fibers.resize(nsd * nsd, mesh.eNoN, mesh.gnEl);
  element_nodal_fibers = 0.0;
  for (int e = 0; e < mesh.gnEl; ++e) {
    for (int a = 0; a < mesh.eNoN; ++a) {
      const int node = mesh.gIEN(a, e);
      for (int direction = 0; direction < nfibers; ++direction) {
        for (int component = 0; component < nsd; ++component) {
          element_nodal_fibers(component + direction * nsd, a, e) =
            input_fibers(component + direction * nsd, node);
        }
      }
    }
  }

  for (int g = 0; g < mesh.nG; ++g) {
    const auto basis_functions = mesh.N.col(g);
    for (int e = 0; e < mesh.gnEl; ++e) {
      const auto fibers =
        get_quadrature_fibers(nsd, basis_functions, e);
      for (int direction = 0; direction < nsd; ++direction) {
        for (int component = 0; component < nsd; ++component) {
          quadrature_fibers(component + direction * nsd, g, e) =
            fibers(component, direction);
        }
      }
    }
  }
}


void Fibers::complete_fiber_basis(Array<double>& fibers, const int nsd)
{
  // The input directions are already checked to be orthonormal, so the
  // remaining basis vector is just the cross product (3D) or a 90-degree
  // rotation (2D). Point data still uses Gram-Schmidt after interpolation.
  const int num_locations = fibers.ncols();
  for (int location = 0; location < num_locations; ++location) {
    Eigen::Map<Eigen::MatrixXd> basis(
      fibers.data() + location * fibers.nrows(), nsd, nsd);

    if (nsd == 3) {
      const Eigen::Vector3d fiber = basis.col(0);
      const Eigen::Vector3d sheet = basis.col(1);
      basis.col(2) = fiber.cross(sheet);
    } else if (nsd == 2) {
      basis(0, 1) = -basis(1, 0);
      basis(1, 1) = basis(0, 0);
    } else if (nsd != 1) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
        "Fiber basis completion is only supported for spatial dimensions 1, 2, or 3.");
    }
  }
}
