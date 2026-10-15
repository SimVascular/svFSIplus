// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Fibers.h"

#include "ComMod.h"
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
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void normalize_vector(std::vector<double>& vector)
{
  double length = 0.0;
  for (const double value : vector) {
    length += value * value;
  }
  length = std::sqrt(length);

  if (utils::is_zero(length)) {
    throw std::runtime_error("A zero-length fiber direction was found.");
  }
  for (double& value : vector) {
    value /= length;
  }
}

Eigen::MatrixXd orthonormalize_fiber_matrix(
    const Eigen::Ref<const Eigen::MatrixXd>& input,
    const std::string& location)
{
  const int nsd = input.rows();
  const int num_input_fibers = input.cols();
  if (num_input_fibers == 1) {
    return input;
  }

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
      throw std::runtime_error(
        "Fiber direction " + std::to_string(direction + 1) +
        " is linearly dependent " + location + ".");
    }
    basis.col(direction) = candidate / norm;
  }

  // Complete an incomplete basis by applying Gram-Schmidt to the Cartesian
  // direction with the largest remaining orthogonal component.
  for (int direction = num_input_fibers; direction < nsd; ++direction) {
    Eigen::VectorXd best_candidate(nsd);
    double best_norm = 0.0;

    for (int axis = 0; axis < nsd; ++axis) {
      Eigen::VectorXd candidate = Eigen::VectorXd::Unit(nsd, axis);
      for (int pass = 0; pass < 2; ++pass) {
        for (int previous = 0; previous < direction; ++previous) {
          candidate -= candidate.dot(basis.col(previous)) * basis.col(previous);
        }
      }

      const double norm = candidate.norm();
      if (norm > best_norm) {
        best_candidate = candidate;
        best_norm = norm;
      }
    }

    if (best_norm <= linear_dependence_tolerance) {
      throw std::runtime_error("Failed to complete the fiber basis " + location + ".");
    }
    basis.col(direction) = best_candidate / best_norm;
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
    throw std::runtime_error("Nodal fibers must be a non-empty array.");
  }
  if (shape_functions.size() != num_nodes) {
    throw std::runtime_error("The number of shape functions must match the number of nodal fibers.");
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
    throw std::runtime_error("At least one shape function must be nonzero.");
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(interpolated_dyad);
  if (eigen_solver.info() != Eigen::Success) {
    throw std::runtime_error("Failed to compute the principal direction of the interpolated fiber dyad.");
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

  if (data_location == DataLocation::explicit_direction) {
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
    throw std::runtime_error("Fiber input data has no cell, point, or explicit direction association.");
  }

  Array3<double> reordered;
  if (cm.mas(cm_mod)) {
    reordered.resize(num_components, nG, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      reordered.set_slice(mesh.otnIEN[e], quadrature_fibers.slice(e));
    }
    quadrature_fibers.clear();
  }

  quadrature_fibers.resize(num_components, nG, nEl);
  for (int i = 0; i < num_proc; ++i) {
    disp[i] = mesh.eDist[i] * nG * num_components;
    scount[i] = (mesh.eDist[i + 1] - mesh.eDist[i]) * nG * num_components;
  }
  MPI_Scatterv(reordered.data(), scount.data(), disp.data(), cm_mod::mpreal,
      quadrature_fibers.data(), nEl * nG * num_components,
      cm_mod::mpreal, cm_mod.master, cm.com());
}


void Fibers::clear()
{
  input_fibers.clear();
  quadrature_fibers.clear();
  nfibers = 0;
  data_location = DataLocation::none;
}


void Fibers::read(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  const int num_paths = mesh_parameters.fiber_direction_file_paths.size();
  clear();

  const int num_dirs = mesh_parameters.fiber_directions.size();
  if (num_paths == 0 && num_dirs == 0) {
    return;
  }

  const int num_input_fibers = num_paths != 0 ? num_paths : num_dirs;
  if (num_input_fibers > nsd) {
    throw std::runtime_error(
      "The number of input fiber directions (" +
      std::to_string(num_input_fibers) +
      ") exceeds the spatial dimension (" + std::to_string(nsd) + ").");
  }

  if (num_paths != 0) {
    read_from_file(mesh_parameters, mesh, nsd);
  } else {
    read_directions(mesh_parameters, nsd);
  }

  if (data_location == DataLocation::point) {
    populate_quadrature_fibers(mesh, nsd);
  } else if (data_location == DataLocation::cell ||
             data_location == DataLocation::explicit_direction) {
    orthonormalize_fibers(input_fibers, nsd);
  } else {
    throw std::runtime_error(
      "Fiber input data has no cell, point, or explicit direction association.");
  }
}


Array<double> Fibers::get_quadrature_fibers(
  const int nsd, const int g, const int e) const
{
  if (nsd <= 0) {
    throw std::runtime_error("The spatial dimension must be positive when accessing fiber directions.");
  }

  double* data = nullptr;
  int num_components = 0;

  if (data_location == DataLocation::explicit_direction) {
    num_components = input_fibers.nrows();
    data = input_fibers.data();
  } else if (data_location == DataLocation::cell) {
    if (e < 0 || e >= input_fibers.ncols()) {
      throw std::runtime_error(
        "The requested fiber element index " + std::to_string(e) +
        " is out of bounds.");
    }
    num_components = input_fibers.nrows();
    data = input_fibers.data() + e * num_components;
  } else if (data_location == DataLocation::point) {
    if (g < 0 || g >= quadrature_fibers.ncols() ||
        e < 0 || e >= quadrature_fibers.nslices()) {
      throw std::runtime_error(
        "The requested quadrature fiber indices (" + std::to_string(g) +
        ", " + std::to_string(e) + ") are out of bounds.");
    }
    num_components = quadrature_fibers.nrows();
    const auto offset =
      (e * quadrature_fibers.ncols() + g) * num_components;
    data = quadrature_fibers.data() + offset;
  } else {
    throw std::runtime_error("Fiber input data has not been initialized.");
  }

  if (num_components == 0 || num_components % nsd != 0) {
    throw std::runtime_error(
      "The fiber data have " + std::to_string(num_components) +
      " components, which is incompatible with spatial dimension " +
      std::to_string(nsd) + ".");
  }
  return {nsd, num_components / nsd, data};
}


Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
Fibers::get_eigen_fibers(const int nsd, const int g, const int e) const
{
  const auto fibers = get_quadrature_fibers(nsd, g, e);
  return {fibers.data(), fibers.nrows(), fibers.ncols()};
}


void Fibers::read_from_file(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  const int num_paths = mesh_parameters.fiber_direction_file_paths.size();
  nfibers = num_paths;

  const auto paths = mesh_parameters.fiber_direction_file_paths();
  for (int idx = 0; idx < nfibers; ++idx) {
    const auto& file_name = paths[idx];

    if (FILE* file = fopen(file_name.c_str(), "r")) {
      fclose(file);
    } else {
      throw std::runtime_error("The fiber direction VTK file '" + file_name + "' can't be read.");
    }

    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file_name.c_str());
    reader->Update();
    auto grid = reader->GetOutput();

    // Check that the number of elements in the fiber direction VTK file is equal to the number of elements in the mesh.
    if (mesh.gnEl != grid->GetNumberOfCells()) {
      throw std::runtime_error("The number of elements (" + std::to_string(grid->GetNumberOfCells()) +
          ") in the fiber direction VTK file '" + file_name + "' is not equal to the number of elements (" +
          std::to_string(mesh.gnEl) + ") for the mesh named '" + mesh.name + "'.");
    }

    // Get the fiber orientation data and determine whether it is associated
    // with cells or points.
    auto array = grid->GetCellData()->GetArray("FIB_DIR");
    DataLocation current_location = DataLocation::cell;
    vtkIdType num_fibers = grid->GetNumberOfCells();
    if (array == nullptr) {
      array = grid->GetPointData()->GetArray("FIB_DIR");
      current_location = DataLocation::point;
      num_fibers = grid->GetNumberOfPoints();
    }
    if (array == nullptr) {
      throw std::runtime_error("No 'FIB_DIR' data found in the cell data or point data of the fiber direction VTK file '" +
          file_name + "'");
    }

    // Check that the number of tuples and components in the fiber orientation data is consistent with the mesh.
    if (array->GetNumberOfTuples() != num_fibers) {
      throw std::runtime_error("The number of 'FIB_DIR' tuples (" + std::to_string(array->GetNumberOfTuples()) +
          ") in the fiber direction VTK file '" + file_name + "' does not match the number of associated cells or points (" +
          std::to_string(num_fibers) + ").");
    }
    if (array->GetNumberOfComponents() < nsd) {
      throw std::runtime_error("The 'FIB_DIR' data in the fiber direction VTK file '" + file_name + "' has fewer than " +
          std::to_string(nsd) + " components.");
    }

    const vtkIdType expected_num_fibers = current_location == DataLocation::cell ? mesh.gnEl : mesh.gnNo;
    if (num_fibers != expected_num_fibers) {
      const std::string association = current_location == DataLocation::cell ? "cells" : "points";
      throw std::runtime_error("The number of " + association + " (" + std::to_string(num_fibers) +
          ") in the fiber direction VTK file '" + file_name + "' does not match the mesh (" +
          std::to_string(expected_num_fibers) + ").");
    }

    if (idx == 0) {
      data_location = current_location;
      if (data_location == DataLocation::cell) {
        const int num_directions = nfibers == 1 ? 1 : nsd;
        input_fibers.resize(num_directions * nsd, num_fibers);
      } else {
        input_fibers.resize(num_paths * nsd, num_fibers);
      }
      input_fibers = 0.0;
    } else if (data_location != current_location) {
      throw std::runtime_error("The 'FIB_DIR' arrays in the fiber direction VTK files must all be associated with "
          "either cells or points.");
    }

    // Set the fiber orientations data.
    const int offset = idx * nsd;
    std::vector<double> direction(nsd);
    for (vtkIdType fiber_index = 0; fiber_index < num_fibers; ++fiber_index) {
      const auto tuple = array->GetTuple(fiber_index);
      for (int component = 0; component < nsd; ++component) {
        direction[component] = tuple[component];
      }
      normalize_vector(direction);
      for (int component = 0; component < nsd; ++component) {
        input_fibers(offset + component, fiber_index) = direction[component];
      }
    }
  }
}


void Fibers::read_directions(const MeshParameters& mesh_parameters, const int nsd)
{
  const int num_dirs = mesh_parameters.fiber_directions.size();
  nfibers = num_dirs;
  data_location = DataLocation::explicit_direction;

  const int num_directions = nfibers == 1 ? 1 : nsd;
  input_fibers.resize(num_directions * nsd, 1);
  input_fibers = 0.0;
  const auto directions = mesh_parameters.fiber_directions;
  for (int i = 0; i < nfibers; ++i) {
    auto direction = directions[i]();
    normalize_vector(direction);
    for (int j = 0; j < nsd; ++j) {
      input_fibers(i * nsd + j, 0) = direction[j];
    }
  }
}

void Fibers::populate_quadrature_fibers(const mshType& mesh, const int nsd)
{
  if (mesh.nG <= 0) {
    throw std::runtime_error(
      "Point fiber directions cannot be populated before mesh quadrature points are initialized.");
  }
  const int num_directions = nfibers == 1 ? 1 : nsd;
  quadrature_fibers.resize(num_directions * nsd, mesh.nG, mesh.gnEl);
  quadrature_fibers = 0.0;
  interpolate_point_fibers(mesh, nsd);
  orthonormalize_quadrature_fibers(nsd);
}


void Fibers::orthonormalize_quadrature_fibers(const int nsd)
{
  // One direction does not define a full basis and must remain unchanged.
  if (nfibers == 1) {
    return;
  }

  const int nG = quadrature_fibers.ncols();
  const int nEl = quadrature_fibers.nslices();

  for (int e = 0; e < nEl; ++e) {
    for (int g = 0; g < nG; ++g) {
      const int offset = (e * nG + g) * nsd * nsd;
      const Eigen::Map<const Eigen::MatrixXd> input(
        quadrature_fibers.data() + offset, nsd, nfibers);
      const auto basis = orthonormalize_fiber_matrix(
        input, "at quadrature point " + std::to_string(g) +
        " of element " + std::to_string(e));
      Eigen::Map<Eigen::MatrixXd> output(
        quadrature_fibers.data() + offset, nsd, nsd);
      output = basis;
    }
  }
}


void Fibers::orthonormalize_fibers(Array<double>& fibers, const int nsd)
{
  // One direction does not define a full basis and is already normalized.
  if (nfibers == 1) {
    return;
  }

  const int num_locations = fibers.ncols();
  for (int location = 0; location < num_locations; ++location) {
    const Eigen::Map<const Eigen::MatrixXd> input(
      fibers.data() + location * fibers.nrows(), nsd, nfibers);
    const std::string description = data_location == DataLocation::cell
      ? "in element " + std::to_string(location)
      : "in the explicit directions";
    const auto basis = orthonormalize_fiber_matrix(input, description);
    Eigen::Map<Eigen::MatrixXd> output(
      fibers.data() + location * fibers.nrows(), nsd, nsd);
    output = basis;
  }
}


void Fibers::interpolate_point_fibers(const mshType& mesh, const int nsd)
{
  if (mesh.N.nrows() != mesh.eNoN || mesh.N.ncols() != mesh.nG) {
    throw std::runtime_error("Mesh shape functions are inconsistent with its nodes and quadrature points.");
  }

  Array<double> element_nodal_fibers(nsd, mesh.eNoN);
  Vector<double> shape_functions(mesh.eNoN);

  // Loop over the quadrature points.
  for (int g = 0; g < mesh.nG; ++g) {
    shape_functions = mesh.N.col(g);

    // Loop over the elements.
    for (int e = 0; e < mesh.gnEl; ++e) {

      // Loop over the fiber directions.
      for (int iFn = 0; iFn < nfibers; ++iFn) {

        // Build the element nodal fibers.
        for (int a = 0; a < mesh.eNoN; ++a) {
          const int node = mesh.gIEN(a, e);
          for (int i = 0; i < nsd; ++i) {
            element_nodal_fibers(i, a) = input_fibers(i + iFn * nsd, node);
          }
        }

        // Interpolate the fiber directions to the quadrature point.
        const auto interpolated = interpolate_dyadic(element_nodal_fibers, shape_functions);
        for (int i = 0; i < nsd; ++i) {
          quadrature_fibers(i + iFn * nsd, g, e) = interpolated[i];
        }
      }
    }
  }
}
