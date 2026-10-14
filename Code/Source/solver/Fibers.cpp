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

} // namespace


void Fibers::distribute(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm,
    mshType& mesh)
{
  const int nsd = com_mod.nsd;
  const int nG = mesh.nG;
  const int num_proc = cm.np();
  const int nEl = mesh.nEl;
  bool present = false;
  Array3<double> reordered;
  Vector<int> scount(num_proc), disp(num_proc);

  if (cm.mas(cm_mod) && has_fibers()) {
    present = true;
    reordered.resize(nfibers * nsd, nG, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      reordered.set_slice(mesh.otnIEN[e], quadrature_fibers.slice(e));
    }
    quadrature_fibers.clear();
  }
  cm.bcast(cm_mod, &present);
  if (!present) {
    return;
  }

  quadrature_fibers.resize(nfibers * nsd, nG, nEl);
  for (int i = 0; i < num_proc; ++i) {
    disp[i] = mesh.eDist[i] * nG * nfibers * nsd;
    scount[i] = (mesh.eDist[i + 1] - mesh.eDist[i]) * nG * nfibers * nsd;
  }
  MPI_Scatterv(reordered.data(), scount.data(), disp.data(), cm_mod::mpreal,
      quadrature_fibers.data(), nEl * nG * nfibers * nsd, cm_mod::mpreal, cm_mod.master, cm.com());
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

  if (num_paths != 0) {
    read_from_file(mesh_parameters, mesh, nsd);
  } else {
    read_directions(mesh_parameters, nsd);
  }

  populate_quadrature_fibers(mesh, nsd);
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
      input_fibers.resize(num_paths * nsd, num_fibers);
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

  input_fibers.resize(num_dirs * nsd, 1);
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
    throw std::runtime_error("Fiber directions cannot be populated before mesh quadrature points are initialized.");
  }

  const int num_components = nfibers * nsd;
  quadrature_fibers.resize(num_components, mesh.nG, mesh.gnEl);

  // Fiber directions are associated with points.
  if (data_location == DataLocation::point) {
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
  // Fiber directions are associated with cells.
  else if (data_location == DataLocation::cell) {
    for (int e = 0; e < mesh.gnEl; ++e) {
      for (int g = 0; g < mesh.nG; ++g) {
        for (int i = 0; i < num_components; ++i) {
          quadrature_fibers(i, g, e) = input_fibers(i, e);
        }
      }
    }
  } 
  // Fiber directions are explicitly defined in the input file
  else if (data_location == DataLocation::explicit_direction) {
    for (int e = 0; e < mesh.gnEl; ++e) {
      for (int g = 0; g < mesh.nG; ++g) {
        for (int i = 0; i < num_components; ++i) {
          quadrature_fibers(i, g, e) = input_fibers(i, 0);
        }
      }
    }
  } else {
    throw std::runtime_error("Fiber input data has no cell, point, or explicit direction association.");
  }
}


Vector<double> Fibers::interpolate_dyadic(
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
