// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Fibers.h"

#include "ComMod.h"
#include "Parameters.h"
#include "utils.h"

#include <mpi.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <cmath>
#include <stdexcept>

void Fibers::clear()
{
  fibers.clear();
}

void Fibers::read(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  const int num_paths = mesh_parameters.fiber_direction_file_paths.size();
  const int num_dirs = mesh_parameters.fiber_directions.size();
  if (num_paths == 0 && num_dirs == 0) {
    return;
  }

  if (num_paths != 0) {
    read_from_file(mesh_parameters, mesh, nsd);
  } else {
    read_directions(mesh_parameters, mesh, nsd);
  }
}

void Fibers::read_from_file(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  const int num_paths = mesh_parameters.fiber_direction_file_paths.size();
  mesh.nFn = num_paths;

  fibers.resize(num_paths * nsd, mesh.gnEl);
  fibers = 0.0;
  const auto paths = mesh_parameters.fiber_direction_file_paths();
  for (int idx = 0; idx < mesh.nFn; ++idx) {
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

    // @todo[javijv4]: should this check the number of nodes of the mesh?
    if (grid->GetNumberOfPoints() == 0) {
      throw std::runtime_error("Failed reading the VTK file '" + file_name + "'.");
    }

    if (mesh.gnEl != grid->GetNumberOfCells()) {
      throw std::runtime_error("The number of elements (" + std::to_string(grid->GetNumberOfCells()) +
          ") in the fiber direction VTK file '" + file_name + "' is not equal to the number of elements (" +
          std::to_string(mesh.gnEl) + ") for the mesh named '" + mesh.name + "'.");
    }

    // Get the 3-component fiber orientation data.
    auto array = vtkDoubleArray::SafeDownCast(grid->GetCellData()->GetArray("FIB_DIR"));
    if (array == nullptr) {
      throw std::runtime_error("No 'FIB_DIR' data found in the fiber direction VTK file '" + file_name + "'");
    }

    // Set the fiber orientations data.
    const int offset = idx * nsd;
    for (int e = 0; e < mesh.gnEl; ++e) {
      const auto direction = array->GetTuple(e);
      for (int i = 0; i < nsd; ++i) {
        fibers(offset + i, e) = direction[i];
      }
    }
  }
}

void Fibers::read_directions(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  const int num_dirs = mesh_parameters.fiber_directions.size();
  mesh.nFn = num_dirs;

  fibers.resize(num_dirs * nsd, mesh.gnEl);
  fibers = 0.0;
  const auto directions = mesh_parameters.fiber_directions;
  for (int i = 0; i < mesh.nFn; ++i) {
    auto direction = directions[i]();
    const double length = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
    if (!utils::is_zero(length)) {
      direction[0] /= length;
      direction[1] /= length;
      direction[2] /= length;
    }
    for (int e = 0; e < mesh.gnEl; ++e) {
      for (int j = 0; j < nsd; ++j) {
        fibers(i * nsd + j, e) = direction[j];
      }
    }
  }
}

void Fibers::distribute(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm,
    mshType& mesh)
{
  const int nsd = com_mod.nsd;
  const int nFn = mesh.nFn;
  const int num_proc = cm.np();
  const int nEl = mesh.nEl;
  bool present = false;
  Array<double> reordered;
  Vector<int> scount(num_proc), disp(num_proc);

  if (cm.mas(cm_mod) && has_fibers()) {
    present = true;
    reordered.resize(nFn * nsd, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      reordered.set_col(mesh.otnIEN[e], fibers.col(e));
    }
    clear();
  }
  cm.bcast(cm_mod, &present);
  if (!present) {
    return;
  }

  fibers.resize(nFn * nsd, nEl);
  for (int i = 0; i < num_proc; ++i) {
    disp[i] = mesh.eDist[i] * nFn * nsd;
    scount[i] = (mesh.eDist[i + 1] - mesh.eDist[i]) * nFn * nsd;
  }
  MPI_Scatterv(reordered.data(), scount.data(), disp.data(), cm_mod::mpreal,
      fibers.data(), nEl * nFn * nsd, cm_mod::mpreal, cm_mod.master, cm.com());
}
