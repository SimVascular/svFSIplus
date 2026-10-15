// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIBERS_H
#define FIBERS_H

#include "Array.h"
#include "Array3.h"

#include "eigen3/Eigen/Core"

class ComMod;
class CmMod;
class MeshParameters;
class cmType;
class mshType;

/// @brief Fiber directions associated with one mesh.
class Fibers {
  public:

    // Number of fiber directions provided as input.
    int nfibers = 0;

    // Reads fiber directions from file paths or constant vectors.
    void read(const MeshParameters& mesh_parameters, mshType& mesh, int nsd);

    // Parallel distribution of fiber directions.
    void distribute(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm,
      mshType& mesh);

    // Clears the fiber directions.
    void clear();

    /// Returns an Array view of the fiber directions for a quadrature point
    /// and element. Explicit directions are the single column of the input
    /// array, cell directions are the requested cell column, and point
    /// directions are their interpolated quadrature-point value. The view
    /// shares this object's storage and must not outlive it.
    Array<double> get_quadrature_fibers(int nsd, int g, int e) const;

    /// Returns a read-only Eigen view of get_quadrature_fibers().
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
    get_eigen_fibers(int nsd, int g, int e) const;

    // Returns true if there are fiber directions.
    bool has_fibers() const { return nfibers != 0; }

  protected:
    // Fiber directions read from the input. The column count depends on the
    // data location: one column for explicit directions, one column per cell
    // for cell data, and one column per point for point data.
    Array<double> input_fibers;

    // Fiber directions at the quadrature points.
    Array3<double> quadrature_fibers;

  private:
    // Location of the fiber directions in the input file.
    enum class DataLocation {
      none,
      // The fiber directions are associated with the mesh nodes.
      point,
      // The fiber directions are associated with the mesh cells.
      cell,
      // The fiber directions are explicit vectors defined in the input file.
      explicit_direction
    };
    DataLocation data_location = DataLocation::none;

    // Reads fiber directions from a file. The fiber directions can be associated with
    // cells or points.
    void read_from_file(const MeshParameters& mesh_parameters, mshType& mesh, int nsd);

    // Reads constant fiber directions as specified in the input file.
    void read_directions(const MeshParameters& mesh_parameters, int nsd);

    // Interpolates point fibers to quadrature points and orthonormalizes them.
    void populate_quadrature_fibers(const mshType& mesh, int nsd);

    // Orthonormalizes multiple fiber directions at every quadrature point and
    // completes them to an nsd-dimensional basis. A single direction is left
    // unchanged.
    void orthonormalize_quadrature_fibers(int nsd);

    // Orthonormalizes explicit or cell-based fibers in place.
    void orthonormalize_fibers(Array<double>& fibers, int nsd);

    // Interpolates nodal fiber directions onto the quadrature points.
    void interpolate_point_fibers(const mshType& mesh, int nsd);

};

#endif
