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

    /// Returns an Array view of the complete nsd x nsd orthonormal fiber basis
    /// for a quadrature point and element. Columns are fiber,
    /// sheet/transverse, and (in 3D) sheet normal. Explicit directions use the
    /// single stored basis, cell directions use the requested cell basis, and
    /// point directions use their interpolated quadrature-point basis. The
    /// view shares this object's storage and must not outlive it.
    Array<double> get_quadrature_fibers(int nsd, int g, int e) const;

    /// Returns a read-only Eigen view of the complete basis returned by
    /// get_quadrature_fibers().
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
    get_eigen_fibers(int nsd, int g, int e) const;

    // Returns true if there are fiber directions.
    bool has_fibers() const { return nfibers != 0; }

  protected:
    // Fiber directions read from the input. Each column stores an nsd x nsd
    // basis (nsd*nsd components). Column count depends on the data location:
    // one for explicit directions, one per cell for cell data, and one per
    // point for point data.
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

    // Reads fiber directions from VTK files. expected_location must be point
    // or cell and selects which FIB_DIR association is required.
    void read_from_file(const MeshParameters& mesh_parameters, mshType& mesh,
        int nsd, DataLocation expected_location);

    // Reads constant fiber directions as specified in the input file.
    void read_directions(const MeshParameters& mesh_parameters, int nsd);

    // Checks that the input fiber directions form an orthonormal set at every
    // stored location. Raises if any direction is not unit length or if any
    // pair is not orthogonal.
    void check_input_fibers_orthonormal(int nsd) const;

    // Interpolates point fibers to quadrature points and orthonormalizes them.
    void populate_quadrature_fibers(const mshType& mesh, int nsd);

    // Orthonormalizes input fiber directions at every quadrature point and
    // completes them to an nsd-dimensional basis.
    void orthonormalize_quadrature_fibers(int nsd);

    // Completes already-orthonormal explicit or cell-based fibers to an
    // nsd-dimensional basis via cross product (3D) or 90-degree rotation (2D).
    void complete_fiber_basis(Array<double>& fibers, int nsd);

    // Interpolates nodal fiber directions onto the quadrature points.
    void interpolate_point_fibers(const mshType& mesh, int nsd);

};

#endif
