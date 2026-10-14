// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIBERS_H
#define FIBERS_H

#include "Array.h"
#include "Array3.h"
#include "Vector.h"

class ComMod;
class CmMod;
class MeshParameters;
class cmType;
class mshType;

/// @brief Fiber directions associated with one mesh.
class Fibers {
  public:

    // Number of fiber directions.
    int nfibers = 0;

    // Reads fiber directions from file paths or constant vectors.
    void read(const MeshParameters& mesh_parameters, mshType& mesh, int nsd);

    // Parallel distribution of fiber directions.
    void distribute(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm,
      mshType& mesh);

    // Clears the fiber directions.
    void clear();

    // Returns the fiber directions at the quadrature points. The array layout
    // is (number_of_fibers * spatial_dimension, number_of_quadrature_points,
    // number_of_elements).
    const Array3<double>& data() const { return quadrature_fibers; }
    Array3<double>& data() { return quadrature_fibers; }

    // Returns the normalized fiber directions read from the input.
    const Array<double>& input_data() const { return input_fibers; }

    // Returns true if there are fiber directions.
    bool has_fibers() const { return quadrature_fibers.size() != 0; }

    /// Interpolates one nodal fiber direction to a quadrature point.
    ///
    /// The nodal fibers are the columns of `nodal_fibers`, and
    /// `shape_functions` contains their values at the quadrature point. The
    /// dyads are interpolated before the principal direction is recovered, so
    /// opposite representations of the same fiber direction do not cancel.
    static Vector<double> interpolate_dyadic(
      const Array<double>& nodal_fibers,
      const Vector<double>& shape_functions);


  protected:
    // Input fibers array (either nodal or cell-based).
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

    // Populates the quadrature-point fibers from the input fibers.
    void populate_quadrature_fibers(const mshType& mesh, int nsd);

};

#endif
