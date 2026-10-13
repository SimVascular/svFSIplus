// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIBERS_H
#define FIBERS_H

#include "Array.h"

class ComMod;
class CmMod;
class MeshParameters;
class cmType;
class mshType;

/// @brief Fiber directions associated with one mesh.
///
/// The stored array has the same layout previously used by mshType::fN:
/// (number_of_fibers * spatial_dimension, number_of_elements).
class Fibers {
  public:

    // Reads fiber directions from file paths or constant vectors.
    void read(const MeshParameters& mesh_parameters, mshType& mesh, int nsd);

    // Parallel distribution of fiber directions.
    void distribute(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm,
      mshType& mesh);

    // Clears the fiber directions.
    void clear();

    // Returns the fiber directions.
    const Array<double>& data() const { return fibers; }
    Array<double>& data() { return fibers; }

    // Returns true if there are fiber directions.
    bool has_fibers() const { return fibers.size() != 0; }


  protected:
    Array<double> fibers;

  private:
    // Reads fiber directions from a file.
    void read_from_file(const MeshParameters& mesh_parameters, mshType& mesh, int nsd);

    // Reads constant fiber directions.
    void read_directions(const MeshParameters& mesh_parameters, mshType& mesh, int nsd);
};

#endif
