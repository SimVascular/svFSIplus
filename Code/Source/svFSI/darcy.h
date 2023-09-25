//
// This code implements the Darcy Equation for 2D and 3D
// problems in perfusion of porus media.
//

#ifndef DARCY_H
#define DARCY_H

#include "ComMod.h"

namespace darcy {

    void b_darcy(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, Array<double>& lR);

    void construct_darcy(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Dg);

    void darcy_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
        const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK, Vector<double>& local_source,
        Vector<double>& local_sink, Vector<double>& local_b0, Vector<double>& local_b1);

    void darcy_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
        const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK);
}

#endif //DARCY_H
