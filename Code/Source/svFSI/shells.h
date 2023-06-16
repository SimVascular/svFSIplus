
#ifndef SHELLS_H 
#define SHELLS_H 

#include "ComMod.h"
#include "Simulation.h"

namespace shells {

void construct_shell(ComMod& com_mod, const mshType& lM, const Array<double>& Ag,
    const Array<double>& Yg, const Array<double>& Dg);

void shell_fp(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& dl, const Array<double>& xl, const Vector<double>& tfl, Array<double>& lR, Array3<double>& lK);

void shl_strs_res(const ComMod& com_mod, const dmnType& lDmn, const int nFn, const Array<double>& fNa0,
    const double aa_0[2][2], const double aa_x[2][2], const double bb_0[2][2], const double bb_x[2][2],
    double& lam3, Array<double>& Sm, Array3<double>& Dm);

};

#endif

