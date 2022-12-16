
#ifndef SHELLS_H 
#define SHELLS_H 

#include "ComMod.h"
#include "Simulation.h"

namespace shells {

void shell_fp(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& dl, const Array<double>& xl, const Vector<double>& tfl, Array<double>& lR, Array3<double>& lK);

};

#endif

