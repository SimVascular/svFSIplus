
#ifndef LS_H
#define LS_H

#include "ComMod.h"
#include "Simulation.h"

namespace ls_ns {

void ls_alloc(ComMod& com_mod, eqType& lEq);

void ls_solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL,
              const Vector<double>& res);

};  // namespace ls_ns

#endif
