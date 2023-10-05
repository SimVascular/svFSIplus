#ifndef CEP_ION_H
#define CEP_ION_H

#include <string>

#include "Array.h"
#include "ComMod.h"
#include "Simulation.h"
#include "all_fun.h"
#include "consts.h"

namespace cep_ion {

void cep_init(Simulation* simulation);

void cep_init_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG,
                Vector<double>& X, Vector<double>& Xg);

void cep_integ(Simulation* simulation, const int iEq, const int iDof,
               const Array<double>& Dg);

void cep_integ_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG,
                 Vector<double>& X, Vector<double>& Xg, const double t1,
                 double& yl, const double I4f, const double dt);

};  // namespace cep_ion

#endif
