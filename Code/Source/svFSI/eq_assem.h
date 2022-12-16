
#ifndef EQ_ASSEM_H 
#define EQ_ASSEM_H 

#include "ComMod.h"
#include "Simulation.h"

namespace eq_assem {

void b_assem_neu_bc(ComMod& com_mod, const faceType& lFa, const Vector<double>& hg, const Array<double>& Yg);

void b_neu_folw_p(ComMod& com_mod, const faceType& lFa, const Vector<double>& hg, const Array<double>& Dg);

void global_eq_assem(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg);

};

#endif

