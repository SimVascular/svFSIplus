#ifndef FSI_H 
#define FSI_H 

#include "ComMod.h"

#include "consts.h"

namespace fsi {

void construct_fsi(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg);

};

#endif

