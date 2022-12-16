#ifndef CONTACTF_H 
#define CONTACT_H 

#include "ComMod.h"

namespace contact {

void contact_forces(ComMod& com_mod, CmMod& cm_mod, const Array<double>& Dg);

};

#endif

