#ifndef CONTACTF_H
#define CONTACT_H

#include "ComMod.h"

namespace contact {

void construct_contact_pnlty(ComMod& com_mod, CmMod& cm_mod,
                             const Array<double>& Dg);

};

#endif
