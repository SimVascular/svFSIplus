#ifndef MESH_H
#define MESH_H

#include "ComMod.h"
#include "consts.h"

namespace mesh {

void construct_mesh(ComMod& com_mod, CepMod& cep_mod, const mshType& lM,
                    const Array<double>& Ag, const Array<double>& Dg);

};

#endif
