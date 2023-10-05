#ifndef REMESH_H
#define REMESH_H

#include "Simulation.h"

namespace remesh {

void remesh_restart(Simulation* simulation);

void set_face_ebc(ComMod& com_mod, CmMod& cm_mod, faceType& lFa, mshType& lM);

};  // namespace remesh

#endif
