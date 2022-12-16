
#include "Simulation.h"

#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

void distribute(Simulation* simulation);

void dist_bc(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bcType& lBc, const std::vector<mshType>& tMs,
             const Vector<int>& gmtl);

void dist_bf(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bfType& lBf);

void dist_eq(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const std::vector<mshType>& tMs,
             const Vector<int>& gmtl, CepMod& cep_mod, eqType& lEq);

void dist_mat_consts(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, stModelType& lStM);

void dist_visc_model(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, viscModelType& lVis);

void part_face(Simulation* simulation, mshType& lM, faceType& lFa, faceType& gFa, Vector<int>& gmtl);

void part_msh(Simulation* simulation, mshType& lM, Vector<int>& mtl, int nP, Vector<float>& wgt);

#endif

