/**
 * Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "eq_assem.h"

#include "all_fun.h"
#include "consts.h"
#include "lhsa.h"
#include "nn.h"
#include "utils.h"

#include "cep.h"
#include "cmm.h"
#include "fluid.h"
#include "fsi.h"
#include "heatf.h"
#include "heats.h"
#include "l_elas.h"
#include "mesh.h"
#include "shells.h"
#include "stokes.h"
#include "sv_struct.h"
#include "ustruct.h"
#include "darcy.h"

#include <math.h>

#ifdef WITH_TRILINOS
#include "trilinos_linear_solver.h"
#endif

namespace eq_assem {

void b_assem_neu_bc(ComMod& com_mod, const faceType& lFa, const Vector<double>& hg, const Array<double>& Yg) 
{
  #define n_debug_b_assem_neu_bc
  #ifdef debug_b_assem_neu_bc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;

  const int iM = lFa.iM;
  const int eNoN = lFa.eNoN;
  const auto& msh = com_mod.msh[iM];
  auto& cDmn = com_mod.cDmn;

  for (int e = 0; e < lFa.nEl; e++) {
    int Ec = lFa.gE(e);
    cDmn = all_fun::domain(com_mod, msh, cEq, Ec);
    auto cPhys = eq.dmn[cDmn].phys;

    Vector<int> ptr(eNoN); 
    Vector<double> N(eNoN), hl(eNoN); 
    Array<double> yl(tDof,eNoN), lR(dof,eNoN); 
    Array3<double> lK(dof*dof,eNoN,eNoN);

    for (int a = 0; a < eNoN; a++) {
      int Ac = lFa.IEN(a,e);
      ptr(a) = Ac;
      hl(a) = hg(Ac);
      for (int i = 0; i < tDof; i++) {
        yl(i,a) = Yg(i,Ac);
      }
    }

    // Updating the shape functions, if neccessary
    if (lFa.eType == ElementType::NRB) {
      //CALL NRBNNXB(msh(iM), lFa, e)
    }

    for (int g = 0; g < lFa.nG; g++) {
      Vector<double> nV(nsd);
      auto Nx = lFa.Nx.slice(g);
      nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, eNoN, Nx, nV);
      double Jac = sqrt(utils::norm(nV));
      nV = nV / Jac;
      double w = lFa.w(g)*Jac;
      N  = lFa.N.col(g);

      double h = 0.0;
      Vector<double> y(tDof);

      for (int a = 0; a < eNoN; a++) {
        h = h + N(a)*hl(a);
        y = y + N(a)*yl.col(a);
      }

      switch ( cPhys) {
        case EquationType::phys_fluid:
          fluid::b_fluid(com_mod, eNoN, w, N, y, h, nV, lR, lK);
        break;

        case EquationType::phys_CMM:
          fluid::b_fluid(com_mod, eNoN, w, N, y, h, nV, lR, lK);
        break;

        case EquationType::phys_heatS:
          heats::b_heats(com_mod, eNoN, w, N, h, lR);
        break;

        case EquationType::phys_darcy:
            darcy::b_darcy(com_mod, eNoN, w, N, h, lR);
            break;

        case EquationType::phys_heatF:
          heatf::b_heatf(com_mod, eNoN, w, N, y, h, nV, lR, lK);
        break;

        case EquationType::phys_lElas:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_struct:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_ustruct:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_shell:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_mesh:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_stokes:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_CEP:
          cep::b_cep(com_mod, eNoN, w, N, h, lR);
        break;

        default:
          throw std::runtime_error("[b_assem_neu_bc] Undefined physics selection for assembly");
      }
    }

    // Now doing the assembly part

#ifdef WITH_TRILINOS
    if (eq.assmTLS) {
      trilinos_doassem_(const_cast<int&>(eNoN), ptr.data(), lK.data(), lR.data());
    } else {
      lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
    }
#else
    lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
#endif
  }
}


/// @brief For struct/ustruct - construct follower pressure load.
///
/// We use Nanson's formula to take change in normal direction with
/// deformation into account. Additional calculations based on mesh
/// need to be performed.
///
/// Reproduces 'SUBROUTINE BNEUFOLWP(lFa, hg, Dg)'
//
void b_neu_folw_p(ComMod& com_mod, const faceType& lFa, const Vector<double>& hg, const Array<double>& Dg) 
{
  using namespace consts;

  #define n_debug_b_neu_folw_p
  #ifdef debug_b_neu_folw_p 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lFa.name: " << lFa.name;
  #endif

  const int cEq = com_mod.cEq;
  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;

  const int iM = lFa.iM;
  const auto& msh = com_mod.msh[iM];
  const int eNoN = msh.eNoN;
  const int eNoNb = lFa.eNoN;

  auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  #ifdef debug_b_neu_folw_p 
  dmsg << "nsd: " << nsd;
  dmsg << "eNoN: " << nsd;
  #endif

  for (int e = 0; e < lFa.nEl; e++) {
    int Ec = lFa.gE(e);
    cDmn = all_fun::domain(com_mod, msh, cEq, Ec);  // Changes global
    auto cPhys = eq.dmn[cDmn].phys;

    Vector<int> ptr(eNoN); 
    Vector<double> hl(eNoN); 
    Array<double> xl(nsd,eNoN); 
    Array<double> dl(tDof,eNoN);
    Vector<double> N(eNoN); 
    Array<double> Nxi(nsd,eNoN); 
    Array<double> Nx(nsd,eNoN); 
    Array<double> lR(dof,eNoN);
    Array3<double> lK(dof*dof,eNoN,eNoN);
    Array3<double> lKd;

    if (cPhys == EquationType::phys_ustruct) {
      lKd.resize(dof*nsd,eNoN,eNoN);
    }

    // Create local copies
    for (int a = 0; a < eNoN; a++) {
      int Ac = msh.IEN(a,Ec);
      ptr(a) = Ac;
      hl(a) = hg(Ac);

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        dl(i,a) = Dg(i,Ac);
      }
    }

    // Initialize parameteric coordinate for Newton's iterations
    Vector<double> xi0(nsd);
    for (int g = 0; g < msh.nG; g++) {
      xi0 = xi0 + msh.xi.col(g);
    }
    xi0 = xi0 / static_cast<double>(msh.nG);

    for (int g = 0; g < lFa.nG; g++) {
      Vector<double> xp(nsd);
      double Jac;

      for (int a = 0; a < eNoNb; a++) {
        int Ac = lFa.IEN(a,e);
        xp = xp + com_mod.x.col(Ac) * lFa.N(a,g);
      }

      auto xi = xi0;
      nn::get_nnx(nsd, msh.eType, eNoN, xl, msh.xib, msh.Nb, xp, xi, N, Nxi);

      if (g == 0 || !msh.lShpF) {
        Array<double> ksix(nsd,nsd);
        nn::gnn(eNoN, nsd, nsd, Nxi, xl, Nx, Jac, ksix);
      }

      Vector<double> nV(nsd);
      auto Nx_g = lFa.Nx.slice(g);
      nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, eNoNb, Nx_g, nV);
      Jac = sqrt(utils::norm(nV));
      nV = nV / Jac;
      double w = lFa.w(g)*Jac;

      if (cPhys == EquationType::phys_ustruct) {
        if (nsd == 3) {
          ustruct::b_ustruct_3d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK, lKd);
        } else {
          ustruct::b_ustruct_2d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK, lKd);
        }

      } else if (cPhys == EquationType::phys_struct) {
        if (nsd == 3) {
          struct_ns::b_struct_3d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK);
        } else {
          struct_ns::b_struct_2d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK);
        }
      }
    }

    // Now doing the assembly part
    //
#ifdef WITH_TRILINOS
    if (eq.assmTLS) {
      trilinos_doassem_(const_cast<int&>(eNoN), const_cast<int*>(ptr.data()), lK.data(), lR.data());
    } else {
#endif
      if (cPhys == EquationType::phys_ustruct) {
        ustruct::ustruct_do_assem(com_mod, eNoN, ptr, lKd, lK, lR);
      } else if (cPhys == EquationType::phys_struct) {
        lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
      }
#ifdef WITH_TRILINOS
    }
#endif
  }
}


/// @brief This routine assembles the equation on a given mesh.
///
/// Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo)
//
void global_eq_assem(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg)
{
  #define n_debug_global_eq_assem
  #ifdef debug_global_eq_assem
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.name: " << lM.name;
  com_mod.timer.set_time();
  #endif

  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  #ifdef debug_global_eq_assem
  dmsg << "cEq: " << cEq;
  dmsg << "eq.sym: " << eq.sym;
  dmsg << "eq.phys: " << eq.phys;
  #endif

  switch (eq.phys) {

    case EquationType::phys_fluid:
      fluid::construct_fluid(com_mod, lM, Ag, Yg);
    break;

    case EquationType::phys_heatF:
      heatf::construct_heatf(com_mod, lM, Ag, Yg);
    break;

    case EquationType::phys_heatS:
      heats::construct_heats(com_mod, lM, Ag, Yg);
    break;

    case EquationType::phys_darcy:
        darcy::construct_darcy(com_mod, lM, Ag, Yg);
        break;

    case EquationType::phys_lElas:
      l_elas::construct_l_elas(com_mod, lM, Ag, Dg);
    break;

    case EquationType::phys_struct:
      struct_ns::construct_dsolid(com_mod, cep_mod, lM, Ag, Yg, Dg);
    break;

    case EquationType::phys_ustruct:
      ustruct::construct_usolid(com_mod, cep_mod, lM, Ag, Yg, Dg);
    break;

    case EquationType::phys_CMM:
      cmm::construct_cmm(com_mod, lM, Ag, Yg, Dg);
    break;

    case EquationType::phys_shell:
      shells::construct_shell(com_mod, lM, Ag, Yg, Dg);
    break;

    case EquationType::phys_FSI:
      fsi::construct_fsi(com_mod, cep_mod, lM, Ag, Yg, Dg);
    break;

    case EquationType::phys_mesh:
      mesh::construct_mesh(com_mod, cep_mod, lM, Ag, Dg);
    break;

    case EquationType::phys_CEP:
      cep::construct_cep(com_mod, cep_mod, lM, Ag, Yg, Dg);
    break;

    case EquationType::phys_stokes:
      stokes::construct_stokes(com_mod, lM, Ag, Yg);
    break;

    default:
      throw std::runtime_error("[global_eq_assem] Undefined physics selection for assembly");
  } 

  #ifdef debug_global_eq_assem
  double elapsed_time = com_mod.timer.get_elapsed_time();
  dmsg << "elapsed_time: " << elapsed_time;
  #endif
}

};


