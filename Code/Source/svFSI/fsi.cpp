/* Copyright (c) Stanford University, The Regents of the University of California, and others.
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

#include "fsi.h"

#include "all_fun.h"
#include "consts.h"
#include "fluid.h"
#include "fs.h"
#include "lhsa.h"
#include "nn.h"
#include "sv_struct.h"
#include "utils.h"

#include <array>
#include <iomanip>
#include <math.h>

#ifdef WITH_TRILINOS
#include "trilinos_linear_solver.h"
#endif

namespace fsi {

void construct_fsi(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg)
{
  #define n_debug_construct_fsi 
  #ifdef debug_construct_fsi
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  int num_alloc = Array<double>::num_allocated;
  com_mod.timer.set_time();
  double elapsed_time = 0.0; 
  #endif

  using namespace consts;

  int eNoN = lM.eNoN;
  int nFn  = lM.nFn;
  if (nFn == 0) nFn = 1;

  bool  vmsStab = false;
  if (lM.nFs == 1) {
     vmsStab = true;
  }

  // l = 3, if nsd==2 ; else 6;
  auto& cem = cep_mod.cem;
  const int l = com_mod.nsymd;
  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;
  const int nsymd = com_mod.nsymd;
  auto& pS0 = com_mod.pS0;
  #ifdef debug_construct_fsi
  dmsg << "lM.nEl: " << lM.nEl << nsd;
  dmsg << "nsd: " << nsd;
  dmsg << "eNoN: " << eNoN;
  dmsg << "vmsStab: " << vmsStab;
  dmsg << "pS0.size(): " << pS0.size();
  #endif

  Vector<int> ptr(eNoN); 
  Array3<double> lK(dof*dof,eNoN,eNoN), lKd(dof*nsd,eNoN,eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), bfl(nsd,eNoN), 
      fN(nsd,nFn), pS0l(nsymd,eNoN), lR(dof,eNoN);
  Vector<double> pSl(nsymd), ya_l(eNoN);

  std::array<fsType,2> fs_1;
  fs::get_thood_fs(com_mod, fs_1, lM, vmsStab, 1);

  std::array<fsType,2> fs_2;
  fs::get_thood_fs(com_mod, fs_2, lM, vmsStab, 2);

  // Loop over all elements of mesh
  //
  double struct_3d_time = 0.0;
  double fluid_3d_time = 0.0;

  for (int e = 0; e < lM.nEl; e++) {
    // setting globals
    cDmn = all_fun::domain(com_mod, lM, cEq, e); 
    auto cPhys = eq.dmn[cDmn].phys;

    if ((cPhys != Equation_fluid) && (cPhys != Equation_lElas)  && 
        (cPhys != Equation_struct) && (cPhys != Equation_ustruct)) {
      continue;
    }

    // Update shape functions for NURBS
    //if (lM.eType == eType_NRB) CALL NRBNNX(lM, e)

    // Create local copies
    fN  = 0.0;
    pS0l = 0.0;
    ya_l = 0.0;

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < xl.nrows(); i++) {
        xl(i,a) = com_mod.x(i,Ac);
        bfl(i,a) = com_mod.Bf(i,Ac);
     }
      for (int i = 0; i < al.nrows(); i++) {
        al(i,a) = Ag(i,Ac);
        yl(i,a) = Yg(i,Ac);
        dl(i,a) = Dg(i,Ac);
      }

      if (lM.fN.size() != 0) {
        for (int iFn = 0; iFn < nFn; iFn++) {
          for (int i = 0; i < nsd; i++) {
            fN(i,iFn) = lM.fN(i+nsd*iFn,e);
          }
        }
      }

      if (pS0.size() != 0) {
        pS0l.set_col(a, pS0.col(Ac));
      }

      if (cem.cpld) {
        ya_l(a) = cem.Ya(Ac);
      }
    }

    // For FSI, fluid domain should be in the current configuration
    //
    if (cPhys == Equation_fluid) {
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < eNoN; j++) {
          xl(i,j) = xl(i,j) + dl(nsd+i+1,j);
        }
      }
    }

    // Initialize residual and tangents
    lR = 0.0;
    lK = 0.0;
    lKd = 0.0;

    //  Define element coordinates appropriate for function spaces
    Array<double> xwl(nsd,fs_1[0].eNoN);
    Array<double> Nwx(nsd,fs_1[0].eNoN);
    Array<double> Nwxx(l,fs_1[0].eNoN);
    Array<double> xql(nsd,fs_1[1].eNoN);
    Array<double> Nqx(nsd,fs_1[1].eNoN);

    xwl = xl;

    for (int i = 0; i < xql.nrows(); i++) {
      for (int j = 0; j < fs_1[1].eNoN; j++) {
        xql(i,j) = xl(i,j);
      }
    }

    // Gauss integration 1
    //
    double Jac{0.0};
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < fs_1[0].nG; g++) {
      if (g == 0 || !fs_1[1].lShpF) {
        auto Nx = fs_1[1].Nx.rslice(g);
        nn::gnn(fs_1[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fsi] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs_1[0].lShpF) {
        auto Nx = fs_1[0].Nx.rslice(g);
        nn::gnn(fs_1[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fsi] Jacobian for element " + std::to_string(e) + " is < 0.");
        }

        if (!vmsStab) {
          auto Nx = fs_1[0].Nx.rslice(g);
          auto Nxx = fs_1[0].Nxx.rslice(g);
          nn::gn_nxx(l, fs_1[0].eNoN, nsd, nsd, Nx, Nxx, xwl, Nwx, Nwxx);
        }
      }

      double w = fs_1[0].w(g) * Jac;

      if (nsd == 3) {
        switch (cPhys) {
          case Equation_fluid: {
            auto N0 = fs_1[0].N.col(g);
            auto N1 = fs_1[1].N.col(g);
            fluid::fluid_3d_m(com_mod, vmsStab, fs_1[0].eNoN, fs_1[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);
          } break;

          case Equation_struct: {
            auto N0 = fs_1[0].N.col(g);
            struct_ns::struct_3d(com_mod, cep_mod, fs_1[0].eNoN, nFn, w, N0, Nwx, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK);
          } break;
          case Equation_lElas:
            throw std::runtime_error("[construct_fsi] LELAS3D not implemented");
            //CALL LELAS3D(fs(1).eNoN, w, fs(1).N(:,g), Nwx, al, dl, bfl, pS0l, pSl, lR, lK)
          break;

          case Equation_ustruct:
            throw std::runtime_error("[construct_fsi] USTRUCT3D_M not implemented");
            //CALL USTRUCT3D_M(vmsStab, fs(1).eNoN, fs(2).eNoN, nFn, w, Jac, fs(1).N(:,g), fs(2).N(:,g), Nwx, al, yl, 
            //                 dl, bfl, fN, ya_l, lR, lK, lKd)
          break;
          }

      } else if (nsd == 2) {
        switch (cPhys) {
          case Equation_fluid: {
            auto N0 = fs_1[0].N.col(g);
            auto N1 = fs_1[1].N.col(g);
            fluid::fluid_2d_m(com_mod, vmsStab, fs_1[0].eNoN, fs_1[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);
          } break;

          case Equation_lElas:
            throw std::runtime_error("[construct_fsi] LELAS2D not implemented");
            //CALL LELAS2D(fs(1).eNoN, w, fs(1).N(:,g), Nwx, al, dl, bfl, pS0l, pSl, lR, lK)
          break;

          case Equation_struct: {
            auto N0 = fs_1[0].N.col(g);
            struct_ns::struct_2d(com_mod, cep_mod, fs_1[0].eNoN, nFn, w, N0, Nwx, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK);
          } break;

          case Equation_ustruct:
            throw std::runtime_error("[construct_fsi] USTRUCT2D_M not implemented");
            //CALL USTRUCT2D_M(vmsStab, fs(1).eNoN, fs(2).eNoN, nFn, w, Jac, fs(1).N(:,g), fs(2).N(:,g), Nwx, al, yl, dl, bfl, fN, ya_l, lR, lK, lKd)
          break;
        }
      }
    } // g: loop

    // Gauss integration 2
    //
    for (int g = 0; g < fs_2[1].nG; g++) {
      if (g == 0 || !fs_2[0].lShpF) {
        auto Nx = fs_2[0].Nx.rslice(g);
        nn::gnn(fs_2[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fsi] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs_2[1].lShpF) {
        auto Nx = fs_2[1].Nx.rslice(g);
        nn::gnn(fs_2[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fsi] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }
      double w = fs_2[1].w(g) * Jac;

      if (nsd == 3) {
        switch (cPhys) {
          case Equation_fluid: {
            auto N0 = fs_2[0].N.col(g);
            auto N1 = fs_2[1].N.col(g);
            fluid::fluid_3d_c(com_mod, vmsStab, fs_2[0].eNoN, fs_2[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);
          } break;

          case Equation_ustruct:
            throw std::runtime_error("[construct_fsi] USTRUCT3D_C not implemented");
            //CALL USTRUCT3D_C(vmsStab, fs(1).eNoN, fs(2).eNoN, w, Jac, fs(1).N(:,g), fs(2).N(:,g), Nwx, Nqx, al, yl, dl, bfl, lR, lK, lKd)
          break;
        }

      } else if (nsd == 2) {
        switch (cPhys) {
          case Equation_fluid: {
            auto N0 = fs_2[0].N.col(g);
            auto N1 = fs_2[1].N.col(g);
            fluid::fluid_2d_c(com_mod, vmsStab, fs_2[0].eNoN, fs_2[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);
          } break;

          case Equation_ustruct:
            throw std::runtime_error("[construct_fsi] USTRUCT2D_C not implemented");
            //CALL USTRUCT2D_C(vmsStab, fs(1).eNoN, fs(2).eNoN, w, Jac, fs(1).N(:,g), fs(2).N(:,g), Nwx, Nqx, al, yl, dl, bfl, lR, lK, lKd)
          break;
        }
      }
    } // g: loop

    // Assembly
#ifdef WITH_TRILINOS
    if (eq.assmTLS) {
      if (cPhys == Equation_ustruct) {
        throw std::runtime_error("[construct_fsi] Cannot assemble USTRUCT using Trilinos");
      }
      trilinos_doassem_(const_cast<int&>(eNoN), ptr.data(), lK.data(), lR.data());
    } else {
#endif
      if (cPhys == Equation_ustruct) {
        //CALL USTRUCT_DOASSEM(eNoN, ptr, lKd, lK, lR)
        throw std::runtime_error("[construct_fsi] USTRUCT_DOASSEM not implemented");
      } else {
        lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
      }
#ifdef WITH_TRILINOS
    }
#endif
  } // e: loop

  #ifdef debug_construct_fsi
  elapsed_time = com_mod.timer.get_elapsed_time();
  dmsg << "elapsed_time: " << elapsed_time;
  dmsg << "num_allocated: " << Array<double>::num_allocated - num_alloc;
  #endif
}

};
