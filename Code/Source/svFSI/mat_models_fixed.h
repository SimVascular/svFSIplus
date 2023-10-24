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

// These are template functions reproducing mat_model functions using
// fixed size arrays. 

#ifndef MAT_MODELS_FIXED_H 
#define MAT_MODELS_FIXED_H 

#include "Array.h"
#include "CepMod.h"
#include "ComMod.h"
#include "Tensor4.h"

#include "mat_fun.h"

namespace mat_models {

template <size_t N>
void cc_to_voigt(const double CC[N][N][N][N], double Dm[2*N][2*N])
{
  if (N == 3) {
    Dm[0][0] = CC[0][0][0][0];
    Dm[0][1] = CC[0][0][1][1];
    Dm[0][2] = CC[0][0][2][2];
    Dm[0][3] = CC[0][0][0][1];
    Dm[0][4] = CC[0][0][1][2];
    Dm[0][5] = CC[0][0][2][0];

    Dm[1][1] = CC[1][1][1][1];
    Dm[1][2] = CC[1][1][2][2];
    Dm[1][3] = CC[1][1][0][1];
    Dm[1][4] = CC[1][1][1][2];
    Dm[1][5] = CC[1][1][2][0];

    Dm[2][2] = CC[2][2][2][2];
    Dm[2][3] = CC[2][2][0][1];
    Dm[2][4] = CC[2][2][1][2];
    Dm[2][5] = CC[2][2][2][0];

    Dm[3][3] = CC[0][1][0][1];
    Dm[3][4] = CC[0][1][1][2];
    Dm[3][5] = CC[0][1][2][0];

    Dm[4][4] = CC[1][2][1][2];
    Dm[4][5] = CC[1][2][2][0];

    Dm[5][5] = CC[2][0][2][0];

    for (int i = 1; i < 6; i++) {
      for (int j = 0; j <= i-1; j++) {
        Dm[i][j] = Dm[j][i];
      }
    }

  } else if (N == 2) { 
     Dm[0][0] = CC[0][0][0][0];
     Dm[0][1] = CC[0][0][1][1];
     Dm[0][2] = CC[0][0][0][1];

     Dm[1][1] = CC[1][1][1][1];
     Dm[1][2] = CC[1][1][0][1];

     Dm[2][2] = CC[0][1][0][1];

     Dm[1][0] = Dm[0][1];
     Dm[2][0] = Dm[0][2];
     Dm[2][1] = Dm[1][2];
  } 
}

template <size_t N>
void get_pk2cc(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const double F[N][N], const int nfd,
    const Array<double>& fl, const double ya, double S[N][N], double Dm[2*N][2*N])
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_get_pk2cc
  int task_id = com_mod.cm.idcm();
  std::string msg_prefix;
  msg_prefix = std::string("[get_pk2cc(fixed):") + std::to_string(task_id) + "] ";
  #ifdef debug_get_pk2cc
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== get_pk2cc fixed ==========" << std::endl;
  #endif

  int nsd = com_mod.nsd;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      S[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      Dm[i][j] = 0.0;
    }
  }

  // Some preliminaries
  const auto& stM = lDmn.stM;
  double nd = static_cast<double>(nsd);
  double Kp = stM.Kpen;
  //std::cout << msg_prefix << "Kp: " << Kp << std::endl;

  // Fiber-reinforced stress
  double Tfa = 0.0;
  get_fib_stress(com_mod, cep_mod, stM.Tf, Tfa);
  //CALL GETIBSTRESS(stM.Tf, Tfa)

  // Electromechanics coupling - active stress
  if (cep_mod.cem.aStress) {
    Tfa = Tfa + ya;
  }
  //std::cout << msg_prefix << "Tfa: " << Tfa << std::endl;

  // Electromechanics coupling - active strain
  double Fe[N][N];
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Fe[i][j]  = F[i][j];
    }
  }

  //Array<double> Fa(nsd,nsd);
  double Fa[N][N]{0.0};
  double Fai[N][N]{0.0};
  //auto Fa = mat_id(nsd);
  //auto Fai = Fa;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) {
        Fa[i][j] = 1.0;
        Fai[i][j] = 1.0;
      }
    }
  }

  //print<N>("[get_pk2cc] F", F);
  //print<N>("[get_pk2cc] Fa", Fa);

  /*
  if (cep_mod.cem.aStrain) {
    actv_strain(com_mod, cep_mod, ya, nfd, fl, Fa);
    //CALL ACTVSTRAIN(ya, nfd, fl, Fa)
    Fai = mat_inv(Fa, nsd);
    Fe = mat_mul(F, Fai);
  }
  */

  double J = mat_det<N>(Fe);
  double J2d = pow(J, (-2.0/nd));
  double J4d = J2d*J2d;

  double Idm[N][N];
  mat_id<N>(Idm);

  double Fe_t[N][N];
  transpose<N>(Fe, Fe_t);
  //print<N>("[get_pk2cc] Fe_t", Fe_t);

  double C[N][N];
  mat_mul<N>(Fe_t, Fe, C);

  double E[N][N];
  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < nsd; j++) {
      E[i][j] = 0.50 * (C[i][j] - Idm[i][j]);
    }
  }

  double Ci[N][N];
  mat_inv<N>(C, Ci);
  //std::cout << msg_prefix << "num_allocated 2: " << Array<double>::num_allocated - num_alloc << std::endl;

  double Cprod[N][N];
  mat_mul<N>(C,C,Cprod);
  double trE = mat_trace<N>(E);
  double Inv1 = J2d * mat_trace<N>(C);
  double Inv2 = 0.50 * (Inv1*Inv1 - J4d * mat_trace<N>(Cprod));

  // Contribution of dilational penalty terms to S and CC
  double p  = 0.0;
  double pl = 0.0;

  if (!utils::is_zero(Kp)) {
    get_svol_p(com_mod, cep_mod, stM, J, p, pl);
    // CALL GETSVOLP(stM, J, p, pl)
  }
  //std::cout << msg_prefix << "p: " << p << std::endl;
  //std::cout << msg_prefix << "pl: " << pl << std::endl;

  // Now, compute isochoric and total stress, elasticity tensors
  //
  double CC[N][N][N][N];
  double Idm_prod[N][N][N][N];
  double Ids[N][N][N][N];

  switch (stM.isoType) {
    case ConstitutiveModelType::stIso_lin: {
      double g1 = stM.C10;    // mu
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] = g1 * Idm[i][j];
        }
      }
      return; 
    } break;

    // St.Venant-Kirchhoff
    case ConstitutiveModelType::stIso_StVK: {
      double g1 = stM.C10;         // lambda
      double g2 = stM.C01 * 2.0;   // 2*mu

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] = g1*trE*Idm[i][j] + g2*E[i][j];
        }
      }

      ten_dyad_prod<N>(Idm, Idm, Idm_prod);
      ten_ids<N>(Ids);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] = g1 * Idm_prod[i][j][k][l] + g2*Ids[i][j][k][l];
            }
          }
        }
      }
      //CC = g1 * ten_dyad_prod(Idm, Idm, nsd) + g2*ten_ids(nsd);
    } break;

    // modified St.Venant-Kirchhoff
    case ConstitutiveModelType::stIso_mStVK: {
      double g1 = stM.C10; // kappa
      double g2 = stM.C01;  // mu

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] = g1*log(J)*Ci[i][j] + g2*(C[i][j]-Idm[i][j]);
        }
      }
      //CC = g1 * ( -2.0*log(J)*ten_symm_prod(Ci, Ci, nsd) + ten_dyad_prod(Ci, Ci, nsd) ) + 2.0*g2*ten_ids(nsd);
    } break;

    // NeoHookean model
    case ConstitutiveModelType::stIso_nHook: {
      //std::cout << msg_prefix << "NeoHookean model  " << std::endl;
      double g1 = 2.0 * stM.C10;

      double Sb[N][N];
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] = g1*Idm[i][j];
        }
      }

      // Fiber reinforcement/active stress
      double prod[N][N];
      mat_dyad_prod<N>(fl.col(0), fl.col(0), prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] += Tfa * prod[i][j];
        }
      }
      //Sb += Tfa * mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      double r1 = g1 * Inv1 / nd;
      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          S[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }
      //S = J2d*Sb - r1*Ci;

      double Ci_S_prod[N][N][N][N];
      double S_Ci_prod[N][N][N][N];
      ten_dyad_prod<N>(Ci, S, Ci_S_prod);
      ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          S[i][j] += p * J * Ci[i][j];
        }
      }
      //S += p*J*Ci;

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] = (-2.0/nsd) * Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l];
            }
          }
        }
      }
      //CC = (-2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd));

      double Ci_sym_prod[N][N][N][N];
      ten_symm_prod<N>(Ci, Ci, Ci_sym_prod);

      double Ci_Ci_prod[N][N][N][N];
      ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] += 2.0*(r1 - p*J) * Ci_sym_prod[i][j][k][l] +  (pl*J - 2.0*r1/nd) * Ci_Ci_prod[i][j][k][l];
            }
          }
        }
      }
      //CC += 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd)  +  (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);

    } break;

    //  Mooney-Rivlin model
    case ConstitutiveModelType::stIso_MR: {
      double g1  = 2.0 * (stM.C10 + Inv1*stM.C01);
      double g2  = -2.0 * stM.C01;

      double Sb[N][N];
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] = g1*Idm[i][j] + g2*J2d*C[i][j];
        }
      }
      //auto Sb = g1*Idm + g2*J2d*C;

      // Fiber reinforcement/active stress
      double prod[N][N];
      mat_dyad_prod<N>(fl.col(0), fl.col(0), prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] += Tfa*prod[i][j];
        }
      }
      //Sb = Sb + Tfa*mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      g1 = 4.0 * J4d * stM.C01;

      double CCb[N][N][N][N];
      ten_dyad_prod<N>(Idm, Idm, CCb);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CCb[i][j][k][l] = g1 * (Idm_prod[i][j][k][l] - Ids[i][j][k][l]);
            }
          }
        }
      }

      double r1 = J2d * mat_ddot<N>(C, Sb) / nd;
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }
      //S = J2d*Sb - r1*Ci;

      double Ci_Ci_prod[N][N][N][N];
      double Ci_C_prod[N][N][N][N];
      ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);
      ten_dyad_prod<N>(Ci, C, Ci_C_prod);

      double PP[N][N][N][N];

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              PP[i][j][k][l] = Ids[i][j][k][l] - (1.0/nd) * Ci_C_prod[i][j][k][l];
            }
          }
        }
      }
      //auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);

      ten_ddot<N>(CCb, PP, CC);
      //CC = ten_ddot(CCb, PP, nsd);

      double CC_t[N][N][N][N];
      ten_transpose<N>(CC, CC_t);
      //CC = ten_transpose(CC, nsd);

      ten_ddot<N>(PP, CC_t, CC);
      //CC = ten_ddot(PP, CC, nsd);

      double Ci_S_prod[N][N][N][N];
      double S_Ci_prod[N][N][N][N];
      ten_dyad_prod<N>(Ci, S, Ci_S_prod);
      ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] -= (2.0/nd) * ( Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l]);
            }
          }
        }
      }
      //CC = CC - (2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd) );

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] += p * J * Ci[i][j];
        }
      }
      //S  = S + p*J*Ci;

      double Ci_sym_prod[N][N][N][N];
      ten_symm_prod<N>(Ci, Ci, Ci_sym_prod);

      //double Ci_Ci_prod[N][N][N][N];
      //ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] += 2.0*(r1 - p*J) * Ci_sym_prod[i][j][k][l]  + (pl*J - 2.0*r1/nd) * Ci_Ci_prod[i][j][k][l];
            }
          }
        }
      }
      //CC = CC + 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd) + (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);

    } break;

    // HGO (Holzapfel-Gasser-Ogden) model with additive splitting of
    // the anisotropic fiber-based strain-energy terms
    case ConstitutiveModelType::stIso_HGO: {
      if (nfd != 2) {
        //err = "Min fiber directions not defined for HGO material model (2)"
      }
      double kap = stM.kap;

      double C_fl[N];
      mat_mul(C, fl.rcol(0), C_fl);
      double Inv4 = J2d * mat_fun::norm<N>(fl.rcol(0), C_fl);

      mat_mul(C, fl.rcol(1), C_fl);
      double Inv6 = J2d * mat_fun::norm<N>(fl.rcol(1), C_fl);

      //double Inv4 = J2d*utils::norm(fl.col(0), mat_mul(C, fl.col(0)));
      //double Inv6 = J2d*utils::norm(fl.col(1), mat_mul(C, fl.col(1)));

      double Eff = kap*Inv1 + (1.0-3.0*kap)*Inv4 - 1.0;
      double Ess = kap*Inv1 + (1.0-3.0*kap)*Inv6 - 1.0;

/*

      auto Hff = mat_dyad_prod(fl.col(0), fl.col(0), nsd);
      Hff = kap*Idm + (1.0-3.0*kap)*Hff;
      auto Hss = mat_dyad_prod(fl.col(1), fl.col(1), nsd);
      Hss = kap*Idm + (1.0-3.0*kap)*Hss;

      double g1 = stM.C10;
      double g2 = stM.aff * Eff * exp(stM.bff*Eff*Eff);
      double g3 = stM.ass * Ess * exp(stM.bss*Ess*Ess);
      auto Sb = 2.0*(g1*Idm + g2*Hff + g3*Hss);

      // Fiber reinforcement/active stress
      Sb = Sb + Tfa*mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      g1 = stM.aff*(1.0 + 2.0*stM.bff*Eff*Eff)*exp(stM.bff*Eff*Eff);
      g2 = stM.ass*(1.0 + 2.0*stM.bss*Ess*Ess)*exp(stM.bss*Ess*Ess);
      g1 = 4.0*J4d*g1;
      g2 = 4.0*J4d*g2;

      auto CCb = g1 * ten_dyad_prod(Hff, Hff, nsd) + g2 * ten_dyad_prod(Hss, Hss, nsd);
      double r1  = J2d*mat_ddot(C, Sb, nsd) / nd;
      S = J2d*Sb - r1*Ci;

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC = ten_ddot(CCb, PP, nsd);
      CC = ten_transpose(CC, nsd);
      CC = ten_ddot(PP, CC, nsd);
      CC = CC - (2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd) );

      S = S + p*J*Ci;
      CC = CC + 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd) + (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);
*/
    } break;

    // Guccione (1995) transversely isotropic model
    case ConstitutiveModelType::stIso_Gucci: {
/*
      #ifdef debug_get_pk2cc
      std::cout << msg_prefix << "stIso_Gucci" << std::endl;
      #endif
      if (nfd != 2) {
        //err = "Min fiber directions not defined for Guccione material model (2)"
      }

      // Compute isochoric component of E
      auto E = 0.50 * (J2d*C - Idm);

      // Transform into local orthogonal coordinate system
      Array<double> Rm(nsd,nsd);

      Rm.set_col(0, fl.col(0));
      Rm.set_col(1, fl.col(1));
      Rm.set_col(2, cross(fl));

      // Project E to local orthogocal coordinate system
      auto Es = mat_mul(E, Rm);
      Es = mat_mul(transpose(Rm), Es);

      double g1 = stM.bff;
      double g2 = stM.bss;
      double g3 = stM.bfs;

      double QQ = g1 *  Es(0,0)*Es(0,0) + 
                  g2 * (Es(1,1)*Es(1,1) + Es(2,2)*Es(2,2) + Es(1,2)*Es(1,2) + Es(2,1)*Es(2,1)) +
                  g3 * (Es(0,1)*Es(0,1) + Es(1,0)*Es(1,0) + Es(0,2)*Es(0,2) + Es(2,0)*Es(2,0));

      double r2 = stM.C10 * exp(QQ);

      #ifdef debug_get_pk2cc
      //std::cout << msg_prefix << "E: " << E << std::endl;
      //std::cout << msg_prefix << "Rm: " << Rm << std::endl;
      //std::cout << msg_prefix << "Es: " << Es << std::endl;
      //std::cout << msg_prefix << "g1: " << g1 << std::endl;
      //std::cout << msg_prefix << "g2: " << g2 << std::endl;
      //std::cout << msg_prefix << "g3: " << g3 << std::endl;
      //std::cout << msg_prefix << "QQ: " << QQ << std::endl;
      //std::cout << msg_prefix << "r2: " << r2<< std::endl;
      #endif

      // Fiber stiffness contribution := (dE*_ab / dE_IJ)
      Array3<double> RmRm(nsd,nsd,6);

      RmRm.set_slice(0, mat_dyad_prod(Rm.col(0), Rm.col(0), nsd));
      RmRm.set_slice(1, mat_dyad_prod(Rm.col(1), Rm.col(1), nsd));
      RmRm.set_slice(2, mat_dyad_prod(Rm.col(2), Rm.col(2), nsd));

      RmRm.set_slice(3, mat_symm_prod(Rm.col(0), Rm.col(1), nsd));
      RmRm.set_slice(4, mat_symm_prod(Rm.col(1), Rm.col(2), nsd));
      RmRm.set_slice(5, mat_symm_prod(Rm.col(2), Rm.col(0), nsd));

      //RmRm(:,:,1) = mat_dyad_prod(Rm(:,1), Rm(:,1), nsd)
      //RmRm(:,:,2) = mat_dyad_prod(Rm(:,2), Rm(:,2), nsd)
      //RmRm(:,:,3) = mat_dyad_prod(Rm(:,3), Rm(:,3), nsd)
      //RmRm(:,:,4) = mat_symm_prod(Rm(:,1), Rm(:,2), nsd)
      //RmRm(:,:,5) = mat_symm_prod(Rm(:,2), Rm(:,3), nsd)
      //RmRm(:,:,6) = mat_symm_prod(Rm(:,3), Rm(:,1), nsd)

      auto Sb = g1 *  Es(0,0) * RmRm.slice(0) + 
                g2 * (Es(1,1) * RmRm.slice(1) + Es(2,2)*RmRm.slice(2) + 2.0*Es(1,2)*RmRm.slice(4)) +
          2.0 * g3 * (Es(0,1) * RmRm.slice(3) + Es(0,2)*RmRm.slice(5));

      //Sb = g1*Es(1,1)*RmRm(:,:,1) + g2*(Es(2,2)*RmRm(:,:,2) +
      //  Es(3,3)*RmRm(:,:,3) + 2.0*Es(2,3)*RmRm(:,:,5)) +
      //  2.0*g3*(Es(1,2)*RmRm(:,:,4) + Es(1,3)*RmRm(:,:,6))

      auto CCb = 2.0*ten_dyad_prod(Sb, Sb, nsd);
      Sb = Sb * r2;
      #ifdef debug_get_pk2cc
      //std::cout << msg_prefix << "Sb: " << Sb << std::endl;
      //std::cout << msg_prefix << "CCb: " << CCb << std::endl;
      //std::cout << msg_prefix << "fl: " << fl << std::endl;
      #endif

      // Fiber reinforcement/active stress
      Sb = Sb + Tfa*mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      double r1 = J2d*mat_ddot(C, Sb, nsd) / nd;
      S = J2d*Sb - r1*Ci;
      //std::cout << msg_prefix << "Sb: " << Sb << std::endl;
      //std::cout << msg_prefix << "r1: " << r1 << std::endl;
      //std::cout << msg_prefix << "S: " << S << std::endl;

      r2  = r2*J4d;
      //std::cout << msg_prefix << "r2: " << r2 << std::endl;

      CCb = r2*(CCb + g1 * ten_dyad_prod(RmRm.slice(0), RmRm.slice(0), nsd) + 
                      g2 * (ten_dyad_prod(RmRm.slice(1), RmRm.slice(1), nsd) +
                           ten_dyad_prod(RmRm.slice(2), RmRm.slice(2), nsd) +
                           ten_dyad_prod(RmRm.slice(4), RmRm.slice(4), nsd)*2.0) +
                2.0 * g3 * (ten_dyad_prod(RmRm.slice(3), RmRm.slice(3), nsd) +
                ten_dyad_prod(RmRm.slice(5), RmRm.slice(5), nsd)));


#if 0

      for (int i = 0; i < nsd; i++) { 
        for (int j = 0; j < nsd; j++) { 
          for (int k = 0; k < nsd; k++) { 
            for (int l = 0; l < nsd; l++) { 
              std::cout << msg_prefix << "CCb: " << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << " " << CCb(i,j,k,l) << std::endl;
            }
          }
        }
      }
#endif

      #ifdef debug_get_pk2cc
      //std::cout << msg_prefix << "r2: " << r2 << std::endl;
      //std::cout << msg_prefix << "RmRm: " << RmRm << std::endl;
      //std::cout << msg_prefix << "Sb: " << Sb << std::endl;
      //std::cout << msg_prefix << "CCb: " << CCb << std::endl;
      #endif
      

      // CCb = r2*(CCb + g1*ten_dyad_prod(RmRm(:,:,1), RmRm(:,:,1), nsd)
      // + g2*(ten_dyad_prod(RmRm(:,:,2), RmRm(:,:,2), nsd)
      // + ten_dyad_prod(RmRm(:,:,3), RmRm(:,:,3), nsd)
      // + ten_dyad_prod(RmRm(:,:,5), RmRm(:,:,5), nsd)*2.0)
      // + 2.0*g3*(ten_dyad_prod(RmRm(:,:,4), RmRm(:,:,4), nsd)
      // + ten_dyad_prod(RmRm(:,:,6), RmRm(:,:,6), nsd)))

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC = ten_ddot(CCb, PP, nsd);
      CC  = ten_transpose(CC, nsd);
      CC  = ten_ddot(PP, CC, nsd);
      CC  = CC - (2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd) );

      S  = S + p*J*Ci;
      CC = CC + 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd) + (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);
*/
    } break;

    //  HO (Holzapfel-Ogden) model for myocardium (2009)
    case ConstitutiveModelType::stIso_HO: {
      if (nfd != 2) {
        //err = "Min fiber directions not defined for Holzapfel material model (2)"
      }
/*
      double Inv4 = J2d*utils::norm(fl.col(0), mat_mul(C, fl.col(0)));
      double Inv6 = J2d*utils::norm(fl.col(1), mat_mul(C, fl.col(1)));
      double Inv8 = J2d*utils::norm(fl.col(0), mat_mul(C, fl.col(1)));

      double Eff = Inv4 - 1.0;
      double Ess = Inv6 - 1.0;
      double Efs = Inv8;

      double g1 = stM.a * exp(stM.b*(Inv1-3.0));
      double g2 = 2.0 * stM.afs * Efs * exp(stM.bfs*Efs*Efs);
      auto Hfs = mat_symm_prod(fl.col(0), fl.col(1), nsd);
      auto Sb = g1*Idm + g2*Hfs;

      Efs = Efs * Efs;
      g1 = 2.0*J4d*stM.b*g1;
      g2 = 4.0*J4d*stM.afs*(1.0 + 2.0*stM.bfs*Efs)* exp(stM.bfs*Efs);
      auto CCb  = g1 * ten_dyad_prod(Idm, Idm, nsd) + g2 * ten_dyad_prod(Hfs, Hfs, nsd);

      //  Fiber reinforcement/active stress
      if (Eff > 0.0) {
        g1 = Tfa;

        g1 = g1 + 2.0 * stM.aff * Eff * exp(stM.bff*Eff*Eff);
        auto Hff = mat_dyad_prod(fl.col(0), fl.col(0), nsd);
        Sb  = Sb + g1*Hff;

        Eff = Eff * Eff;
        g1  = 4.0*J4d*stM.aff*(1.0 + 2.0*stM.bff*Eff)*exp(stM.bff*Eff);
        CCb = CCb + g1*ten_dyad_prod(Hff, Hff, nsd);
      }

      if (Ess > 0.0) {
        g2 = 2.0 * stM.ass * Ess * exp(stM.bss*Ess*Ess);
        auto Hss = mat_dyad_prod(fl.col(1), fl.col(1), nsd);
        Sb  = Sb + g2*Hss;

        Ess = Ess * Ess;
        g2  = 4.0*J4d*stM.ass*(1.0 + 2.0*stM.bss*Ess)*exp(stM.bss*Ess);
        CCb = CCb + g2*ten_dyad_prod(Hss, Hss, nsd);
      }

      double r1 = J2d*mat_ddot(C, Sb, nsd) / nd;
      S  = J2d*Sb - r1*Ci;

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC = ten_ddot(CCb, PP, nsd);
      CC  = ten_transpose(CC, nsd);
      CC  = ten_ddot(PP, CC, nsd);
      CC  = CC - (2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd) );

      S   = S + p*J*Ci;
      CC  = CC + 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd) + (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);

      if (cep_mod.cem.aStrain) {
        S = mat_mul(Fai, S);
        S = mat_mul(S, transpose(Fai));
        CCb = 0.0;
        CCb = ten_dyad_prod(Fai, Fai, nsd);
        CC = ten_ddot_3424(CC, CCb, nsd);
        CC = ten_ddot_2412(CCb, CC, nsd);
      }
*/
    } break;

    default:
      throw std::runtime_error("Undefined material constitutive model.");
  } 

  // Convert to Voigt Notation
  cc_to_voigt<N>(CC, Dm);
  //cc_to_voigt(nsd, CC, Dm);
}




};

#endif

