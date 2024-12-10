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

// These are template functions reproducing mat_model functions using C++ arrays. 

#ifndef MAT_MODELS_CARRAY_H 
#define MAT_MODELS_CARRAY_H 

#include "Array.h"
#include "CepMod.h"
#include "ComMod.h"
#include "Tensor4.h"

#include "mat_fun.h"
#include "mat_fun_carray.h"

namespace mat_models_carray {

//--------------
// cc_to_voigt
//--------------
//
template <size_t N>
void cc_to_voigt_carray(const double CC[N][N][N][N], double Dm[2*N][2*N])
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
void voigt_to_cc_carray(const double Dm[2*N][2*N], double CC[N][N][N][N]) {
    if (N == 3) {
        // Initialize the CC array with zeros
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                for (size_t k = 0; k < N; ++k) {
                    for (size_t l = 0; l < N; ++l) {
                        CC[i][j][k][l] = 0.0;
                    }
                }
            }
        }

        // Voigt indices mapping
        const int index_map[6][2] = {
            {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {2, 0}
        };

        // Fill in the CC array based on the Dm matrix
        for (int I = 0; I < 6; ++I) {
            for (int J = 0; J < 6; ++J) {
                int i = index_map[I][0], j = index_map[I][1];
                int k = index_map[J][0], l = index_map[J][1];
                CC[i][j][k][l] = Dm[I][J];
                CC[j][i][k][l] = Dm[I][J];
                CC[i][j][l][k] = Dm[I][J];
                CC[j][i][l][k] = Dm[I][J];
            }
        }
    } else if (N == 2) {
        // Initialize the CC array with zeros
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                for (size_t k = 0; k < N; ++k) {
                    for (size_t l = 0; l < N; ++l) {
                        CC[i][j][k][l] = 0.0;
                    }
                }
            }
        }

        // Voigt indices mapping for 2D
        const int index_map[3][2] = {
            {0, 0}, {1, 1}, {0, 1}
        };

        // Fill in the CC array based on the Dm matrix
        for (int I = 0; I < 3; ++I) {
            for (int J = 0; J < 3; ++J) {
                int i = index_map[I][0], j = index_map[I][1];
                int k = index_map[J][0], l = index_map[J][1];
                CC[i][j][k][l] = Dm[I][J];
                CC[i][j][l][k] = Dm[I][J];
                CC[j][i][k][l] = Dm[I][J];
                CC[j][i][l][k] = Dm[I][J];
            }
        }
    }
}

/**
 * @brief Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
 * including both dilational and isochoric components.
 *
 * Reproduces the Fortran 'GETPK2CC' subroutine.
 *
 * @tparam N Template parameter for number of spatial dimensions.
 * @param[in] com_mod Object containing global common variables.
 * @param[in] cep_mod Object containing electrophysiology-specific common variables.
 * @param[in] lDmn Domain object.
 * @param[in] F Deformation gradient tensor.
 * @param[in] nfd Number of fiber directions.
 * @param[in] fl Fiber directions.
 * @param[in] ya Electrophysiology active stress.
 * @param[out] S 2nd Piola-Kirchhoff stress tensor (modified in place).
 * @param[out] Dm Material stiffness tensor (modified in place).
 * @return None, but modifies S and Dm in place.
 */
template <size_t N>
void get_pk2cc(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const double F[N][N], const int nfd,
    const Array<double>& fl, const double ya, double S[N][N], double Dm[2*N][2*N])
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  using CArray2 = double[N][N];
  using CArray4 = double[N][N][N][N];
  int nsd = com_mod.nsd;

  #define n_debug_get_pk2cc
  int task_id = com_mod.cm.idcm();
  std::string msg_prefix;
  msg_prefix = std::string("[get_pk2cc(carray):") + std::to_string(task_id) + "] ";
  #ifdef debug_get_pk2cc
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== get_pk2cc carray ==========" << std::endl;
  std::cout << msg_prefix << "N: " << N << std::endl;
  std::cout << msg_prefix << "nsd: " << nsd << std::endl;
  #endif

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

  // Fiber-reinforced stress
  double Tfa = 0.0;
  mat_models::get_fib_stress(com_mod, cep_mod, stM.Tf, Tfa);
  double Tsa = Tfa*stM.Tf.eta_s;

  // Electromechanics coupling - active stress
  if (cep_mod.cem.aStress) {
    Tfa = Tfa + ya;
  }

  // Electromechanics coupling - active strain
  double Fe[N][N];
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Fe[i][j]  = F[i][j];
    }
  }

  double Fa[N][N];
  mat_fun_carray::mat_zero(Fa);
  double Fai[N][N];
  mat_fun_carray::mat_zero(Fai);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) {
        Fa[i][j] = 1.0;
        Fai[i][j] = 1.0;
      }
    }
  }

  double J = mat_fun_carray::mat_det<N>(Fe);
  double J2d = pow(J, (-2.0/nd));
  double J4d = J2d*J2d;

  double Idm[N][N];
  mat_fun_carray::mat_id<N>(Idm);

  double Fe_t[N][N];
  mat_fun_carray::transpose<N>(Fe, Fe_t);

  double C[N][N];
  mat_fun_carray::mat_mul<N>(Fe_t, Fe, C);

  double E[N][N];
  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < nsd; j++) {
      E[i][j] = 0.50 * (C[i][j] - Idm[i][j]);
    }
  }

  double Ci[N][N];
  mat_fun_carray::mat_inv<N>(C, Ci);

  double Cprod[N][N];
  mat_fun_carray::mat_mul<N>(C,C,Cprod);
  double trE = mat_fun_carray::mat_trace<N>(E);
  double Inv1 = J2d * mat_fun_carray::mat_trace<N>(C);
  double Inv2 = 0.50 * (Inv1*Inv1 - J4d * mat_fun_carray::mat_trace<N>(Cprod));

  // Contribution of dilational penalty terms to S and CC
  double p  = 0.0;
  double pl = 0.0;

  if (!utils::is_zero(Kp)) {
    mat_models::get_svol_p(com_mod, cep_mod, stM, J, p, pl);
  }

  // Now, compute isochoric and total stress, elasticity tensors
  //
  double CC[N][N][N][N];
  double Idm_prod[N][N][N][N];
  double Ids[N][N][N][N];
  mat_fun_carray::ten_ids<N>(Ids);

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
    //
    case ConstitutiveModelType::stIso_StVK: {
      double g1 = stM.C10;         // lambda
      double g2 = stM.C01 * 2.0;   // 2*mu

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          S[i][j] = g1 * trE * Idm[i][j]  +  g2 * E[i][j];
        }
      }

      mat_fun_carray::ten_dyad_prod<N>(Idm, Idm, Idm_prod);
      mat_fun_carray::ten_ids<N>(Ids);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CC[i][j][k][l] = g1 * Idm_prod[i][j][k][l]  +  g2 * Ids[i][j][k][l];
            }
          }
        }
      }
    } break;

    // modified St.Venant-Kirchhoff
    //
    case ConstitutiveModelType::stIso_mStVK: {
      double g1 = stM.C10; // kappa
      double g2 = stM.C01;  // mu

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] = g1 * log(J) * Ci[i][j]  +  g2 * (C[i][j] - Idm[i][j]);
        }
      }

      mat_fun_carray::ten_ids<N>(Ids);

      double Ci_Ci_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);

      double Ci_Ci_symm_prod[N][N][N][N];
      mat_fun_carray::ten_symm_prod<N>(Ci, Ci, Ci_Ci_symm_prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CC[i][j][k][l] = g1 * ( -2.0 * log(J) * Ci_Ci_symm_prod[i][j][k][l] + Ci_Ci_prod[i][j][k][l] ) + 
                  2.0 * g2 * Ids[i][j][k][l];
            }
          }
        }
      }

    } break;

    // NeoHookean model
    //
    case ConstitutiveModelType::stIso_nHook: {
      double g1 = 2.0 * stM.C10;

      double Sb[N][N];
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] = g1*Idm[i][j];
        }
      }

      // Fiber reinforcement/active stress
      double prod[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(0), fl.col(0), prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] += Tfa * prod[i][j];
        }
      }

      double r1 = g1 * Inv1 / nd;
      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          S[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }

      double Ci_S_prod[N][N][N][N];
      double S_Ci_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, S, Ci_S_prod);
      mat_fun_carray::ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          S[i][j] += p * J * Ci[i][j];
        }
      }

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] = (-2.0/nd) * (Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l]);
            }
          }
        }
      }

      double Ci_sym_prod[N][N][N][N];
      mat_fun_carray::ten_symm_prod<N>(Ci, Ci, Ci_sym_prod);

      double Ci_Ci_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] += 2.0*(r1 - p*J) * Ci_sym_prod[i][j][k][l] +  (pl*J - 2.0*r1/nd) * Ci_Ci_prod[i][j][k][l];
            }
          }
        }
      }

    } break;

    //  Mooney-Rivlin model
    //
    case ConstitutiveModelType::stIso_MR: {
      double g1  = 2.0 * (stM.C10 + Inv1*stM.C01);
      double g2  = -2.0 * stM.C01;

      double Sb[N][N];
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] = g1*Idm[i][j] + g2*J2d*C[i][j];
        }
      }

      // Fiber reinforcement/active stress
      double prod[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(0), fl.col(0), prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          Sb[i][j] += Tfa*prod[i][j];
        }
      }

      g1 = 4.0 * J4d * stM.C01;

      double CCb[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Idm, Idm, Idm_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CCb[i][j][k][l] = g1 * (Idm_prod[i][j][k][l] - Ids[i][j][k][l]);
            }
          }
        }
      }

      double r1 = J2d * mat_fun_carray::mat_ddot<N>(C, Sb) / nd;
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }

      double Ci_Ci_prod[N][N][N][N];
      double Ci_C_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);
      mat_fun_carray::ten_dyad_prod<N>(Ci, C, Ci_C_prod);

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

      mat_fun_carray::ten_ddot<N>(CCb, PP, CC);

      double CC_t[N][N][N][N];
      mat_fun_carray::ten_transpose<N>(CC, CC_t);

      mat_fun_carray::ten_ddot<N>(PP, CC_t, CC);

      double Ci_S_prod[N][N][N][N];
      double S_Ci_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, S, Ci_S_prod);
      mat_fun_carray::ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] -= (2.0/nd) * ( Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l]);
            }
          }
        }
      }

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] += p * J * Ci[i][j];
        }
      }

      double Ci_sym_prod[N][N][N][N];
      mat_fun_carray::ten_symm_prod<N>(Ci, Ci, Ci_sym_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] += 2.0*(r1 - p*J) * Ci_sym_prod[i][j][k][l]  + (pl*J - 2.0*r1/nd) * Ci_Ci_prod[i][j][k][l];
            }
          }
        }
      }

    } break;

    // HGO (Holzapfel-Gasser-Ogden) model with additive splitting of
    // the anisotropic fiber-based strain-energy terms
    //
    case ConstitutiveModelType::stIso_HGO: {
      if (nfd != 2) {
        //err = "Min fiber directions not defined for HGO material model (2)"
      }
      double kap = stM.kap;

      double C_fl[N];
      mat_fun_carray::mat_mul(C, fl.rcol(0), C_fl);
      double Inv4 = J2d * mat_fun_carray::norm<N>(fl.rcol(0), C_fl);

      mat_fun_carray::mat_mul(C, fl.rcol(1), C_fl);
      double Inv6 = J2d * mat_fun_carray::norm<N>(fl.rcol(1), C_fl);

      double Eff = kap*Inv1 + (1.0-3.0*kap)*Inv4 - 1.0;
      double Ess = kap*Inv1 + (1.0-3.0*kap)*Inv6 - 1.0;

      double Hff[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(0), fl.col(0), Hff);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Hff[i][j] = kap*Idm[i][j] + (1.0-3.0*kap)*Hff[i][j];
        }
      }

      double Hss[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(1), fl.col(1), Hss);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Hss[i][j] = kap*Idm[i][j] + (1.0-3.0*kap)*Hss[i][j];
        }
      }

      double g1 = stM.C10;
      double g2 = stM.aff * Eff * exp(stM.bff*Eff*Eff);
      double g3 = stM.ass * Ess * exp(stM.bss*Ess*Ess);

      double Sb[N][N];
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] = 2.0 * (g1*Idm[i][j] + g2*Hff[i][j] + g3*Hss[i][j]);
        }
      }

      // Fiber reinforcement/active stress
      double prod[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(0), fl.col(0), prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] += Tfa*prod[i][j];
        }
      }

      g1 = stM.aff*(1.0 + 2.0*stM.bff*Eff*Eff)*exp(stM.bff*Eff*Eff);
      g2 = stM.ass*(1.0 + 2.0*stM.bss*Ess*Ess)*exp(stM.bss*Ess*Ess);
      g1 = 4.0*J4d*g1;
      g2 = 4.0*J4d*g2;

      double CCb_hff[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Hff, Hff, CCb_hff);

      double CCb_hss[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Hss, Hss, CCb_hss);

      double CCb[N][N][N][N];

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CCb[i][j][k][l] = g1 * CCb_hff[i][j][k][l]  +  g2 * CCb_hss[i][j][k][l];
            }
          }
        }
      }

      double r1 = J2d * mat_fun_carray::mat_ddot<N>(C, Sb) / nd;

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          S[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }

      double Ci_C_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, C, Ci_C_prod);
      double PP[N][N][N][N];

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              PP[i][j][k][l] = Ids[i][j][k][l] - (1.0/nd) * Ci_C_prod[i][j][k][l];
            }
          }
        }
      }

      mat_fun_carray::ten_ddot<N>(CCb, PP, CC);

      double CC_t[N][N][N][N];
      mat_fun_carray::ten_transpose<N>(CC, CC_t);

      mat_fun_carray::ten_ddot<N>(PP, CC_t, CC);

      double Ci_S_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, S, Ci_S_prod);

      double S_Ci_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CC[i][j][k][l] -= (2.0/nd) * ( Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l]);
            }
          }
        }
      }

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] += p * J * Ci[i][j];
        }
      }

      double Ci_sym_prod[N][N][N][N];
      mat_fun_carray::ten_symm_prod<N>(Ci, Ci, Ci_sym_prod);

      double Ci_Ci_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] += 2.0*(r1 - p*J) * Ci_sym_prod[i][j][k][l]  + (pl*J - 2.0*r1/nd) * Ci_Ci_prod[i][j][k][l];
            }
          }
        }
      }

    } break;

    // Guccione (1995) transversely isotropic model
    //
    case ConstitutiveModelType::stIso_Gucci: {
      #ifdef debug_get_pk2cc
      std::cout << msg_prefix << "stIso_Gucci" << std::endl;
      #endif
      if (nfd != 2) {
        //err = "Min fiber directions not defined for Guccione material model (2)"
      }

      // Compute isochoric component of E
      double E[N][N];
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          E[i][j] = 0.50 * (J2d*C[i][j] - Idm[i][j]);
        }
      }

      // Transform into local orthogonal coordinate system
      //
      double Rm[N][N];
      auto fl_0 = fl.col(0);
      auto fl_1 = fl.col(1);
      auto fl_cross = cross(fl);

      for (int i = 0; i < N; i++) {
        Rm[i][0] = fl_0[i];
        Rm[i][1] = fl_1[i];
        Rm[i][2] = fl_cross[i];
      }

      // Project E to local orthogocal coordinate system
      //
      double Es_1[N][N];
      mat_fun_carray::mat_mul(E, Rm, Es_1);

      double Rm_t[N][N];
      mat_fun_carray::transpose<N>(Rm, Rm_t);

      double Es[N][N];
      mat_fun_carray::mat_mul(Rm_t, Es_1, Es);

      double g1 = stM.bff;
      double g2 = stM.bss;
      double g3 = stM.bfs;

      double QQ = g1 *  Es[0][0]*Es[0][0] + 
                  g2 * (Es[1][1]*Es[1][1] + Es[2][2]*Es[2][2] + Es[1][2]*Es[1][2] + Es[2][1]*Es[2][1]) +
                  g3 * (Es[0][1]*Es[0][1] + Es[1][0]*Es[1][0] + Es[0][2]*Es[0][2] + Es[2][0]*Es[2][0]);

      double r2 = stM.C10 * exp(QQ);

      // Fiber stiffness contribution := (dE*_ab / dE_IJ)
      //
      // Use Array in the following since it is a bit tricky.
      //
      Array<double> Rm_a(N,N);
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Rm_a(i,j) = Rm[i][j];
        }
      }

      Array3<double> RmRm(N,N,6);
     
      RmRm.set_slice(0, mat_dyad_prod(Rm_a.col(0), Rm_a.col(0), nsd));
      RmRm.set_slice(1, mat_dyad_prod(Rm_a.col(1), Rm_a.col(1), nsd));
      RmRm.set_slice(2, mat_dyad_prod(Rm_a.col(2), Rm_a.col(2), nsd));

      RmRm.set_slice(3, mat_symm_prod(Rm_a.col(0), Rm_a.col(1), nsd));
      RmRm.set_slice(4, mat_symm_prod(Rm_a.col(1), Rm_a.col(2), nsd));
      RmRm.set_slice(5, mat_symm_prod(Rm_a.col(2), Rm_a.col(0), nsd));

      double Sb[N][N];

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] = g1 *  Es[0][0] * RmRm.rslice(0)(i,j) + 
                     g2 * (Es[1][1] * RmRm.rslice(1)(i,j) + Es[2][2]*RmRm.rslice(2)(i,j) + 
                           2.0*Es[1][2]*RmRm.rslice(4)(i,j)) +
               2.0 * g3 * (Es[0][1] * RmRm.rslice(3)(i,j) + Es[0][2]*RmRm.rslice(5)(i,j));

        }
      }

      double Sb_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Sb, Sb, Sb_prod);
      double CCb[N][N][N][N];

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CCb[i][j][k][l] = 2.0 * Sb_prod[i][j][k][l];
            }
          }
        }
      }

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] = Sb[i][j] * r2;
        }
      }

      // Fiber reinforcement/active stress
      //
      double prod[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(0), fl.col(0), prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] += Tfa * prod[i][j];
        }
      }

      double r1 = J2d * mat_fun_carray::mat_ddot<N>(C, Sb) / nd;

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }

      CArray4 prod_0;
      CArray4 prod_1;
      CArray4 prod_2;
      CArray4 prod_3;
      CArray4 prod_4;
      CArray4 prod_5;

      mat_fun_carray::ten_dyad_prod<N>(RmRm.slice(0), RmRm.slice(0), prod_0);
      mat_fun_carray::ten_dyad_prod<N>(RmRm.slice(1), RmRm.slice(1), prod_1);
      mat_fun_carray::ten_dyad_prod<N>(RmRm.slice(2), RmRm.slice(2), prod_2);
      mat_fun_carray::ten_dyad_prod<N>(RmRm.slice(3), RmRm.slice(3), prod_3);
      mat_fun_carray::ten_dyad_prod<N>(RmRm.slice(4), RmRm.slice(4), prod_4);
      mat_fun_carray::ten_dyad_prod<N>(RmRm.slice(5), RmRm.slice(5), prod_5);

      r2 = r2 * J4d;

      for (int i = 0; i < N; i++) { 
        for (int j = 0; j < N; j++) { 
          for (int k = 0; k < N; k++) { 
            for (int l = 0; l < N; l++) { 
              CCb[i][j][k][l] = r2 * (CCb[i][j][k][l] + g1 * prod_0[i][j][k][l] +
                         g2 * ( prod_1[i][j][k][l] + prod_2[i][j][k][l] + 2.0*prod_4[i][j][k][l] ) +
                         2.0 * g3 * (prod_3[i][j][k][l] + prod_5[i][j][k][l] ));
            }
          }
        }
      }

      double Ci_C_prod[N][N][N][N];
      mat_fun_carray::ten_dyad_prod<N>(Ci, C, Ci_C_prod);
      CArray4 PP;

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              PP[i][j][k][l] = Ids[i][j][k][l] - (1.0/nd) * Ci_C_prod[i][j][k][l];
            }
          }
        }
      }

      mat_fun_carray::ten_ddot<N>(CCb, PP, CC);

      double CC_t[N][N][N][N];
      mat_fun_carray::ten_transpose<N>(CC, CC_t);

      mat_fun_carray::ten_ddot<N>(PP, CC_t, CC);

      CArray4 Ci_S_prod;
      mat_fun_carray::ten_dyad_prod<N>(Ci, S, Ci_S_prod);

      CArray4 S_Ci_prod;
      mat_fun_carray::ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CC[i][j][k][l] -= (2.0/nd) * ( Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l]);
            }
          }
        }
      }

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] += p * J * Ci[i][j];
        }
      }

      CArray4 Ci_sym_prod;
      mat_fun_carray::ten_symm_prod<N>(Ci, Ci, Ci_sym_prod);

      CArray4 Ci_Ci_prod;
      mat_fun_carray::ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] += 2.0*(r1 - p*J) * Ci_sym_prod[i][j][k][l]  + (pl*J - 2.0*r1/nd) * Ci_Ci_prod[i][j][k][l];
            }
          }
        }
      }

    } break;

    //  HO (Holzapfel-Ogden) model for myocardium (2009)
    case ConstitutiveModelType::stIso_HO: {
      if (nfd != 2) {
        //err = "Min fiber directions not defined for Holzapfel material model (2)"
      }

      // Compute fiber-based isochoric invariants
      double C_fl[N];
      mat_fun_carray::mat_mul(C, fl.rcol(0), C_fl);
      double Inv4 = J2d * mat_fun_carray::norm<N>(fl.rcol(0), C_fl);

      mat_fun_carray::mat_mul(C, fl.rcol(1), C_fl);
      double Inv6 = J2d * mat_fun_carray::norm<N>(fl.rcol(1), C_fl);

      mat_fun_carray::mat_mul(C, fl.rcol(1), C_fl);
      double Inv8 = J2d * mat_fun_carray::norm<N>(fl.rcol(0), C_fl);

      double Eff = Inv4 - 1.0;
      double Ess = Inv6 - 1.0;
      double Efs = Inv8;

      // Smoothed Heaviside function: 1 / (1 + exp(-kx)) = 1 - 1 / (1 + exp(kx))
      double k = stM.khs;
      double one_over_exp_plus_one_f = 1.0 / (exp(k * Eff) + 1.0);
      double one_over_exp_plus_one_s = 1.0 / (exp(k * Ess) + 1.0);
      double c4f  = 1.0 - one_over_exp_plus_one_f;
      double c4s  = 1.0 - one_over_exp_plus_one_s;
      
      // Approx. derivative of smoothed heaviside function
      // double dc4f = 0.25*stM.khs*exp(-stM.khs*abs(Eff));
      // double dc4s = 0.25*stM.khs*exp(-stM.khs*abs(Ess));

      // Exact first derivative of smoothed heaviside function (from Wolfram Alpha)
      double dc4f = k * (one_over_exp_plus_one_f - pow(one_over_exp_plus_one_f,2));
      double dc4s = k * (one_over_exp_plus_one_s - pow(one_over_exp_plus_one_s,2));

      // Exact second derivative of smoothed heaviside function (from Wolfram Alpha)
      double ddc4f = pow(k,2) * (-one_over_exp_plus_one_f + 3.0*pow(one_over_exp_plus_one_f,2) - 2.0*pow(one_over_exp_plus_one_f,3));
      double ddc4s = pow(k,2) * (-one_over_exp_plus_one_s + 3.0*pow(one_over_exp_plus_one_s,2) - 2.0*pow(one_over_exp_plus_one_s,3));
      
      // Isotropic + fiber-sheet interaction stress
      double g1 = stM.a * exp(stM.b*(Inv1-3.0));
      double g2 = 2.0 * stM.afs * exp(stM.bfs*Efs*Efs);

      double Hfs[N][N];
      mat_fun_carray::mat_symm_prod<N>(fl.rcol(0), fl.rcol(1), Hfs);

      double Sb[N][N];
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] = g1*Idm[i][j] + g2*Efs*Hfs[i][j];
        }
      }

      // Isotropic + fiber-sheet interaction stiffness
      g1 = g1 * 2.0 * J4d * stM.b;
      g2 = g2 * 2.0 * J4d * (1.0 + 2.0*stM.bfs*Efs*Efs);

      CArray4 Idm_prod;
      mat_fun_carray::ten_dyad_prod<N>(Idm, Idm, Idm_prod);

      CArray4 Hfs_prod;
      mat_fun_carray::ten_dyad_prod<N>(Hfs, Hfs, Hfs_prod);

      CArray4 CCb;

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CCb[i][j][k][l] = g1 * Idm_prod[i][j][k][l]  +  g2 * Hfs_prod[i][j][k][l];
            }
          }
        }
      }

      //  Fiber-fiber interaction stress + additional reinforcement (Tfa)
      double rexp = exp(stM.bff*Eff*Eff);
      g1 = c4f * Eff * rexp;
      g1 = g1 + (0.5*dc4f/stM.bff) * (rexp - 1.0);
      g1 = 2.0 * stM.aff * g1 + Tfa;
      double Hff[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(0), fl.col(0), Hff);
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            Sb[i][j] += g1 * Hff[i][j]; 
          }
        }

      // Fiber-fiber interaction stiffness
      g1 = c4f * (1.0 + 2.0*stM.bff*Eff*Eff);
      g1 = (g1 + 2.0*dc4f*Eff) * rexp;
      g1 = g1 + (0.5*ddc4f/stM.bff)*(rexp - 1.0);
      g1 = 4.0 * J4d * stM.aff * g1;
      CArray4 Hff_prod;
      mat_fun_carray::ten_dyad_prod<N>(Hff, Hff, Hff_prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CCb[i][j][k][l] += g1 * Hff_prod[i][j][k][l];
            }
          }
        }
      }

      // Sheet-sheet interaction stress + additional cross-fiber stress
      rexp = exp(stM.bss*Ess*Ess);
      g2 = c4s * Ess * rexp;
      g2 = g2 + (0.5*dc4s/stM.bss) * (rexp - 1.0);
      g2 = 2.0 * stM.ass * g2 + Tsa;
      double Hss[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.col(1), fl.col(1), Hss);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] += g2 * Hss[i][j];
        }
      }

      // Sheet-sheet interaction stiffness
      g2 = c4s * (1.0 + 2.0 * stM.bss * Ess * Ess);
      g2 = (g2 + 2.0*dc4s*Ess) * rexp;
      g2 = g2 + (0.5*ddc4s/stM.bss)*(rexp - 1.0);
      g2 = 4.0 * J4d * stM.ass * g2;

      CArray4 Hss_prod;
      mat_fun_carray::ten_dyad_prod<N>(Hss, Hss, Hss_prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CCb[i][j][k][l] += g2 * Hss_prod[i][j][k][l];
            }
          }
        }
      }

      // Isochoric 2nd-Piola-Kirchhoff stress and stiffness tensors
      double r1 = J2d * mat_fun_carray::mat_ddot<N>(C, Sb) / nd;

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }

      CArray4 Ci_C_prod;
      mat_fun_carray::ten_dyad_prod<N>(Ci, C, Ci_C_prod);
      double PP[N][N][N][N];

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              PP[i][j][k][l] = Ids[i][j][k][l] - (1.0/nd) * Ci_C_prod[i][j][k][l];
            }
          }
        }
      }

      mat_fun_carray::ten_ddot<N>(CCb, PP, CC);

      CArray4 CC_t;
      mat_fun_carray::ten_transpose<N>(CC, CC_t);

      mat_fun_carray::ten_ddot<N>(PP, CC_t, CC);

      CArray4 Ci_S_prod;
      mat_fun_carray::ten_dyad_prod<N>(Ci, S, Ci_S_prod);

      CArray4 S_Ci_prod;
      mat_fun_carray::ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CC[i][j][k][l] -= (2.0/nd) * (Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l]);
            }
          }
        }
      }

      // Add pressure contribution
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          S[i][j] += p * J * Ci[i][j];
        }
      }

      CArray4 Ci_sym_prod;
      mat_fun_carray::ten_symm_prod<N>(Ci, Ci, Ci_sym_prod);

      CArray4 Ci_Ci_prod;
      mat_fun_carray::ten_dyad_prod<N>(Ci, Ci, Ci_Ci_prod);

      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          for (int k = 0; k < nsd; k++) {
            for (int l = 0; l < nsd; l++) {
              CC[i][j][k][l] += 2.0*(r1 - p*J) * Ci_sym_prod[i][j][k][l]  + (pl*J - 2.0*r1/nd) * Ci_Ci_prod[i][j][k][l];
            }
          }
        }
      }

      if (cep_mod.cem.aStrain) {
        double S_prod[N][N];
        mat_fun_carray::mat_mul<N>(Fai, S, S_prod);
        //S = mat_mul(Fai, S);

        double Fai_t[N][N];
        mat_fun_carray::transpose<N>(Fai, Fai_t);
        mat_fun_carray::mat_mul<N>(S_prod, Fai_t, S);
        //S = mat_mul(S, transpose(Fai));

        mat_fun_carray::ten_dyad_prod<N>(Fai, Fai, CCb);
        //CCb = ten_dyad_prod(Fai, Fai, nsd);

        CArray4 CC_dot;
        mat_fun_carray::ten_ddot_3424<N>(CC, CCb, CC_dot);
        //CC = ten_ddot_3424(CC, CCb, nsd);

        mat_fun_carray::ten_ddot_2412<N>(CCb, CC_dot, CC);
        //CC = ten_ddot_2412(CCb, CC, nsd);
      }

    } break;

//  HO (Holzapfel-Ogden)-MA model for myocardium with full invariants for the anisotropy terms (modified-anisotropy)
    case ConstitutiveModelType::stIso_HO_ma: {
      if (nfd != 2) {
        //err = "Min fiber directions not defined for Holzapfel material model (2)"
      }

      // Compute fiber-based full invariants (not isochoric)
      double C_fl[N];
      mat_fun_carray::mat_mul(C, fl.rcol(0), C_fl);
      double Inv4 = mat_fun_carray::norm<N>(fl.rcol(0), C_fl);

      mat_fun_carray::mat_mul(C, fl.rcol(1), C_fl);
      double Inv6 = mat_fun_carray::norm<N>(fl.rcol(1), C_fl);

      mat_fun_carray::mat_mul(C, fl.rcol(1), C_fl);
      double Inv8 = mat_fun_carray::norm<N>(fl.rcol(0), C_fl);

      double Eff = Inv4 - 1.0;
      double Ess = Inv6 - 1.0;
      double Efs = Inv8;

      // Smoothed Heaviside function: 1 / (1 + exp(-kx)) = 1 - 1 / (1 + exp(kx))
      double k = stM.khs;
      double one_over_exp_plus_one_f = 1.0 / (exp(k * Eff) + 1.0);
      double one_over_exp_plus_one_s = 1.0 / (exp(k * Ess) + 1.0);
      double c4f  = 1.0 - one_over_exp_plus_one_f;
      double c4s  = 1.0 - one_over_exp_plus_one_s;
      
      // Approx. derivative of smoothed heaviside function
      // double dc4f = 0.25*stM.khs*exp(-stM.khs*abs(Eff));
      // double dc4s = 0.25*stM.khs*exp(-stM.khs*abs(Ess));

      // Exact first derivative of smoothed heaviside function (from Wolfram Alpha)
      double dc4f = k * (one_over_exp_plus_one_f - pow(one_over_exp_plus_one_f,2));
      double dc4s = k * (one_over_exp_plus_one_s - pow(one_over_exp_plus_one_s,2));

      // Exact second derivative of smoothed heaviside function (from Wolfram Alpha)
      double ddc4f = pow(k,2) * (-one_over_exp_plus_one_f + 3.0*pow(one_over_exp_plus_one_f,2) - 2.0*pow(one_over_exp_plus_one_f,3));
      double ddc4s = pow(k,2) * (-one_over_exp_plus_one_s + 3.0*pow(one_over_exp_plus_one_s,2) - 2.0*pow(one_over_exp_plus_one_s,3));
      
      // Isochoric stress and stiffness
      double g1 = stM.a * exp(stM.b*(Inv1-3.0));

      double Sb[N][N];
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          Sb[i][j] = g1*Idm[i][j];
        }
      }
      double r1 = J2d/nd*mat_fun_carray::mat_ddot(C,Sb);
      g1 = g1*2.0*J4d*stM.b;

      CArray4 Idm_prod;
      mat_fun_carray::ten_dyad_prod<N>(Idm, Idm, Idm_prod);

      CArray4 CCb;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CCb[i][j][k][l] = g1 * Idm_prod[i][j][k][l];
            }
          }
        }
      }

      // Add isochoric stress and stiffness contribution
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          S[i][j] = J2d*Sb[i][j] - r1*Ci[i][j];
        }
      }
      CArray4 Ci_C_prod;
      mat_fun_carray::ten_dyad_prod<N>(Ci, C, Ci_C_prod);
      double PP[N][N][N][N];

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              PP[i][j][k][l] = Ids[i][j][k][l] - (1.0/nd) * Ci_C_prod[i][j][k][l];
            }
          }
        }
      }

      mat_fun_carray::ten_ddot<N>(CCb, PP, CC);

      CArray4 CC_t;
      mat_fun_carray::ten_transpose<N>(CC, CC_t);

      mat_fun_carray::ten_ddot<N>(PP, CC_t, CC);

      CArray4 Ci_S_prod;
      mat_fun_carray::ten_dyad_prod<N>(Ci, S, Ci_S_prod);

      CArray4 S_Ci_prod;
      mat_fun_carray::ten_dyad_prod<N>(S, Ci, S_Ci_prod);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
              CC[i][j][k][l] -= (2.0/nd) * (Ci_S_prod[i][j][k][l] + S_Ci_prod[i][j][k][l]);
            }
          }
        }
      }

      // Add pressure contribution to stress and stiffness
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            S[i][j] += p*J*Ci[i][j];
          }
        }
    
      CArray4 Ci_Ci_prod;
      mat_fun_carray::ten_dyad_prod(Ci, Ci, Ci_Ci_prod);
      CArray4 Ci_Ci_symprod;
      mat_fun_carray::ten_symm_prod(Ci, Ci, Ci_Ci_symprod);
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
              for (int l = 0; l < N; l++) {
                CC[i][j][k][l] += 2.0*(r1 - p*J)*Ci_Ci_symprod[i][j][k][l] + (pl*J - 2.0*r1/nd)*Ci_Ci_prod[i][j][k][l];
              }
            }
          }
        }

      // Now that both isochoric and volumetric components were added, anisotropic components need to be added

      // Fiber-sheet interaction terms
      g1   = 2.0 * stM.afs * exp(stM.bfs*Efs*Efs);
      double Hfs[N][N];
      mat_fun_carray::mat_symm_prod<N>(fl.rcol(0), fl.rcol(1), Hfs);
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            S[i][j] += g1*Efs*Hfs[i][j];
          }
        }

      g1   = g1 * 2.0*(1.0 + 2.0*stM.bfs*Efs*Efs);

      CArray4 Hfs_Hfs_prod;
      mat_fun_carray::ten_dyad_prod(Hfs, Hfs, Hfs_Hfs_prod);
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
              for (int l = 0; l < N; l++) {
                CC[i][j][k][l] += g1*Hfs_Hfs_prod[i][j][k][l];
              }
            }
          }
        }

      // Fiber-fiber interaction stress + additional reinforcement (Tfa)
      double rexp = exp(stM.bff * Eff * Eff);
      g1   = c4f*Eff*rexp;
      g1   = g1 + (0.5*dc4f/stM.bff)*(rexp - 1.0);
      g1   = (2.0*stM.aff*g1) + Tfa;
      double Hff[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.rcol(0), fl.rcol(0), Hff);
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          S[i][j] += g1*Hff[i][j];
        }
      }
 
      // Fiber-fiber interaction stiffness
      g1   = c4f*(1.0 + (2.0*stM.bff*Eff*Eff));
      g1   = (g1 + (2.0*dc4f*Eff))*rexp;
      g1   = g1 + (0.5*ddc4f/stM.bff)*(rexp - 1.0);
      g1   = 4.0*stM.aff*g1;
      CArray4 Hff_Hff_prod;
      mat_fun_carray::ten_dyad_prod(Hff, Hff, Hff_Hff_prod);
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
              for (int l = 0; l < N; l++) {
                CC[i][j][k][l] += g1*Hff_Hff_prod[i][j][k][l];
              }
            }
          }
        }
         
      // Sheet-sheet interaction stress + additional cross-fiber stress
      rexp = exp(stM.bss * Ess * Ess);
      double g2   = c4s*Ess*rexp;
      g2   = g2 + (0.5*dc4s/stM.bss)*(rexp - 1.0);
      g2   = 2.0*stM.ass*g2 + Tsa;
      double Hss[N][N];
      mat_fun_carray::mat_dyad_prod<N>(fl.rcol(1), fl.rcol(1), Hss);
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          S[i][j] += g2*Hss[i][j];
        }
      }
  
      // Sheet-sheet interaction stiffness
      g2   = c4s*(1.0 + (2.0*stM.bss*Ess*Ess));
      g2   = (g2 + (2.0*dc4s*Ess))*rexp;
      g2   = g2 + (0.5*ddc4s/stM.bss)*(rexp - 1.0);
      g2   = 4.0*stM.ass*g2;
      CArray4 Hss_Hss_prod;
      mat_fun_carray::ten_dyad_prod(Hss,Hss,Hss_Hss_prod);
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
              for (int l = 0; l < N; l++) {
                CC[i][j][k][l] += g2*Hss_Hss_prod[i][j][k][l];
              }
            }
          }
        }
      } break;

    default:
      throw std::runtime_error("Undefined material constitutive model.");
  } 

  // Convert to Voigt Notation
  cc_to_voigt_carray<N>(CC, Dm);
}


/**
 * @brief Get the viscous PK2 stress and corresponding tangent matrix contributions for a solid
 * with a viscous pseudo-potential model.
 * This is defined by a viscous pseuo-potential
 * Psi = mu/2 * tr(E_dot^2)
 * The viscous 2nd Piola-Kirchhoff stress is given by
 * Svis = dPsi/dE_dot 
 *   = mu * E_dot
 *   = mu * 1/2 * F^T * (grad(v) + grad(v)^T) * F
 *   = mu * 1/2 * ( (F^T * Grad(v)) + (F^T * Grad(v))^T )
 * 
 * @tparam nsd Number of spatial dimensions
 * @param mu Solid viscosity parameter
 * @param eNoN Number of nodes in an element
 * @param Nx Shape function gradient w.r.t. reference configuration coordinates (dN/dX)
 * @param vx Velocity gradient matrix w.r.t reference configuration coordinates (dv/dX)
 * @param F Deformation gradient matrix
 * @param Svis Viscous 2nd Piola-Kirchhoff stress matrix
 * @param Kvis_u Viscous tangent matrix contribution due to displacement
 * @param Kvis_v Visous tangent matrix contribution due to velocity
 */
template <size_t nsd>
void get_visc_stress_pot(const double mu, const int eNoN, const Array<double>& Nx, const double vx[nsd][nsd], const double F[nsd][nsd],
                        Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v) {

    
    // Initialize Svis, Kvis_u, Kvis_v to zero
    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            Svis(i,j) = 0.0;
            for (int a = 0; a < eNoN; a++) {
                for (int b = 0; b < eNoN; b++) {
                    Kvis_u(i*nsd+j,a,b) = 0.0;
                    Kvis_v(i*nsd+j,a,b) = 0.0;
                }
            }
        }
    }


    // Required intermediate terms for stress and tangent
    double Ft[nsd][nsd] = {0}, vxt[nsd][nsd] = {0}, F_Ft[nsd][nsd] = {0}, Ft_vx[nsd][nsd] = {0}, F_vxt[nsd][nsd] = {0};
    mat_fun_carray::transpose<nsd>(F, Ft);
    mat_fun_carray::mat_mul<nsd>(F, Ft, F_Ft);
    mat_fun_carray::mat_mul<nsd>(Ft, vx, Ft_vx);
    mat_fun_carray::transpose<nsd>(vx, vxt);
    mat_fun_carray::mat_mul<nsd>(F, vxt, F_vxt);

    //double F_Nx[nsd][eNoN] = {0}, vx_Nx[nsd][eNoN] = {0};
    Array<double> F_Nx(nsd,eNoN), vx_Nx(nsd,eNoN);
    
    for (int a = 0; a < eNoN; ++a) {
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                F_Nx(i,a) += F[i][j] * Nx(j,a);
                vx_Nx(i,a) += vx[i][j] * Nx(j,a);
            }
        }
    }

    // 2nd Piola-Kirchhoff stress due to viscosity
    // Svis = mu * 1/2 * ( (F^T * dv/dX) + (F^T * dv/dX)^T )
    double Ft_vx_symm[nsd][nsd] = {0};
    mat_fun_carray::mat_symm<nsd>(Ft_vx, Ft_vx_symm);
    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            Svis(i,j) = mu * Ft_vx_symm[i][j];
        }
    }

    // Tangent matrix contributions due to viscosity
    for (int b = 0; b < eNoN; ++b) {
        for (int a = 0; a < eNoN; ++a) {
            double Nx_Nx = 0.0;
            for (int i = 0; i < nsd; ++i) {
                Nx_Nx += Nx(i,a) * Nx(i,b);
            }

            for (int i = 0; i < nsd; ++i) {
                for (int j = 0; j < nsd; ++j) {
                    int ii = i * nsd + j;
                    Kvis_u(ii,a,b) = 0.5 * mu * (F_Nx(i,b) * vx_Nx(j,a) + Nx_Nx * F_vxt[i][j]);
                    Kvis_v(ii,a,b) = 0.5 * mu * (Nx_Nx * F_Ft[i][j] + F_Nx(i,b) * F_Nx(j,a));
                }
            }
        }
    }
}

/**
 * @brief Get the viscous PK2 stress and corresponding tangent matrix contributions for a solid
 * with a Newtonian fluid-like viscosity model.
 * The viscous deviatoric Cauchy stress is given by
 * sigma_vis_dev = 2 * mu * d_dev
 * where d_dev = 1/2 * (grad(v) + grad(v)^T) - 1/3 * (div(v)) * I
 * The viscous 2nd Piola-Kirchhoff stress is given by a pull-back operation
 * Svis = 2 * mu * J * F^-1 * d_dev * F^-T
 * 
 * Note, there is likely an error/bug in the tangent contributions that leads to suboptimal nonlinear convergence
 * 
 * @tparam nsd Number of spatial dimensions
 * @param mu Solid viscosity parameter
 * @param eNoN Number of nodes in an element
 * @param Nx Shape function gradient w.r.t. reference configuration coordinates (dN/dX)
 * @param vx Velocity gradient matrix w.r.t reference configuration coordinates (dv/dX)
 * @param F Deformation gradient matrix
 * @param Svis Viscous 2nd Piola-Kirchhoff stress matrix
 * @param Kvis_u Viscous tangent matrix contribution due to displacement
 * @param Kvis_v Visous tangent matrix contribution due to velocity
 */
template <size_t nsd>
void get_visc_stress_newt(const double mu, const int eNoN, const Array<double>& Nx, const double vx[nsd][nsd], const double F[nsd][nsd],
                           Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v) {
    

    // Initialize Svis, Kvis_u, Kvis_v to zero
    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            Svis(i,j) = 0.0;
            for (int a = 0; a < eNoN; a++) {
                for (int b = 0; b < eNoN; b++) {
                    Kvis_u(i*nsd+j,a,b) = 0.0;
                    Kvis_v(i*nsd+j,a,b) = 0.0;
                }
            }
        }
    }

    // Get identity matrix, Jacobian, and F^-1
    double Idm[nsd][nsd] = {0};
    mat_fun_carray::mat_id<nsd>(Idm);
    double J = mat_fun_carray::mat_det<nsd>(F);
    double Fi[nsd][nsd] = {0};
    mat_fun_carray::mat_inv<nsd>(F, Fi); 

    // Required intermediate terms for stress and tangent
    double vx_Fi[nsd][nsd] = {0}, vx_Fi_symm[nsd][nsd] = {0}, ddev[nsd][nsd] = {0};
    // vx_Fi: Velocity gradient in current configuration          
    mat_fun_carray::mat_mul<nsd>(vx, Fi, vx_Fi);        
    mat_fun_carray::mat_symm<nsd>(vx_Fi, vx_Fi_symm);
    // ddev: Deviatoric part of rate of strain tensor
    mat_fun_carray::mat_dev<nsd>(vx_Fi_symm, ddev);
    //double Nx_Fi[nsd][eNoN] = {0}, ddev_Nx_Fi[nsd][eNoN] = {0}, vx_Fi_Nx_Fi[nsd][eNoN] = {0};
    Array<double> Nx_Fi(nsd,eNoN), ddev_Nx_Fi(nsd,eNoN), vx_Fi_Nx_Fi(nsd,eNoN);
    for (int a = 0; a < eNoN; ++a) {
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                Nx_Fi(i,a) += Nx(j,a) * Fi[j][i];
            }
        }
        for (int i = 0; i < nsd; ++i) {
          for (int j = 0; j < nsd; ++j) {
              ddev_Nx_Fi(i,a) += ddev[i][j] * Nx_Fi(j,a);
              vx_Fi_Nx_Fi(i,a) += vx_Fi[j][i] * Nx_Fi(j,a);
          }
        }
    }

    // 2nd Piola-Kirchhoff stress due to viscosity
    // Svis = 2 * mu * J * F^-1 * d_dev * F^-T
    double Fit[nsd][nsd] = {0};
    mat_fun_carray::transpose<nsd>(Fi, Fit);
    double ddev_Fit[nsd][nsd] = {0};
    mat_fun_carray::mat_mul<nsd>(ddev, Fit, ddev_Fit);
    double Fi_ddev_Fit[nsd][nsd] = {0};
    mat_fun_carray::mat_mul<nsd>(Fi, ddev_Fit, Fi_ddev_Fit);
    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            Svis(i,j) = 2.0 * mu * J * Fi_ddev_Fit[i][j];
        }
    }

    // Tangent matrix contributions due to viscosity
    double r2d = 2.0 / nsd;
    for (int b = 0; b < eNoN; ++b) {
        for (int a = 0; a < eNoN; ++a) {
            double Nx_Fi_Nx_Fi = 0.0;
            for (int i = 0; i < nsd; ++i) {
                Nx_Fi_Nx_Fi += Nx_Fi(i,a) * Nx_Fi(i,b);
            }

            for (int i = 0; i < nsd; ++i) {
                for (int j = 0; j < nsd; ++j) {
                    int ii = i * nsd + j;

                    // Derivative of the residual w.r.t displacement
                    Kvis_u(ii,a,b) = mu * J * (2.0 * 
                                    (ddev_Nx_Fi(i,a) * Nx_Fi(j,b) - ddev_Nx_Fi(i,b) * Nx_Fi(j,a)) -
                                    (Nx_Fi_Nx_Fi * vx_Fi[i][j] + Nx_Fi(i,b) * vx_Fi_Nx_Fi(j,a) -
                                    r2d * Nx_Fi(i,a) * vx_Fi_Nx_Fi(j,b)));

                    // Derivative of the residual w.r.t velocity
                    Kvis_v(ii,a,b) = mu * J * (Nx_Fi_Nx_Fi * Idm[i][j] +
                                    Nx_Fi(i,b) * Nx_Fi(j,a) - r2d * Nx_Fi(i,a) * Nx_Fi(j,b));
                }
            }
        }
    }
}


/**
 * @brief Get the solid viscous PK2 stress and corresponding tangent matrix contributions
 * Calls the appropriate function based on the viscosity type, either viscous 
 * pseudo-potential or Newtonian viscosity model.
 * 
 * @tparam nsd Number of spatial dimensions
 * @param[in] lDmn Domain object
 * @param[in] eNoN Number of nodes in an element
 * @param[in] Nx Shape function gradient w.r.t. reference configuration coordinates (dN/dX)
 * @param[in] vx Velocity gradient matrix w.r.t reference configuration coordinates (dv/dX)
 * @param[in] F Deformation gradient matrix
 * @param[out] Svis Viscous 2nd Piola-Kirchhoff stress matrix
 * @param[out] Kvis_u Viscous tangent matrix contribution due to displacement
 * @param[out] Kvis_v Viscous tangent matrix contribution due to velocity
 */
template <size_t nsd>
void get_visc_stress_and_tangent(const dmnType& lDmn, const int eNoN, const Array<double>& Nx, const double vx[nsd][nsd], const double F[nsd][nsd],
                                 Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v) {

    switch (lDmn.solid_visc.viscType) {
      case consts::SolidViscosityModelType::viscType_Newtonian:
        get_visc_stress_newt<nsd>(lDmn.solid_visc.mu, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);
      break;

      case consts::SolidViscosityModelType::viscType_Potential:
        get_visc_stress_pot<nsd>(lDmn.solid_visc.mu, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);
      break;
    }
}

/**
 * @brief Get the solid viscous PK2 stress and corresponding tangent matrix contributions
 * Calls the appropriate function based on the viscosity type, either viscous 
 * pseudo-potential or Newtonian viscosity model.
 * 
 * Same as above, except takes vx and F as Array objects instead of C-style arrays.
 * 
 * @tparam nsd Number of spatial dimensions
 * @param[in] lDmn Domain object
 * @param[in] eNoN Number of nodes in an element
 * @param[in] Nx Shape function gradient w.r.t. reference configuration coordinates (dN/dX)
 * @param[in] vx Velocity gradient matrix w.r.t reference configuration coordinates (dv/dX)
 * @param[in] F Deformation gradient matrix
 * @param[out] Svis Viscous 2nd Piola-Kirchhoff stress matrix
 * @param[out] Kvis_u Viscous tangent matrix contribution due to displacement
 * @param[out] Kvis_v Viscous tangent matrix contribution due to velocity
 */
template <size_t nsd>
void get_visc_stress_and_tangent(const dmnType& lDmn, const int eNoN, const Array<double>& Nx, Array<double>& vx_Array, Array<double>& F_Array,
                                 Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v) {

    // Convert vx_Array and F_Array to C-style arrays
    double vx[nsd][nsd], F[nsd][nsd];
    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            vx[i][j] = vx_Array(i,j);
            F[i][j] = F_Array(i,j);
        }
    }
    switch (lDmn.solid_visc.viscType) {
      case consts::SolidViscosityModelType::viscType_Newtonian:
        get_visc_stress_newt<nsd>(lDmn.solid_visc.mu, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);
      break;

      case consts::SolidViscosityModelType::viscType_Potential:
        get_visc_stress_pot<nsd>(lDmn.solid_visc.mu, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);
      break;
    }
}

};

#endif

