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

// Reproduces subroutines defined in MATMODELS.f.

#include "mat_models.h"

#include "fft.h"
#include "mat_fun.h"
#include "mat_fun_carray.h"
#include "utils.h"

#include <math.h>

namespace mat_models {

/// @brief Compute active component of deformation gradient tensor for
/// electromechanics coupling based on active strain formulation
///
/// Reproduces Fortran 'ACTVSTRAIN'.
//
void actv_strain(const ComMod& com_mod, const CepMod& cep_mod, const double gf, 
    const int nfd, const Array<double>& fl, Array<double>& Fa) 
{
  using namespace mat_fun;

  int nsd = com_mod.nsd;
  auto af = fl.col(0);
  auto as = fl.col(1);
  auto an = utils::cross(fl);

  double gn = 4.0 * gf;
  double gs = 1.0 / ((1.0+gf) * (1.0+gn)) - 1.0;

  auto IDm = mat_id(nsd);
  auto Hf = mat_dyad_prod(af, af, nsd);
  auto Hs = mat_dyad_prod(as, as, nsd);
  auto Hn = mat_dyad_prod(an, an, nsd);

  Fa = IDm + gf*Hf + gs*Hs + gn*Hn;
}

void cc_to_voigt(const int nsd, const Tensor4<double>& CC, Array<double>& Dm)
{
  if (nsd == 3) {
    Dm(0,0) = CC(0,0,0,0);
    Dm(0,1) = CC(0,0,1,1);
    Dm(0,2) = CC(0,0,2,2);
    Dm(0,3) = CC(0,0,0,1);
    Dm(0,4) = CC(0,0,1,2);
    Dm(0,5) = CC(0,0,2,0);

    Dm(1,1) = CC(1,1,1,1);
    Dm(1,2) = CC(1,1,2,2);
    Dm(1,3) = CC(1,1,0,1);
    Dm(1,4) = CC(1,1,1,2);
    Dm(1,5) = CC(1,1,2,0);

    Dm(2,2) = CC(2,2,2,2);
    Dm(2,3) = CC(2,2,0,1);
    Dm(2,4) = CC(2,2,1,2);
    Dm(2,5) = CC(2,2,2,0);

    Dm(3,3) = CC(0,1,0,1);
    Dm(3,4) = CC(0,1,1,2);
    Dm(3,5) = CC(0,1,2,0);

    Dm(4,4) = CC(1,2,1,2);
    Dm(4,5) = CC(1,2,2,0);

    Dm(5,5) = CC(2,0,2,0);

    for (int i = 1; i < 6; i++) {
      for (int j = 0; j <= i-1; j++) {
        Dm(i,j) = Dm(j,i);
      }
    }

  } else if (nsd == 2) { 
     Dm(0,0) = CC(0,0,0,0);
     Dm(0,1) = CC(0,0,1,1);
     Dm(0,2) = CC(0,0,0,1);

     Dm(1,1) = CC(1,1,1,1);
     Dm(1,2) = CC(1,1,0,1);

     Dm(2,2) = CC(0,1,0,1);

     Dm(1,0) = Dm(0,1);
     Dm(2,0) = Dm(0,2);
     Dm(2,1) = Dm(1,2);
  } 
}

template <int nsd>
void cc_to_voigt_eigen(const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& CC, Eigen::Matrix<double, 2*nsd, 2*nsd>& Dm)
{
  if (nsd == 3) {
    Dm(0,0) = CC(0,0,0,0);
    Dm(0,1) = CC(0,0,1,1);
    Dm(0,2) = CC(0,0,2,2);
    Dm(0,3) = CC(0,0,0,1);
    Dm(0,4) = CC(0,0,1,2);
    Dm(0,5) = CC(0,0,2,0);

    Dm(1,1) = CC(1,1,1,1);
    Dm(1,2) = CC(1,1,2,2);
    Dm(1,3) = CC(1,1,0,1);
    Dm(1,4) = CC(1,1,1,2);
    Dm(1,5) = CC(1,1,2,0);

    Dm(2,2) = CC(2,2,2,2);
    Dm(2,3) = CC(2,2,0,1);
    Dm(2,4) = CC(2,2,1,2);
    Dm(2,5) = CC(2,2,2,0);

    Dm(3,3) = CC(0,1,0,1);
    Dm(3,4) = CC(0,1,1,2);
    Dm(3,5) = CC(0,1,2,0);

    Dm(4,4) = CC(1,2,1,2);
    Dm(4,5) = CC(1,2,2,0);

    Dm(5,5) = CC(2,0,2,0);

    for (int i = 1; i < 6; i++) {
      for (int j = 0; j <= i-1; j++) {
        Dm(i,j) = Dm(j,i);
      }
    }

  } else if (nsd == 2) {
    Dm(0,0) = CC(0,0,0,0);
    Dm(0,1) = CC(0,0,1,1);
    Dm(0,2) = CC(0,0,0,1);

    Dm(1,1) = CC(1,1,1,1);
    Dm(1,2) = CC(1,1,0,1);

    Dm(2,2) = CC(0,1,0,1);

    Dm(1,0) = Dm(0,1);
    Dm(2,0) = Dm(0,2);
    Dm(2,1) = Dm(1,2);
  }
}

void voigt_to_cc(const int nsd, const Array<double>& Dm, Tensor4<double>& CC)
{
  if (nsd == 3) {
    // Initialize the CC array with zeros
    CC = 0.0;

    // Voigt indices mapping
    const int index_map[6][2] = {
      {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {2, 0}
    };

    // Fill in the CC array based on the Dm matrix
    for (int I = 0; I < 6; ++I) {
      for (int J = 0; J < 6; ++J) {
        int i = index_map[I][0], j = index_map[I][1];
        int k = index_map[J][0], l = index_map[J][1];
        CC(i,j,k,l) = Dm(I,J);
        CC(j,i,k,l) = Dm(I,J);
        CC(i,j,l,k) = Dm(I,J);
        CC(j,i,l,k) = Dm(I,J);
      }
    }

  } else if (nsd == 2) { 
    // Initialize the CC array with zeros
    CC = 0.0;

    // Voigt indices mapping for 2D
    const int index_map[3][2] = {
      {0, 0}, {1, 1}, {0, 1}
    };

    // Fill in the CC array based on the Dm matrix
    for (int I = 0; I < 3; ++I) {
      for (int J = 0; J < 3; ++J) {
        int i = index_map[I][0], j = index_map[I][1];
        int k = index_map[J][0], l = index_map[J][1];
        CC(i,j,k,l) = Dm(I,J);
        CC(i,j,l,k) = Dm(I,J);
        CC(j,i,k,l) = Dm(I,J);
        CC(j,i,l,k) = Dm(I,J);
      }
    }
  }
}




/// @brief Compute additional fiber-reinforcement stress.
///
/// Reproduces Fortran 'GETFIBSTRESS' subroutine.
//
void get_fib_stress(const ComMod& com_mod, const CepMod& cep_mod, const fibStrsType& Tfl, double& g)
{
  using namespace consts;

  g = 0.0;

  if (utils::btest(Tfl.fType, iBC_std)) {
    g = Tfl.g;
  } else if (utils::btest(Tfl.fType, iBC_ustd)) { 
    Vector<double> gv(1), tv(1);
    ifft(com_mod, Tfl.gt, gv, tv);
    g = gv[0];
  }
}




/**
 * @brief Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
 * including both dilational and isochoric components.
 *
 * Reproduces the Fortran 'GETPK2CC' subroutine.
 *
 * @param[in] com_mod Object containing global common variables.
 * @param[in] cep_mod Object containing electrophysiology-specific common variables.
 * @param[in] lDmn Domain object.
 * @param[in] F Deformation gradient tensor.
 * @param[in] nfd Number of fiber directions.
 * @param[in] fl Fiber directions.
 * @param[in] ya Electrophysiology active stress.
 * @param[out] S 2nd Piola-Kirchhoff stress tensor (modified in place).
 * @param[out] Dm Material stiffness tensor (modified in place).
 * @param[out] Ja Jacobian for active strain
 * @return None, but modifies S, Dm, and Ja in place.
 */
template<size_t nsd>
void _get_pk2cc(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const Eigen::Matrix<double, nsd, nsd>& F, const int nfd,
    const Eigen::Matrix<double, nsd, Eigen::Dynamic> fl, const double ya, Eigen::Matrix<double, nsd, nsd>& S, Eigen::Matrix<double, 2*nsd, 2*nsd>& Dm, double& Ja)
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_get_pk2cc
  #ifdef debug_get_pk2cc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Check that this function is only called for struct, ustruct, or fsi physics
  if (lDmn.phys != EquationType::phys_struct && lDmn.phys != EquationType::phys_ustruct && lDmn.phys != EquationType::phys_FSI) {
    throw std::runtime_error("[get_pk2cc] This function is only valid for struct, ustruct, or fsi physics.");
  }

  // ustruct flag
  bool ustruct = (lDmn.phys == EquationType::phys_ustruct);

  S.setZero();
  Dm.setZero();

  // Some preliminaries
  const auto& stM = lDmn.stM;
  double nd = static_cast<double>(nsd);
  double Kp = stM.Kpen;

  // Fiber-reinforced stress
  double Tfa = 0.0;
  get_fib_stress(com_mod, cep_mod, stM.Tf, Tfa);
  double Tsa = Tfa*stM.Tf.eta_s;

  // Electromechanics coupling - active stress
  if (cep_mod.cem.aStress) {
    Tfa = Tfa + ya;
  }

  // Electromechanics coupling - active strain
  auto Fe  = F;
  auto Fa = Eigen::Matrix<double, nsd, nsd>::Identity();
  auto Fai = Fa;

  // if (cep_mod.cem.aStrain) {
  //   actv_strain(com_mod, cep_mod, ya, nfd, fl, Fa);
  //   Fai = Fa.inverse();
  //   Fe = F * Fai;
  // }

  Ja = Fa.determinant();
  double J = Fe.determinant();
  double J2d = pow(J, (-2.0/nd));
  double J4d = J2d*J2d;

  auto Idm = Eigen::Matrix<double, nsd, nsd>::Identity();
  auto C = Fe.transpose() * Fe;
  auto E = 0.50 * (C - Idm);

  auto Ci = C.inverse();
  double trE = E.trace();
  double Inv1 = J2d * C.trace();
  double Inv2 = 0.50 * (Inv1*Inv1 - J4d * (C*C).trace());

  // Contribution of dilational penalty terms to S and CC
  double p  = 0.0;
  double pl = 0.0;

  if (!ustruct) {
    if (!utils::is_zero(Kp)) {
      get_svol_p(com_mod, cep_mod, stM, J, p, pl);
    }
  }
  

  // Now, compute isochoric and total stress, elasticity tensors
  //
  Eigen::TensorFixedSize<double, Eigen::Sizes<nsd,nsd,nsd,nsd>> CC;

  switch (stM.isoType) {

    // NeoHookean model
    case ConstitutiveModelType::stIso_nHook: {
      double g1 = 2.0 * stM.C10;
      Eigen::Matrix<double, nsd, nsd> Sb = g1*Idm;

      // Fiber reinforcement/active stress
      Sb += Tfa * (fl.col(0) * fl.col(0).transpose());

      double r1 = g1 * Inv1 / nd;
      S = J2d*Sb - r1*Ci;

      CC = (-2.0/nd) * ( ten_dyad_prod_eigen<nsd>(Ci, S)  + ten_dyad_prod_eigen<nsd>(S, Ci));
      S += p*J*Ci;
      CC += 2.0*(r1 - p*J) * ten_symm_prod_eigen<nsd>(Ci, Ci)  +  (pl*J - 2.0*r1/nd) * ten_dyad_prod_eigen<nsd>(Ci, Ci);

    } break;

      default:
      throw std::runtime_error("Undefined material constitutive model.");
  } 

  // Convert to Voigt Notation
  cc_to_voigt_eigen<nsd>(CC, Dm);
}

/**
 * @brief Get the 2nd Piola-Kirchhoff stress tensor and material elasticity tensor.
 * 
 * This is a wrapper function for the templated function _get_pk2cc.
 * 
 */
void get_pk2cc(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const Array<double>& F, const int nfd,
    const Array<double>& fl, const double ya, Array<double>& S, Array<double>& Dm, double& Ja)
{
    // Number of spatial dimensions
    int nsd = com_mod.nsd;

    if (nsd == 2) {
        // Copy deformation gradient to Eigen matrix
        auto F_2D = mat_fun::toEigenMatrix<Eigen::Matrix2d>(F);
        
        // Copy fiber directions to Eigen matrix
        Eigen::Matrix<double, 2, Eigen::Dynamic> fl_2D(2, nfd);
        for (int i = 0; i < nfd; i++) {
            fl_2D(0, i) = fl(0, i);
            fl_2D(1, i) = fl(1, i);
        }

        // Initialize stress and elasticity tensors
        Eigen::Matrix2d S_2D = Eigen::Matrix2d::Zero();
        Eigen::Matrix4d Dm_2D = Eigen::Matrix4d::Zero();

        // Call templated function
        _get_pk2cc<2>(com_mod, cep_mod, lDmn, F_2D, nfd, fl_2D, ya, S_2D, Dm_2D, Ja);

        // Copy results back
        mat_fun::toArray(S_2D, S);
        mat_fun::copyDm(Dm_2D, Dm, 4, 4);

    } else if (nsd == 3) {
        // Copy deformation gradient to Eigen matrix
        auto F_3D = mat_fun::toEigenMatrix<Eigen::Matrix3d>(F);

        // Copy fiber directions to Eigen matrix
        Eigen::Matrix<double, 3, Eigen::Dynamic> fl_3D(3, nfd);
        for (int i = 0; i < nfd; i++) {
            fl_3D(0, i) = fl(0, i);
            fl_3D(1, i) = fl(1, i);
            fl_3D(2, i) = fl(2, i);
        }

        // Initialize stress and elasticity tensors
        Eigen::Matrix3d S_3D = Eigen::Matrix3d::Zero();
        Eigen::Matrix<double, 6, 6> Dm_3D;
        Dm_3D.setZero();

        // Call templated function
        _get_pk2cc<3>(com_mod, cep_mod, lDmn, F_3D, nfd, fl_3D, ya, S_3D, Dm_3D, Ja);

        // Copy results back
        mat_fun::toArray(S_3D, S);
        mat_fun::copyDm(Dm_3D, Dm, 6, 6);
    }
}

/// @brief Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
/// for compressible shell elements.
//
void get_pk2cc_shlc(const ComMod& com_mod, const dmnType& lDmn, const int nfd, const Array<double>& fNa0,
    const Array<double>& gg_0, const Array<double>& gg_x, double& g33, Vector<double>& Sml, Array<double>& Dml)
{
  // [NOTE] The tolerance here is a bit larger than Fortran.
  const double ATOL = 1.0e-9;
  //const double ATOL = 1E-10;
  const int MAXITR = 20;

  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_get_pk2cc_shlc
  #ifdef debug_get_pk2cc_shlc 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int nsd = com_mod.nsd;
  Sml = 0.0;
  Dml = 0.0;

  // Initialize tensor operations
  ten_init(3);

  // Some preliminaries
  auto stM = lDmn.stM;
  auto kap = stM.Kpen;
  auto mu = 2.0 * stM.C10;
  auto f13 = 1.0 / 3.0;
  auto f23 = 2.0 / 3.0;
  #ifdef debug_get_pk2cc_shlc 
  dmsg << "kap: " << kap;
  dmsg << "mu: " << mu;
  #endif

  // Inverse of metric coefficients in shell continuum
  auto gi_x = mat_inv(gg_x, 2);

  Array<double> gi_0(3,3);
  auto gg_0_inv = mat_inv(gg_0, 2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      gi_0(i,j) = gg_0_inv(i,j);
    }
  }

  gi_0(2,2) = 1.0;

  // Ratio of inplane Jacobian determinant squared
  auto Jg2 = mat_det(gg_x, 2) / mat_det(gg_0, 2);

  #ifdef debug_get_pk2cc_shlc 
  dmsg << "gi_0: " << gi_0; 
  dmsg << "gi_x: " << gi_x; 
  dmsg << "Jg2: " << Jg2; 
  #endif

  // Begin Newton iterations to satisfy plane-stress condition.
  // The objective is to find C33 that satisfies S33 = 0.
  //
  int itr = 0;
  double C33 = 1.0;
  Array<double> S(3,3);
  Tensor4<double> CC(3,3,3,3);

  while (true) { 
    itr  = itr + 1;
    //dmsg << "------- itr: " << itr;

    // Trace (C)
    auto trC3 = (gg_x(0,0)*gi_0(0,0) + gg_x(0,1)*gi_0(0,1) +  gg_x(1,0)*gi_0(1,0) + 
                 gg_x(1,1)*gi_0(1,1) + C33)*f13;

    // Jacobian-related quantities
    auto J2 = Jg2*C33;
    auto J23 = pow(J2,-f13);
    //dmsg << "trC3: " << trC3;
    //dmsg << "J2: " << J2;
    //dmsg << "J23: " << J23;
    //dmsg << "f13: " << f13;

    // Inverse of curvilinear Cauchy-Green deformation tensor
    //
    Array<double> Ci(3,3);
    Ci(2,2) = 1.0 / C33;

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        Ci(i,j) = gi_x(i,j);
      }
    }
    //dmsg << "Ci: " << Ci;

    // Contribution from dilational penalty terms to S and CC
    auto pJ  = 0.50 * kap * (J2 - 1.0);
    auto plJ = kap * J2;
    #ifdef debug_get_pk2cc_shlc 
    dmsg << "pJ: " << pJ;
    dmsg << "plJ: " << plJ;
    dmsg << "J23: " << J23;
    dmsg << "mu: " << mu;
    dmsg << "trC3: " << trC3;
    dmsg << "stM.isoType: " << stM.isoType;
    #endif

    switch (stM.isoType) {

      case ConstitutiveModelType::stIso_nHook: {
        // 2nd Piola Kirchhoff stress
        S = mu*J23*(gi_0 - trC3*Ci) + pJ*Ci;

        // Elasticity tensor
        CC = (mu*J23*f23*trC3 + plJ)*ten_dyad_prod(Ci, Ci, 3) + (mu*J23*trC3 - pJ)*2.0*ten_symm_prod(Ci, Ci, 3) - 
            f23*mu*J23*(ten_dyad_prod(gi_0, Ci, 3) + ten_dyad_prod(Ci, gi_0, 3));
      } break;

      case ConstitutiveModelType::stIso_MR: {

        // 2nd Piola Kirchhoff stress
        auto C1  = stM.C10;
        auto C2  = stM.C01;
        auto J43 = pow(J2,-f23);

        auto I2ijkl = ten_dyad_prod(gi_0, gi_0, 3) - ten_symm_prod(gi_0, gi_0, 3);
        auto I2ij = ten_mddot(I2ijkl, Ci, 3);

        auto Gi4AS = ten_asym_prod12(gi_0, gi_0, 3);
        auto I2 = mat_ddot(Ci, ten_mddot(Gi4AS, Ci, 3), 3);
        auto Cikl = -1.0 * ten_symm_prod(Ci, Ci, 3);

        S  = C1*J23*(gi_0 - trC3*Ci) + pJ*Ci + C2*J43*(I2ij - f23*I2*Ci);

        //  Elasticity tensor
        CC = (C1*J23*f23*trC3 + plJ)*ten_dyad_prod(Ci, Ci, 3) + (C1*J23*trC3 - pJ)*2.0*ten_symm_prod(Ci, Ci, 3) - 
            f23*C1*J23*(ten_dyad_prod(gi_0, Ci, 3) + ten_dyad_prod(Ci, gi_0, 3));

        CC += 2.0 * f23 * C2 * J43 * ( ten_dyad_prod((f23*I2*Ci-I2ij), Ci, 3) - I2*Cikl - 
            ten_dyad_prod(Ci, I2ij, 3) + I2ijkl);
      } break;

      case ConstitutiveModelType::stIso_HO_ma: {

         if (nfd !=  2) {
           throw std::runtime_error("Min fiber directions not defined for Holzapfel material model (1)");
         }
     
         Array3<double> fl(2,2,nfd);

         for (int iFn = 0; iFn < nfd; iFn++) {
           fl(0,0,iFn) = fNa0(0,iFn)*fNa0(0,iFn);
           fl(0,1,iFn) = fNa0(0,iFn)*fNa0(1,iFn);
           fl(1,0,iFn) = fNa0(1,iFn)*fNa0(0,iFn);
           fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn);
         }

         //  Compute fiber-based invariants
         //
         auto Inv4 = gg_x(0,0)*fl(0,0,0) + gg_x(0,1)*fl(0,1,0) + gg_x(1,0)*fl(1,0,0) + gg_x(1,1)*fl(1,1,0);
         auto Inv6 = gg_x(0,0)*fl(0,0,1) + gg_x(0,1)*fl(0,1,1) + gg_x(1,0)*fl(1,0,1) + gg_x(1,1)*fl(1,1,1);
         auto Inv8 = gg_x(0,0)*fNa0(0,0)*fNa0(0,1) + gg_x(0,1)*fNa0(0,0)*fNa0(1,1) + 
             gg_x(1,0)*fNa0(1,0)*fNa0(0,1) + gg_x(1,1)*fNa0(1,0)*fNa0(1,1);

         auto Eff = Inv4 - 1.0;
         auto Ess = Inv6 - 1.0;
         auto Efs = Inv8;

         // Smoothed heaviside function
         auto c4f = 1.0 / (1.0 + exp(-stM.khs*Eff));
         auto c4s = 1.0 / (1.0 + exp(-stM.khs*Ess));

         // Approx. derivative of smoothed heaviside function
         auto dc4f = 0.250 * stM.khs * exp(-stM.khs*fabs(Eff));
         auto dc4s = 0.250 * stM.khs * exp(-stM.khs*fabs(Ess));

         // Add isochoric stress and stiffness contribution
         //
         // EI1  = I1 + Jg2i - 3.0
         //
         auto d1 = stM.a*J23*exp(2.0*stM.b*(trC3*J23 - 1.0));
         auto SN = (gi_0 - trC3*Ci);

         S  = d1*SN + pJ*Ci;

         CC = f23*trC3*ten_dyad_prod(Ci, Ci, 3)
            + trC3*2.0*ten_symm_prod(Ci, Ci, 3)
            - f23*(ten_dyad_prod(gi_0, Ci, 3)
            + ten_dyad_prod(Ci, gi_0, 3))
            + 2.0*stM.b*J23*ten_dyad_prod(SN, SN, 3);

         CC = d1*CC + plJ*ten_dyad_prod(Ci, Ci, 3) - pJ*2.0*ten_symm_prod(Ci, Ci, 3);

         // Anisotropic part
         // Fiber sheet
         //
         Array<double> Hfs(3,3);
         auto hfs_sym_prod = mat_symm_prod(fNa0.col(0), fNa0.col(1), 2);

         for (int i = 0; i < 2; i++) {
           for (int j = 0; j < 2; j++) {
             Hfs(i,j) = hfs_sym_prod(i,j);;
           }
         }

         auto g1 = 2.0*stM.afs*exp(stM.bfs*Efs*Efs);
         S += g1*Efs*Hfs;

         auto g2 = g1*2.0*(1.0 + 2.0*Efs*Efs);
         CC += g2*ten_dyad_prod(Hfs, Hfs, 3);

         // Fiber
         Array<double> flM(3,3);

         if (Eff > 0.0) {
           for (int i = 0; i < 2; i++) {
             for (int j = 0; j < 2; j++) {
               flM(i,j) = fl(i,j,0);
             }
           }

           // S  = S + 2.0*stM.aff*Eff*flM
           // CC = CC + 4.0*stM.aff*ten_dyad_prod(flM, flM, 2)
           g1 = 2.0*stM.aff*exp(stM.bff*Eff*Eff);
           S += g1 * Eff * flM;
           g2 = g1 * 2.0 * (1.0  +  2.0 * stM.bff * Eff * Eff);
           CC += g2 * ten_dyad_prod(flM, flM, 3);
         }

         // Sheet
         //
         if (Ess > 0.0) {
           Array<double> flM(3,3);
           auto flM_1 = fl.rslice(1);
           for (int i = 0; i < 2; i++) {
             for (int j = 0; j < 2; j++) {
               flM(i,j) = flM_1(i,j); 
             }
           }

           auto g1 = 2.0 * stM.ass * exp(stM.bss*Ess*Ess);
           S += g1*Ess*flM;
           auto g2 = g1 * 2.0 * (1.0  +  2.0 * stM.bss * Ess * Ess);
           CC += g2*ten_dyad_prod(flM, flM, 3);
         }
      } break;

      case ConstitutiveModelType::stIso_LS: {
        Array3<double> fl(2,2,nfd);

         for (int iFn = 0; iFn < nfd; iFn++) {
           fl(0,0,iFn) = fNa0(0,iFn)*fNa0(0,iFn);
           fl(0,1,iFn) = fNa0(0,iFn)*fNa0(1,iFn);
           fl(1,0,iFn) = fNa0(1,iFn)*fNa0(0,iFn);
           fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn);
         }

         //  Compute fiber-based invariants

         auto Inv4 = gg_x(0,0)*fl(0,0,0) + gg_x(0,1)*fl(0,1,0) + gg_x(1,0)*fl(1,0,0) + gg_x(1,1)*fl(1,1,0);
         auto Eff  = Inv4 - 1.0;

         // Isotropic contribution
         //
         auto d1 = 2.0*stM.a0*stM.mu0*stM.b1*J23 * exp(2.0*stM.b1*(trC3*J23 - 1.0));
         auto SN = (gi_0 - trC3*Ci);
         auto CCN = f23*trC3*ten_dyad_prod(Ci, Ci, 3) + trC3*2.0*ten_symm_prod(Ci, Ci, 3) - 
             f23*(ten_dyad_prod(gi_0, Ci, 3) + ten_dyad_prod(Ci, gi_0, 3));

         S  = (stM.a + d1*(2.0*(trC3*J23 - 1.0)))*SN + pJ*Ci;

         CC = (1.0 + 18.0*stM.b1*(trC3*J23 - 1.0) * (trC3*J23 - 1.0))*J23*ten_dyad_prod(SN, SN, 3) + 
             3.0*(trC3*J23 - 1.0)*CCN;
         CC = stM.a*CCN + 2.0*d1*CC + plJ*ten_dyad_prod(Ci, Ci, 3) - pJ*2.0*ten_symm_prod(Ci, Ci, 3);

         // Anisotropic contribution
         Array<double> flM(3,3);

         if (Eff > 0.0) {
           auto flM_0 = fl.rslice(0);
           for (int i = 0; i < 2; i++) {
             for (int j = 0; j < 2; j++) {
               flM(i,j) = flM_0(i,j);
             }
           }

           S += 2.0*stM.a0*(1.0 - stM.mu0) *stM.b2*Eff * exp(stM.b2*Eff*Eff)*flM;

           CC += 4.0*stM.a0*(1.0 - stM.mu0) * stM.b2 * (1.0 + 2.0*stM.b2*Eff*Eff) * 
               exp(stM.b2*Eff*Eff) * ten_dyad_prod(flM, flM, 3);
         }

      } break;

      default:
        //err = "Undefined material constitutive model"
        break;

     }

     if (fabs(S(2,2)) <= ATOL) {
       break;
     }

     if (itr > MAXITR) {
        std::cerr << "[get_pk2cc_shlc] Failed to converge plane-stress condition." << std::endl;
        //exit(0);
        break;
     }

     C33 = C33 - (2.0 * S(2,2) / CC(2,2,2,2));
     //dmsg << "1: C33: " << C33;
     //dmsg << "CC(3,3,3,3): " << CC(2,2,2,2);
     //exit(0);
  }

  g33 = C33;

  // Statically condense CC
  //
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < 2; l++) { 
          C33 = CC(i,j,2,2)*CC(2,2,k,l) / CC(2,2,2,2);
          CC(i,j,k,l) = CC(i,j,k,l) - C33;
        }
      }
    }
  }

  g33 = C33;

  //dmsg << "2: C33: " << C33;
  //exit(0);

  // Convert the in-plane components to Voigt notation
  Sml(0) = S(0,0);
  Sml(1) = S(1,1);
  Sml(2) = S(0,1);

  Dml(0,0) = CC(0,0,0,0);
  Dml(0,1) = CC(0,0,1,1);
  Dml(0,2) = CC(0,0,0,1);

  Dml(1,1) = CC(1,1,1,1);
  Dml(1,2) = CC(1,1,0,1);

  Dml(2,2) = CC(0,1,0,1);

  Dml(1,0) = Dml(0,1);
  Dml(2,0) = Dml(0,2);
  Dml(2,1) = Dml(1,2);
}

/// @brief Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
/// for incompressible shell elements
///
/// Reproduces Fortran GETPK2CC_SHLi
//
void get_pk2cc_shli(const ComMod& com_mod, const dmnType& lDmn, const int nfd, const Array<double>& fNa0, 
    const Array<double>& gg_0, const Array<double>& gg_x, double& g33, Vector<double>& Sml, Array<double>& Dml)
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_get_pk2cc_shli
  #ifdef debug_get_pk2cc_shli 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int nsd = com_mod.nsd;
  Sml = 0.0;
  Dml = 0.0;

  // Some preliminaries
  auto stM  = lDmn.stM;

  // Inverse of metric coefficients in shell continuum
  auto gi_0 = mat_inv(gg_0,2);
  auto gi_x = mat_inv(gg_x,2);

  // Ratio of inplane Jacobian determinants
  auto Jg2i = mat_det(gg_x,2);

  if (is_zero(Jg2i)) {
    throw std::runtime_error(" Divide by zero in-plane Jacobian determinant.");
  }

  Jg2i = mat_det(gg_0,2) / Jg2i;

  double I1 = 0.0;

  for (int a = 0; a < 2; a++) {
    for (int b = 0; b < 2; b++) {
      I1 = I1 + gi_0(a,b) * gg_x(a,b);
    }
  }

  Tensor4<double> CC(nsd,nsd,nsd,nsd);
  Array<double> S;
  Array3<double> fl(2,2,nfd);

  switch (stM.isoType) {

    case ConstitutiveModelType::stIso_nHook: {
      auto mu = 2.0 * stM.C10;
      S = mu*(gi_0 - Jg2i*gi_x);
      auto CC = 2.0*mu*Jg2i*(ten_dyad_prod(gi_x, gi_x,1) + ten_symm_prod(gi_x, gi_x,1));
    } break;

    case ConstitutiveModelType::stIso_MR: {
      auto SN = (gi_0 - Jg2i*gi_x);
      auto CCN = 2.0*Jg2i*(ten_dyad_prod(gi_x, gi_x,1) + ten_symm_prod(gi_x, gi_x,1));
      S = stM.C10*SN + stM.C01*Jg2i* (gi_0 - I1*gi_x) + stM.C01/Jg2i*gi_x;

      CC  = (stM.C10 + stM.C01*I1) * CCN - 2.0*stM.C01 * Jg2i * (ten_dyad_prod(gi_0, gi_x,1) +
            ten_dyad_prod(gi_x, gi_0,1)) + 2.0*stM.C01 / Jg2i *(ten_dyad_prod(gi_x, gi_0,1) -
            ten_symm_prod(gi_x, gi_x,1));
    } break;

    // HO (Holzapfel-Ogden) model for myocardium with full invariants
    // for the anisotropy terms (modified-anisotropy)
    //
    case ConstitutiveModelType::stIso_HO_ma: {

      if (nfd != 2) {
        throw std::runtime_error("Min fiber directions not defined for Holzapfel material model (1)");
      }

      for (int iFn = 0; iFn < nfd; iFn++) {
        fl(0,0,iFn) = fNa0(0,iFn)*fNa0(0,iFn);
        fl(0,1,iFn) = fNa0(0,iFn)*fNa0(1,iFn);
        fl(1,0,iFn) = fNa0(1,iFn)*fNa0(0,iFn);
        fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn);
      }

      // Compute fiber-based invariants
      auto Inv4 = gg_x(0,0)*fl(0,0,0) + gg_x(0,1)*fl(0,1,0) + gg_x(1,0)*fl(1,0,0) + gg_x(1,1)*fl(1,1,0);
      auto Inv6 = gg_x(0,0)*fl(0,0,1) + gg_x(0,1)*fl(0,1,1) + gg_x(1,0)*fl(1,0,1) + gg_x(1,1)*fl(1,1,1);
      auto Inv8 = gg_x(0,0)*fNa0(0,0)*fNa0(0,1) + gg_x(0,1)*fNa0(0,0)*fNa0(1,1) + 
          gg_x(1,0)*fNa0(1,0)*fNa0(0,1) + gg_x(1,1)*fNa0(1,0)*fNa0(1,1);

      auto Eff = Inv4 - 1.0;
      auto Ess = Inv6 - 1.0;
      auto Efs = Inv8;

      // Smoothed heaviside function
      auto c4f = 1.0 / (1.0 + exp(-stM.khs*Eff));
      auto c4s = 1.0 / (1.0 + exp(-stM.khs*Ess));

      // Approx. derivative of smoothed heaviside function
      auto dc4f = 0.250*stM.khs*exp(-stM.khs*fabs(Eff));
      auto dc4s = 0.250*stM.khs*exp(-stM.khs*fabs(Ess));

      // Add isochoric stress and stiffness contribution
      //
      auto EI1 = I1 + Jg2i - 3.0;
      auto SN = (gi_0 - Jg2i*gi_x);
      auto CCN = 2.0*Jg2i*(ten_dyad_prod(gi_x, gi_x,1) + ten_symm_prod(gi_x, gi_x,1));

      auto d1 = stM.a*exp(stM.b*EI1);
      S = d1*SN;
      CC = d1*(CCN + 2.0*stM.b*ten_dyad_prod(SN, SN,1));

      // Anisotropic part
      // Fiber sheet
      auto Hfs = mat_symm_prod(fNa0.col(0), fNa0.col(0),1);
      auto g1 = 2.0*stM.afs*exp(stM.bfs*Efs*Efs);
      S += g1*Efs*Hfs;

      auto g2 = g1*2.0*(1.0 + 2.0*Efs*Efs);
      CC += CC + g2*ten_dyad_prod(Hfs, Hfs,1);

      // Fiber
      if (Eff > 0.0) {
        auto flM = fl.slice(0);
        // S  = S + 2.0*stM.aff*Eff*flM
        // CC = CC + 4.0*stM.aff*ten_dyad_prod(flM, flM,1)
        g1 = 2.0*stM.aff*exp(stM.bff*Eff*Eff);
        S += g1*Eff*flM;
        g2 = g1*2.0*(1.0 + 2.0*stM.bff*Eff*Eff);
        CC += g2*ten_dyad_prod(flM, flM,1);
      }

      // Sheet
      if (Ess > 0.0) {
        auto flM = fl.slice(0);
        g1 = 2.0*stM.ass*exp(stM.bss*Ess*Ess);
        S += g1*Ess*flM;
        g2 = g1*2.0*(1.0 + 2.0*stM.bss*Ess*Ess);
        CC += g2*ten_dyad_prod(flM, flM,1);
      }
    } break;

    // Lee Sacks model for aorta with full invariants
    // for the anisotropy terms (modified-anisotropy)
    //
    case ConstitutiveModelType::stIso_LS: {
      auto SN = (gi_0 - Jg2i*gi_x);
      auto CCN = 2.0*Jg2i*(ten_dyad_prod(gi_x, gi_x,1) + ten_symm_prod(gi_x, gi_x,1));

      for (int iFn = 0; iFn < nfd; iFn++) {
        fl(0,0,iFn) = fNa0(0,iFn)*fNa0(0,iFn);
        fl(0,1,iFn) = fNa0(0,iFn)*fNa0(1,iFn);
        fl(1,0,iFn) = fNa0(1,iFn)*fNa0(0,iFn);
        fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn);
      }

      // Compute fiber-based invariants
      auto Inv4 = gg_x(0,0)*fl(0,0,0) + gg_x(0,1)*fl(0,1,0) + gg_x(1,0)*fl(1,0,0) + gg_x(1,1)*fl(1,1,0);
      auto Eff = Inv4 - 1.0;

      // Isotropic contribution
      auto EI1 = (I1 + Jg2i - 3.0);
      S = stM.a*SN + 2.0*stM.a0*stM.mu0*stM.b1*EI1 * exp(stM.b1*EI1*EI1) * SN;
      CC = stM.a*CCN + 4.0*stM.a0*stM.mu0*stM.b1 * exp(stM.b1*EI1*EI1) * ((1.0 + 2.0*stM.b1*EI1*EI1) * 
          ten_dyad_prod(SN, SN,1) + 0.50*EI1*CCN);

      // Anisotropic contribution
      if (Eff > 0.0) {
        auto flM = fl.rslice(0);
        S += 2.0*stM.a0*(1.0 - stM.mu0)*stM.b2*Eff * exp(stM.b2*Eff*Eff)*flM;
        CC += 4.0*stM.a0*(1.0 - stM.mu0)*stM.b2 * (1.0 + 2.0*stM.b2*Eff*Eff) * 
            exp(stM.b2*Eff*Eff) * ten_dyad_prod(flM, flM,1);
      }
    } break;

    default: 
     //err = "Undefined material constitutive model"
      break;
  }

  g33 = Jg2i;

  // Convert to Voigt notation
  Sml(0) = S(0,0);
  Sml(1) = S(1,1);
  Sml(2) = S(0,1);

  Dml(0,0) = CC(0,0,0,0);
  Dml(0,1) = CC(0,0,1,1);
  Dml(0,2) = CC(0,0,0,1);

  Dml(1,1) = CC(1,1,1,1);
  Dml(1,2) = CC(1,1,0,1);

  Dml(2,2) = CC(0,1,0,1);

  Dml(1,0) = Dml(0,1);
  Dml(2,0) = Dml(0,2);
  Dml(2,1) = Dml(1,2);

}


/// @brief Reproduces Fortran 'GETSVOLP'.
//
void get_svol_p(const ComMod& com_mod, const CepMod& cep_mod, const stModelType& stM, const double J, 
    double& p, double& pl) 
{
  using namespace consts;

  double Kp = stM.Kpen;

  switch (stM.volType) {
    case ConstitutiveModelType::stVol_Quad: 
      p  = Kp*(J-1.0);
      pl = Kp*(2.0*J-1.0);
    break;
    
    case ConstitutiveModelType::stVol_ST91:
      p  = 0.50*Kp*(J-1.0/J);
      pl = Kp*J;
    break;

    case ConstitutiveModelType::stVol_M94:
      p  = Kp*(1.0-1.0/J);
      pl = Kp;
    break;
  } 
}

/// @brief Compute stabilization parameters tauM and tauC.
///
/// Reproduces Fortran 'GETTAU'.
//
void get_tau(const ComMod& com_mod, const dmnType& lDmn, const double detF, const double Je, double& tauM, double& tauC)
{
  using namespace consts;

  double he = 0.50 * pow(Je,1.0/static_cast<double>(com_mod.nsd));
  double rho0 = lDmn.prop.at(PhysicalProperyType::solid_density);
  double Em   = lDmn.prop.at(PhysicalProperyType::elasticity_modulus);
  double nu   = lDmn.prop.at(PhysicalProperyType::poisson_ratio);
  double ctM  = lDmn.prop.at(PhysicalProperyType::ctau_M);
  double ctC  = lDmn.prop.at(PhysicalProperyType::ctau_C);

  double mu = 0.50*Em / (1.0 + nu);
  double c = 0.0;

  if (utils::is_zero(nu-0.50)) {
     c = sqrt(mu / rho0);
  } else { 
     double lam = 2.0*mu*nu / (1.0-2.0*nu);
     c = sqrt((lam + 2.0*mu)/rho0);
  }

  tauM = ctM * (he/c) * (detF/rho0);
  tauC = ctC * (he*c) * (rho0/detF);
}

/**
 * @brief Compute rho, beta, drho/dp, dbeta/dp for volumetric penalty terms in 
 * the ustruct formulation.
 *
 * See ustruct paper (https://doi.org/10.1016/j.cma.2018.03.045) Section 2.4.
 *
 * @param[in] com_mod Object containing global common variables
 * @param[in] lDmn Domain object
 * @param[in] p Pressure.
 * @param[out] ro Solid density, rho.
 * @param[out] bt Isothermal compressibility coefficient, beta.
 * @param[out] dro Derivative of rho with respect to p.
 * @param[out] dbt Derivative of beta with respect to p.
 * @param[out] Ja Active strain Jacobian.
 * @return None, but updates ro, bt, dro, and dbt.
 */
void g_vol_pen(const ComMod& com_mod, const dmnType& lDmn, const double p, 
    double& ro, double& bt, double& dro, double& dbt, const double Ja)
{
  using namespace consts;

  ro = lDmn.prop.at(PhysicalProperyType::solid_density) / Ja;
  bt  = 0.0;
  dbt = 0.0;
  dro = 0.0;

  double Kp = lDmn.stM.Kpen;

  if (utils::is_zero(Kp)) {
    return;
  }

  switch (lDmn.stM.volType) {

    case ConstitutiveModelType::stVol_Quad : {
      double r1  = 1.0/(Kp - p);
      ro  = ro*Kp*r1;
      bt  = r1;
      dro = ro*r1;
      dbt = r1*r1;
    } break;

    case ConstitutiveModelType::stVol_ST91: {
      double r1 = ro/Kp;
      double r2 = sqrt(p*p + Kp*Kp);

      ro  = r1*(p + r2);
      bt  = 1.0/r2;
      dro = ro*bt;
      dbt = -bt*p/(p*p + Kp*Kp);
    } break;

    case ConstitutiveModelType::stVol_M94: {
      double r1  = ro/Kp;
      double r2  = Kp + p;

      ro  = r1*r2;
      bt  = 1.0/r2;
      dro = r1;
      dbt = -bt*bt;
    } break;

    default:
    break;
  }
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
void get_visc_stress_pot(const double mu, const int eNoN, const Array<double>& Nx, const Array<double>& vx, const Array<double>& F,
                        Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v) {

    using namespace consts;
    using namespace mat_fun;
    using namespace utils;

    // Number of spatial dimensions
    int nsd = F.nrows();

    // Initialize Svis, Kvis_u, Kvis_v to zero
    Svis = 0.0;
    Kvis_u = 0.0;
    Kvis_v = 0.0;


    // Required intermediate terms for stress and tangent
    auto Ft = transpose(F);
    auto F_Ft = mat_mul(F, Ft);
    auto Ft_vx = mat_mul(Ft, vx);
    auto vxt = transpose(vx);
    auto F_vxt = mat_mul(F, vxt);

    //double F_Nx[nsd][eNoN] = {0}, vx_Nx[nsd][eNoN] = {0};
    Array<double> F_Nx(nsd,eNoN), vx_Nx(nsd,eNoN);
    
    for (int a = 0; a < eNoN; ++a) {
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                F_Nx(i,a) += F(i,j) * Nx(j,a);
                vx_Nx(i,a) += vx(i,j) * Nx(j,a);
            }
        }
    }

    // 2nd Piola-Kirchhoff stress due to viscosity
    // Svis = mu * 1/2 * ( (F^T * dv/dX) + (F^T * dv/dX)^T )
    Svis = mu * mat_symm(Ft_vx, nsd);

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
                    Kvis_u(ii,a,b) = 0.5 * mu * (F_Nx(i,b) * vx_Nx(j,a) + Nx_Nx * F_vxt(i,j));
                    Kvis_v(ii,a,b) = 0.5 * mu * (Nx_Nx * F_Ft(i,j) + F_Nx(i,b) * F_Nx(j,a));
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
void get_visc_stress_newt(const double mu, const int eNoN, const Array<double>& Nx, const Array<double>& vx, const Array<double>& F,
                           Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v) {
    using namespace consts;
    using namespace mat_fun;
    using namespace utils;

    // Number of spatial dimensions
    int nsd = F.nrows();

    // Initialize Svis, Kvis_u, Kvis_v to zero
    Svis = 0.0;
    Kvis_u = 0.0;
    Kvis_v = 0.0;

    // Get identity matrix, Jacobian, and F^-1
    auto Idm = mat_id(nsd);
    auto J = mat_det(F, nsd);
    auto Fi = mat_inv(F, nsd);

    // Required intermediate terms for stress and tangent
    // vx_Fi: Velocity gradient in current configuration          
    auto vx_Fi = mat_mul(vx, Fi);
    auto vx_Fi_symm = mat_symm(vx_Fi, nsd);
    // ddev: Deviatoric part of rate of strain tensor
    auto ddev = mat_dev(vx_Fi_symm, nsd);
    //double Nx_Fi[nsd][eNoN] = {0}, ddev_Nx_Fi[nsd][eNoN] = {0}, vx_Fi_Nx_Fi[nsd][eNoN] = {0};
    Array<double> Nx_Fi(nsd,eNoN), ddev_Nx_Fi(nsd,eNoN), vx_Fi_Nx_Fi(nsd,eNoN);
    for (int a = 0; a < eNoN; ++a) {
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                Nx_Fi(i,a) += Nx(j,a) * Fi(j,i);
            }
        }
        ddev_Nx_Fi = mat_mul(ddev, Nx_Fi);
        vx_Fi_Nx_Fi = mat_mul(vx_Fi, Nx_Fi);
    }

    // 2nd Piola-Kirchhoff stress due to viscosity
    // Svis = 2 * mu * J * F^-1 * d_dev * F^-T
    auto Fit = transpose(Fi);
    auto ddev_Fit = mat_mul(ddev, Fit);
    auto Fi_ddev_Fit = mat_mul(Fi, ddev_Fit);
    Svis = 2.0 * mu * J * Fi_ddev_Fit;

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
                                    (Nx_Fi_Nx_Fi * vx_Fi(i,j) + Nx_Fi(i,b) * vx_Fi_Nx_Fi(j,a) -
                                    r2d * Nx_Fi(i,a) * vx_Fi_Nx_Fi(j,b)));

                    // Derivative of the residual w.r.t velocity
                    Kvis_v(ii,a,b) = mu * J * (Nx_Fi_Nx_Fi * Idm(i,j) +
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
void get_visc_stress_and_tangent(const dmnType& lDmn, const int eNoN, const Array<double>& Nx, const  Array<double>& vx, const  Array<double>& F,
                                 Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v) {

    switch (lDmn.solid_visc.viscType) {
      case consts::SolidViscosityModelType::viscType_Newtonian:
        get_visc_stress_newt(lDmn.solid_visc.mu, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);
      break;

      case consts::SolidViscosityModelType::viscType_Potential:
        get_visc_stress_pot(lDmn.solid_visc.mu, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);
      break;
    }
}

};
