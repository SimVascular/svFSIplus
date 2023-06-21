
// Reproduces subroutines defined in MATMODELS.f.

#include "mat_models.h"

#include "fft.h"
#include "mat_fun.h"
#include "utils.h"

#include <math.h>

namespace mat_models {

//-------------
// actv_strain
//-------------
// Compute active component of deformation gradient tensor for
// electromechanics coupling based on active strain formulation
//
// Reproduces Fortran 'ACTVSTRAIN'.
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

//-------------
// cc_to_voigt
//-------------
//
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

//----------------
// get_fib_stress
//----------------
// Compute additional fiber-reinforcement stress.
//
// Reproduces Fortran 'GETFIBSTRESS' subroutine.
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

//-----------
// get_pk2cc
//-----------
// Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
// including both dilational and isochoric components.
//
// Reproduces the Fortran 'GETPK2CC' subroutine.
//
void get_pk2cc(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const Array<double>& F, const int nfd,
    const Array<double>& fl, const double ya, Array<double>& S, Array<double>& Dm)
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_get_pk2cc
  #ifdef debug_get_pk2cc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int nsd = com_mod.nsd;
  S = 0.0;
  Dm = 0.0;

  // Some preliminaries
  const auto& stM = lDmn.stM;
  double nd = static_cast<double>(nsd);
  double Kp = stM.Kpen;

  // Fiber-reinforced stress
  double Tfa = 0.0;
  get_fib_stress(com_mod, cep_mod, stM.Tf, Tfa);

  // Electromechanics coupling - active stress
  if (cep_mod.cem.aStress) {
    Tfa = Tfa + ya;
  }

  // Electromechanics coupling - active strain
  auto Fe  = F;
  auto Fa = mat_id(nsd);
  auto Fai = Fa;

  if (cep_mod.cem.aStrain) {
    actv_strain(com_mod, cep_mod, ya, nfd, fl, Fa);
    Fai = mat_inv(Fa, nsd);
    Fe = mat_mul(F, Fai);
  }

  double J = mat_det(Fe, nsd);
  double J2d = pow(J, (-2.0/nd));
  double J4d = J2d*J2d;

  auto Idm = mat_id(nsd);
  auto C = mat_mul(transpose(Fe), Fe);
  Array<double> E(nsd,nsd);
  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < nsd; j++) {
      E(i,j) = 0.50 * (C(i,j) - Idm(i,j));
    }
  }

  auto Ci = mat_inv(C, nsd);
  double trE = mat_trace(E, nsd);
  double Inv1 = J2d * mat_trace(C,nsd);
  double Inv2 = 0.50 * (Inv1*Inv1 - J4d * mat_trace(mat_mul(C,C), nsd));

  // Contribution of dilational penalty terms to S and CC
  double p  = 0.0;
  double pl = 0.0;

  if (!utils::is_zero(Kp)) {
    get_svol_p(com_mod, cep_mod, stM, J, p, pl);
  }

  // Now, compute isochoric and total stress, elasticity tensors
  //
  Tensor4<double> CC(nsd,nsd,nsd,nsd);

  switch (stM.isoType) {
    case ConstitutiveModelType::stIso_lin: {
      double g1 = stM.C10;    // mu
      S = g1*Idm;
      return; 
    } break;

    // St.Venant-Kirchhoff
    case ConstitutiveModelType::stIso_StVK: {
      double g1 = stM.C10;         // lambda
      double g2 = stM.C01 * 2.0;   // 2*mu

      S = g1*trE*Idm + g2*E;
      CC = g1 * ten_dyad_prod(Idm, Idm, nsd) + g2*ten_ids(nsd);
    } break;

    // modified St.Venant-Kirchhoff
    case ConstitutiveModelType::stIso_mStVK: {
      double g1 = stM.C10; // kappa
      double g2 = stM.C01;  // mu

      S = g1*log(J)*Ci + g2*(C-Idm);
      CC = g1 * ( -2.0*log(J)*ten_symm_prod(Ci, Ci, nsd) +
         ten_dyad_prod(Ci, Ci, nsd) ) + 2.0*g2*ten_ids(nsd);
    } break;

    // NeoHookean model
    case ConstitutiveModelType::stIso_nHook: {
      double g1 = 2.0 * stM.C10;
      auto Sb = g1*Idm;

      // Fiber reinforcement/active stress
      Sb += Tfa * mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      double r1 = g1 * Inv1 / nd;
      for (int j = 0; j < S.ncols(); j++) {
        for (int i = 0; i < S.nrows(); i++) {
          S(i,j) = J2d*Sb(i,j) - r1*Ci(i,j);
        }
      }

      CC = (-2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd));
      S += p*J*Ci;
      CC += 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd)  +  (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);

    } break;

    //  Mooney-Rivlin model
    case ConstitutiveModelType::stIso_MR: {
      double g1  = 2.0 * (stM.C10 + Inv1*stM.C01);
      double g2  = -2.0 * stM.C01;
      auto Sb = g1*Idm + g2*J2d*C;

      // Fiber reinforcement/active stress
      Sb = Sb + Tfa*mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      g1  = 4.0*J4d* stM.C01;
      auto CCb = g1 * (ten_dyad_prod(Idm, Idm, nsd) - ten_ids(nsd));

      double r1  = J2d*mat_ddot(C, Sb, nsd) / nd;
      S = J2d*Sb - r1*Ci;

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC = ten_ddot(CCb, PP, nsd);
      CC = ten_transpose(CC, nsd);
      CC = ten_ddot(PP, CC, nsd);
      CC = CC - (2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd) );

      S  = S + p*J*Ci;
      CC = CC + 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd) + (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);
    } break;

    // HGO (Holzapfel-Gasser-Ogden) model with additive splitting of
    // the anisotropic fiber-based strain-energy terms
    case ConstitutiveModelType::stIso_HGO: {
      if (nfd != 2) {
        throw std::runtime_error("[get_pk2cc] Min fiber directions not defined for HGO material model.");
      }
      double kap = stM.kap;
      double Inv4 = J2d*utils::norm(fl.col(0), mat_mul(C, fl.col(0)));
      double Inv6 = J2d*utils::norm(fl.col(1), mat_mul(C, fl.col(1)));

      double Eff = kap*Inv1 + (1.0-3.0*kap)*Inv4 - 1.0;
      double Ess = kap*Inv1 + (1.0-3.0*kap)*Inv6 - 1.0;

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
    } break;

    // Guccione (1995) transversely isotropic model
    case ConstitutiveModelType::stIso_Gucci: {
      if (nfd != 2) {
        throw std::runtime_error("[get_pk2cc] Min fiber directions not defined for Guccione material model.");
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

      // Fiber stiffness contribution := (dE*_ab / dE_IJ)
      Array3<double> RmRm(nsd,nsd,6);

      RmRm.set_slice(0, mat_dyad_prod(Rm.col(0), Rm.col(0), nsd));
      RmRm.set_slice(1, mat_dyad_prod(Rm.col(1), Rm.col(1), nsd));
      RmRm.set_slice(2, mat_dyad_prod(Rm.col(2), Rm.col(2), nsd));

      RmRm.set_slice(3, mat_symm_prod(Rm.col(0), Rm.col(1), nsd));
      RmRm.set_slice(4, mat_symm_prod(Rm.col(1), Rm.col(2), nsd));
      RmRm.set_slice(5, mat_symm_prod(Rm.col(2), Rm.col(0), nsd));

      auto Sb = g1 *  Es(0,0) * RmRm.slice(0) + 
                g2 * (Es(1,1) * RmRm.slice(1) + Es(2,2)*RmRm.slice(2) + 2.0*Es(1,2)*RmRm.slice(4)) +
          2.0 * g3 * (Es(0,1) * RmRm.slice(3) + Es(0,2)*RmRm.slice(5));

      auto CCb = 2.0*ten_dyad_prod(Sb, Sb, nsd);
      Sb = Sb * r2;

      // Fiber reinforcement/active stress
      Sb = Sb + Tfa*mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      double r1 = J2d*mat_ddot(C, Sb, nsd) / nd;
      S = J2d*Sb - r1*Ci;
      r2  = r2*J4d;

      CCb = r2*(CCb + g1 * ten_dyad_prod(RmRm.slice(0), RmRm.slice(0), nsd) + 
                      g2 * (ten_dyad_prod(RmRm.slice(1), RmRm.slice(1), nsd) +
                           ten_dyad_prod(RmRm.slice(2), RmRm.slice(2), nsd) +
                           ten_dyad_prod(RmRm.slice(4), RmRm.slice(4), nsd)*2.0) +
                2.0 * g3 * (ten_dyad_prod(RmRm.slice(3), RmRm.slice(3), nsd) +
                ten_dyad_prod(RmRm.slice(5), RmRm.slice(5), nsd)));

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC = ten_ddot(CCb, PP, nsd);
      CC  = ten_transpose(CC, nsd);
      CC  = ten_ddot(PP, CC, nsd);
      CC  = CC - (2.0/nd) * ( ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd) );

      S  = S + p*J*Ci;
      CC = CC + 2.0*(r1 - p*J) * ten_symm_prod(Ci, Ci, nsd) + (pl*J - 2.0*r1/nd) * ten_dyad_prod(Ci, Ci, nsd);
    } break;

    //  HO (Holzapfel-Ogden) model for myocardium (2009)
    case ConstitutiveModelType::stIso_HO: {
      if (nfd != 2) {
        throw std::runtime_error("[get_pk2cc] Min fiber directions not defined for Holzapfel material model.");
      }
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
    } break;

    default:
      throw std::runtime_error("Undefined material constitutive model.");
  } 

  // Convert to Voigt Notation
  cc_to_voigt(nsd, CC, Dm);
}

//---------------
// get_pk2cc_dev
//---------------
// Compute isochoric (deviatoric) component of 2nd Piola-Kirchhoff stress and material stiffness tensors.
//
// Reproduces 'SUBROUTINE GETPK2CCdev(lDmn, F, nfd, fl, ya, S, Dm, Ja)'. 
//
void get_pk2cc_dev(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const Array<double>& F, const int nfd, 
    const Array<double>& fl, const double ya, Array<double>& S, Array<double>& Dm, double& Ja)
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_get_pk2cc_dev 
  #ifdef debug_get_pk2cc_dev 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int nsd = com_mod.nsd;
  S = 0.0;
  Dm = 0.0;

  // Some preliminaries
  auto& stM = lDmn.stM;
  double nd = static_cast<double>(nsd);

  // Fiber-reinforced stress
  double Tfa = 0.0;
  get_fib_stress(com_mod, cep_mod, stM.Tf, Tfa);

  // Electromechanics coupling - active stress
  if (cep_mod.cem.aStress) {
    Tfa = Tfa + ya;
  }

  // Electromechanics coupling - active strain
  auto Fe = F;
  auto Fa = mat_id(nsd);
  auto Fai = Fa;

  if (cep_mod.cem.aStrain) {
    actv_strain(com_mod, cep_mod, ya, nfd, fl, Fa);
    Fai = mat_inv(Fa, nsd);
    Fe = mat_mul(F, Fai);
  }

  #ifdef debug_get_pk2cc_dev 
  dmsg << "cep_mod.cem.aStress: " << cep_mod.cem.aStress;
  dmsg << "cep_mod.cem.aStrain: " << cep_mod.cem.aStrain;
  #endif

  Ja = mat_det(Fa, nsd);
  double J = mat_det(Fe, nsd);
  double J2d = pow(J, (-2.0/nd));
  double J4d = J2d * J2d;

  auto IDm = mat_id(nsd);
  auto C = mat_mul(transpose(Fe), Fe);
  auto E = 0.5 * (C - IDm);
  auto Ci = mat_inv(C, nsd);

  double trE = mat_trace(E, nsd);
  double Inv1 = J2d * mat_trace(C,nsd);
  double Inv2 = 0.5 * (Inv1*Inv1 - J4d*mat_trace(mat_mul(C,C), nsd));

  // Isochoric part of 2nd Piola-Kirchhoff and elasticity tensors
  //
  Tensor4<double> CC(nsd,nsd,nsd,nsd);

  switch (stM.isoType) {

    // NeoHookean model
    //
    case ConstitutiveModelType::stIso_nHook: {
      double g1 = 2.0 * stM.C10;
      auto Sb = g1 * IDm;

      // Fiber reinforcement/active stress
      auto u = fl.col(0);
      Sb += Tfa * mat_dyad_prod(u, u, nsd);

      double r1 = g1 * Inv1 / nd;
      S = J2d*Sb - r1*Ci;

      CC = 2.0 * r1 * (ten_symm_prod(Ci, Ci, nsd) - 1.0/nd * ten_dyad_prod(Ci, Ci, nsd)) - 
                2.0/nd * (ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd));
  
    } break; 

    // Mooney-Rivlin model
    //
    case ConstitutiveModelType::stIso_MR: {
      double g1 = 2.0 * (stM.C10 + Inv1*stM.C01);
      double g2 = -2.0 * stM.C01;
      auto Sb = g1*IDm + g2*J2d*C;

      // Fiber reinforcement/active stress
      auto u = fl.col(0);
      Sb += Tfa * mat_dyad_prod(u, u, nsd);
      g1 = 4.0 * J4d * stM.C01;
      auto CCb = g1 * (ten_dyad_prod(IDm, IDm, nsd) - ten_ids(nsd));

      double r1  = J2d*mat_ddot(C, Sb, nsd) / nd;
      S = J2d*Sb - r1*Ci;

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);

      CC = ten_ddot(CCb, PP, nsd);
      CC = ten_transpose(CC, nsd);
      CC = ten_ddot(PP, CC, nsd);

      CC = CC + 2.0 * r1 * (ten_symm_prod(Ci, Ci, nsd) - (1.0/nd) * ten_dyad_prod(Ci, Ci, nsd)) - 
                (2.0/nd) * (ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd));

    } break; 

    // HGO (Holzapfel-Gasser-Ogden) model with additive splitting of
    // the anisotropic fiber-based strain-energy terms
    //
    case ConstitutiveModelType::stIso_HGO: {
      if (nfd != 2) {
        throw std::runtime_error("[get_pk2cc_dev] Min fiber directions not defined for HGO material model.");
      }

      double kap = stM.kap;
      double Inv4 = J2d * norm(fl.col(0), mat_mul(C, fl.col(0)));
      double Inv6 = J2d * norm(fl.col(1), mat_mul(C, fl.col(1)));

      double Eff = kap*Inv1 + (1.0 - 3.0*kap) * Inv4 - 1.0;
      double Ess = kap*Inv1 + (1.0 - 3.0*kap) * Inv6 - 1.0;

      auto Hff = mat_dyad_prod(fl.col(0), fl.col(0), nsd);
      Hff = kap*IDm + (1.0 - 3.0*kap) * Hff;

      auto Hss = mat_dyad_prod(fl.col(1), fl.col(1), nsd);
      Hss  = kap*IDm + (1.0 - 3.0*kap) * Hss;

      double g1 = stM.C10;
      double g2 = stM.aff * Eff * exp(stM.bff*Eff*Eff);
      double g3 = stM.ass * Ess * exp(stM.bss*Ess*Ess);
      auto Sb = 2.0 * (g1*IDm + g2*Hff + g3*Hss);

      // Fiber reinforcement/active stress
      //
      Sb += Tfa *  mat_dyad_prod(fl.col(0), fl.col(0), nsd); 

      g1 = stM.aff*(1.0 + 2.0*stM.bff*Eff*Eff) * exp(stM.bff*Eff*Eff);
      g2 = stM.ass*(1.0 + 2.0*stM.bss*Ess*Ess) * exp(stM.bss*Ess*Ess);
      g1 = 4.0*J4d * g1;
      g2 = 4.0*J4d * g2;

      auto CCb = g1 * ten_dyad_prod(Hff, Hff, nsd) + g2 * ten_dyad_prod(Hss, Hss, nsd);

      double r1  = J2d*mat_ddot(C, Sb, nsd) / nd;
      S = J2d*Sb - r1*Ci;

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC = ten_ddot(CCb, PP, nsd);
      CC = ten_transpose(CC, nsd);
      CC = ten_ddot(PP, CC, nsd);

      CC = CC + 2.0*r1 * (ten_symm_prod(Ci, Ci, nsd) - 1.0/nd * ten_dyad_prod(Ci, Ci, nsd)) - 2.0/nd * (ten_dyad_prod(Ci, S, nsd) +
          ten_dyad_prod(S, Ci, nsd));
    } break; 

    // Guccione (1995) transversely isotropic model
    //
    case ConstitutiveModelType::stIso_Gucci: {
      if (nfd != 2) {
        throw std::runtime_error("[get_pk2cc_dev] Min fiber directions not defined for Guccione material model.");
      }

      // Compute isochoric component of E
      auto E = 0.5 * (J2d*C - IDm);

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

      auto QQ = g1 *  Es(0,0)*Es(0,0) + 
                g2 * (Es(1,1)*Es(1,1)  +  Es(2,2)*Es(2,2)  +  Es(1,2)*Es(1,2) + Es(2,1)*Es(2,1)) + 
                g3 * (Es(0,1)*Es(0,1)  +  Es(1,0)*Es(1,0)  +  Es(0,2)*Es(0,2) + Es(2,0)*Es(2,0));

      auto r2 = stM.C10 * exp(QQ);

      // Fiber stiffness contribution := (dE*_ab / dE_IJ)
      //
      Array3<double> RmRm(nsd,nsd,6);

      RmRm.set_slice(0, mat_dyad_prod(Rm.col(0), Rm.col(0), nsd));
      RmRm.set_slice(1, mat_dyad_prod(Rm.col(1), Rm.col(1), nsd));
      RmRm.set_slice(2, mat_dyad_prod(Rm.col(2), Rm.col(2), nsd));

      RmRm.set_slice(3, mat_symm_prod(Rm.col(0), Rm.col(1), nsd));
      RmRm.set_slice(4, mat_symm_prod(Rm.col(1), Rm.col(2), nsd));
      RmRm.set_slice(5, mat_symm_prod(Rm.col(2), Rm.col(0), nsd));
  
      auto Sb = g1*Es(0,0)*RmRm.slice(0) + 
               g2*(Es(1,1)*RmRm.slice(1) + Es(2,2)*RmRm.slice(2) + 2.0*Es(1,2)*RmRm.slice(4)) +
           2.0*g3*(Es(0,1)*RmRm.slice(3) + Es(0,2)*RmRm.slice(5));

      auto CCb = 2.0*ten_dyad_prod(Sb, Sb, nsd);
      Sb += Sb * r2;

      // Fiber reinforcement/active stress
      Sb += Sb + Tfa*mat_dyad_prod(fl.col(0), fl.col(0), nsd);

      double r1 = J2d*mat_ddot(C, Sb, nsd) / nd;
      S = J2d*Sb - r1*Ci;
      r2 = r2*J4d;

      CCb = r2*(CCb + g1*ten_dyad_prod(RmRm.slice(0), RmRm.slice(0), nsd) + 
                g2*(ten_dyad_prod(RmRm.slice(1), RmRm.slice(1), nsd) +
                ten_dyad_prod(RmRm.slice(2), RmRm.slice(2), nsd) +
                ten_dyad_prod(RmRm.slice(4), RmRm.slice(4), nsd)*2.0) +
                2.0*g3*(ten_dyad_prod(RmRm.slice(3), RmRm.slice(3), nsd) +
                ten_dyad_prod(RmRm.slice(5), RmRm.slice(5), nsd)));

      auto PP  = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC  = ten_ddot(CCb, PP, nsd);
      CC = ten_transpose(CC, nsd);
      CC = ten_ddot(PP, CC, nsd);
      CC = CC + 2.0*r1 * (ten_symm_prod(Ci, Ci, nsd) -
                          1.0/nd * ten_dyad_prod(Ci, Ci, nsd)) - 
                          2.0/nd * (ten_dyad_prod(Ci, S, nsd) +
                          ten_dyad_prod(S, Ci, nsd));
    } break; 

    // HO (Holzapfel-Ogden) model for myocardium (2009)
    //
    case ConstitutiveModelType::stIso_HO: {
      if (nfd != 2) {
        throw std::runtime_error("[get_pk2cc_dev] Min fiber directions not defined for Holzapfel material model.");
      }

      double Inv4 = J2d*norm(fl.col(0), mat_mul(C, fl.col(0)));
      double Inv6 = J2d*norm(fl.col(1), mat_mul(C, fl.col(1)));
      double Inv8 = J2d*norm(fl.col(0), mat_mul(C, fl.col(1)));

      double Eff = Inv4 - 1.0;
      double Ess = Inv6 - 1.0;
      double Efs = Inv8;

      double g1 = stM.a * exp(stM.b*(Inv1-3.0));
      double g2 = 2.0 * stM.afs * Efs * exp(stM.bfs*Efs*Efs);
      auto Hfs = mat_symm_prod(fl.col(0), fl.col(1), nsd);
      auto Sb = g1*IDm + g2*Hfs;

      Efs  = Efs * Efs;
      g1 = 2.0*J4d*stM.b*g1;
      g2 = 4.0*J4d*stM.afs*(1.0 + 2.0*stM.bfs*Efs)* exp(stM.bfs*Efs);
      auto CCb = g1 * ten_dyad_prod(IDm, IDm, nsd) + g2 * ten_dyad_prod(Hfs, Hfs, nsd);

      // Fiber reinforcement/active stress
      //
      if (Eff > 0.0) {
        g1 = Tfa;
        g1  = g1 + 2.0 * stM.aff * Eff * exp(stM.bff*Eff*Eff);
        auto Hff = mat_dyad_prod(fl.col(0), fl.col(0), nsd);
        Sb  = Sb + g1*Hff;

        Eff = Eff * Eff;
        g1  = 4.0*J4d*stM.aff*(1.0 + 2.0*stM.bff*Eff)*exp(stM.bff*Eff);
        CCb = CCb + g1*ten_dyad_prod(Hff, Hff, nsd);
      }

      if (Ess >  0.0) {
        g2 = 2.0 * stM.ass * Ess * exp(stM.bss*Ess*Ess);
        auto Hss = mat_dyad_prod(fl.col(1), fl.col(1), nsd);
        Sb  = Sb + g2*Hss;

        Ess = Ess * Ess;
        g2  = 4.0*J4d*stM.ass*(1.0 + 2.0*stM.bss*Ess)*exp(stM.bss*Ess);
        CCb = CCb + g2*ten_dyad_prod(Hss, Hss, nsd);
      }

      double r1 = J2d*mat_ddot(C, Sb, nsd) / nd;
      auto S = J2d*Sb - r1*Ci;

      auto PP = ten_ids(nsd) - (1.0/nd) * ten_dyad_prod(Ci, C, nsd);
      CC = ten_ddot(CCb, PP, nsd);
      CC = ten_transpose(CC, nsd);
      CC = ten_ddot(PP, CC, nsd);
      CC = CC + 2.0*r1 * (ten_symm_prod(Ci, Ci, nsd) - 1.0/nd * ten_dyad_prod(Ci, Ci, nsd)) -
                2.0/nd * (ten_dyad_prod(Ci, S, nsd) + ten_dyad_prod(S, Ci, nsd));

      if (cep_mod.cem.aStrain) {
        S = mat_fun::mat_mul(Fai, S);
        S = mat_fun::mat_mul(S, mat_fun::transpose(Fai));
        CCb = 0.0;
        CCb = ten_dyad_prod(Fai, Fai, nsd);
        CC = ten_ddot(CC, CCb, nsd);
        CC = ten_ddot(CCb, CC, nsd);
      }

    } break;

    default: 
      throw std::runtime_error("Undefined isochoric material constitutive model.");
  } 

  //  Convert to Voigt Notation
  cc_to_voigt(nsd, CC, Dm);
}

//----------------
// get_pk2cc_shlc
//----------------
// Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
// for compressible shell elements.
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

//----------------
// get_pk2cc_shli
//----------------
// Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
// for incompressible shell elements
//
// Reproduces Fortran GETPK2CC_SHLi
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


//------------
// get_svol_p
//------------
// Reproduces Fortran 'GETSVOLP'.
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

//---------
// get_tau
//---------
// Compute stabilization parameters tauM and tauC.
//
// Reproduces Fortran 'GETTAU'.
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

//-----------
// g_vol_pen
//-----------
//
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

};
