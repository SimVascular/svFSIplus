
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
