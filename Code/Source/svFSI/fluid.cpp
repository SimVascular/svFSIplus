
#include "fluid.h"

#include "all_fun.h"
#include "consts.h"
#include "fs.h"
#include "lhsa.h"
#include "nn.h"
#include "utils.h"

#include <array>
#include <iomanip>
#include <math.h>

#ifdef WITH_TRILINOS
#include "trilinos_linear_solver.h"
#endif

namespace fluid {

//---------
// b_fluid
//---------
//
void b_fluid(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Vector<double>& y, 
    const double h, const Vector<double>& nV, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_b_fluid
  #ifdef debug_b_fluid
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const int cDmn = com_mod.cDmn;
  const auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double wl  = w * eq.af * eq.gam * dt;
  double udn = 0.0;
  #ifdef debug_b_fluid
  dmsg << "nsd: " << nsd;
  dmsg << "w: " << w;
  dmsg << "h: " << h;
  dmsg << "nV: " << nV(0) << " " << nV(1) << " " << nV(2);
  dmsg << "wl: " << wl;
  dmsg << "com_mod.mvMsh: " << com_mod.mvMsh;
  #endif 

  Vector<double> u(nsd);

  if (com_mod.mvMsh) {
    for (int i = 0; i < nsd; i++) {
      int j = i + nsd + 1;
      u(i) = y(i) - y(j);
      udn  = udn + u(i)*nV(i);
    }

  } else {
    for (int i = 0; i < nsd; i++) {
      u(i) = y(i);
      udn  = udn + u(i)*nV(i);
    }
  }

  udn = 0.50 * dmn.prop.at(PhysicalProperyType::backflow_stab) * dmn.prop.at(PhysicalProperyType::fluid_density) * (udn - fabs(udn));
  auto hc  = h*nV + udn*u;
  #ifdef debug_b_fluid
  dmsg << "udn: " << udn;
  dmsg << "u: " << u(0) << " " << u(1) << " " << u(2);
  dmsg << "hc: " << hc(0) << " " << hc(1) << " " << hc(2);
  #endif

  // Here the loop is started for constructing left and right hand side
  //
  if (nsd == 2) {
    for (int a = 0; a < eNoN; a++) {
      lR(0,a) = lR(0,a) - w*N(a)*hc(0);
      lR(1,a) = lR(1,a) - w*N(a)*hc(1);

      for (int b = 0; b < eNoN; b++) {
        double T1 = wl*N(a)*N(b)*udn;
        lK(0,a,b) = lK(0,a,b) - T1;
        lK(4,a,b) = lK(4,a,b) - T1;
      }
    }

  } else {
    for (int a = 0; a < eNoN; a++) {
      lR(0,a) = lR(0,a) - w*N(a)*hc(0);
      lR(1,a) = lR(1,a) - w*N(a)*hc(1);
      lR(2,a) = lR(2,a) - w*N(a)*hc(2);

      for (int b = 0; b < eNoN; b++) {
        double T1 = wl*N(a)*N(b)*udn;
        lK(0,a,b)  = lK(0,a,b)  - T1;
        lK(5,a,b)  = lK(5,a,b)  - T1;
        lK(10,a,b) = lK(10,a,b) - T1;
      }
    }
  }
} 

//-------------
// bw_fluid_2d
//-------------
//
void bw_fluid_2d(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, const Vector<double>& Nw, 
    const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& yl, const Vector<double>& ub, 
    const Vector<double>& nV, const Vector<double>& tauB, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_bw_fluid_2d
  #ifdef debug_bw_fluid_2d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double rho = dmn.prop.at(PhysicalProperyType::fluid_density);
  double tauT = tauB(0);
  double tauN = tauB(1);

  double T1 = eq.af * eq.gam * dt;
  double wl = w * T1;

  #ifdef debug_bw_fluid_2d 
  dmsg << "tauT: " << tauT;
  dmsg << "tauN: " << tauN;
  dmsg << "T1: " << T1;
  dmsg << "yl: " << yl;
  dmsg << "ub: " << ub;
  #endif

  Vector<double> u(2), Nxn(eNoNw); 
  Array<double> ux(2,2);

  for (int a = 0; a < eNoNw; a++) {
    u(0) = u(0) + Nw(a)*yl(0,a);
    u(1) = u(1) + Nw(a)*yl(1,a);

    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);

    Nxn(a) = Nwx(0,a)*nV(0) + Nwx(1,a)*nV(1);
  }

  double p = 0.0;
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(2,a);
  }

  Vector<double> uh(2); 

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      uh(0) = uh(0) + Nw(a)*yl(3,a);
      uh(1) = uh(1) + Nw(a)*yl(4,a);
    }
  }

  double un = (u(0)-uh(0))*nV(0) + (u(1)-uh(1))*nV(1);
  un = (fabs(un) - un) * 0.50;

  u = u - ub;
  double ubn  = u(0)*nV(0) + u(1)*nV(1);

  //  Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  Array<double> es(2,2);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,1) = ux(1,1) + ux(1,1);
  es(1,0) = ux(1,0) + ux(0,1);
  es(0,1) = es(1,0);

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(0,1)*es(0,1) + es(1,1)*es(1,1);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  // sigma.n (deviatoric)
  Vector<double> sgmn(2);
  sgmn(0) = mu*(es(0,0)*nV(0) + es(1,0)*nV(1));
  sgmn(1) = mu*(es(0,1)*nV(0) + es(1,1)*nV(1));

  Vector<double> rV(2);
  rV(0) = p*nV(0) - sgmn(0) + (tauT + rho*un)*u(0) + (tauN-tauT)*ubn*nV(0);
  rV(1) = p*nV(1) - sgmn(1) + (tauT + rho*un)*u(1) + (tauN-tauT)*ubn*nV(1);

  Array<double> rM(2,2);
  rM(0,0) = -mu*( u(0)*nV(0) + u(0)*nV(0) );
  rM(1,0) = -mu*( u(0)*nV(1) + u(1)*nV(0) );
  rM(0,1) = -mu*( u(1)*nV(0) + u(0)*nV(1) );
  rM(1,1) = -mu*( u(1)*nV(1) + u(1)*nV(1) );

  // Local residue (momentum)
  //
  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + w*(Nw(a)*rV(0) + Nwx(0,a)*rM(0,0) + Nwx(1,a)*rM(1,0));
    lR(1,a) = lR(1,a) + w*(Nw(a)*rV(1) + Nwx(0,a)*rM(0,1) + Nwx(1,a)*rM(1,1));
  }

  // Local residue (continuity)
  for (int a = 0; a < eNoNq; a++) {
    lR(2,a) = lR(2,a) - w*Nq(a)*ubn;
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T3 = Nw(a)*Nw(b);
      double T1 = (tauT + rho*un)*T3 - mu*(Nw(a)*Nxn(b) + Nxn(a)*Nw(b));
      double T2{0};

     // dRm_a1/du_b1
     T2 = (tauN - tauT)*T3*nV(0)*nV(0) - mu*(Nw(a)*Nwx(0,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(0));
     lK(0,a,b) = lK(0,a,b) + wl*(T2 + T1);

     // dRm_a1/du_b2
     T2 = (tauN - tauT)*T3*nV(0)*nV(1) - mu*(Nw(a)*Nwx(0,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(0));
     lK(1,a,b) = lK(1,a,b) + wl*T2;

     // dRm_a2/du_b1
     T2 = (tauN - tauT)*T3*nV(1)*nV(0) - mu*(Nw(a)*Nwx(1,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(1));
     lK(3,a,b) = lK(3,a,b) + wl*T2;

     // dRm_a2/du_b2
     T2 = (tauN - tauT)*T3*nV(1)*nV(1) - mu*(Nw(a)*Nwx(1,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(1));
     lK(4,a,b) = lK(4,a,b)  + wl*(T2 + T1);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T1 = wl*Nw(a)*Nq(b);

      // dRm_a1/dp_b
      lK(2,a,b) = lK(2,a,b) + T1*nV(0);

      // dRm_a2/dp_b
      lK(5,a,b) = lK(5,a,b) + T1*nV(1);
    }
  }

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      double T1 = wl*Nq(a)*Nw(b);

      // dRc_a/du_b1
      lK(6,a,b) = lK(6,a,b) - T1*nV(0);

      // dRc_a/du_b2
      lK(7,a,b) = lK(7,a,b) - T1*nV(1);
    }
  }
}

//-------------
// bw_fluid_3d
//-------------
//
void bw_fluid_3d(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, const Vector<double>& Nw, 
    const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& yl, const Vector<double>& ub, 
    const Vector<double>& nV, const Vector<double>& tauB, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double rho = dmn.prop.at(PhysicalProperyType::fluid_density);
  double tauT = tauB(0);
  double tauN = tauB(1);

  double T1 = eq.af * eq.gam * dt;
  double wl = w * T1;

  Vector<double> u(3), Nxn(eNoNw); 
  Array<double> ux(3,3);

  for (int a = 0; a < eNoNw; a++) {
    u(0) = u(0) + Nw(a)*yl(0,a);
    u(1) = u(1) + Nw(a)*yl(1,a);
    u(2) = u(2) + Nw(a)*yl(2,a);

    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    ux(1,0) = ux(1,1) + Nwx(1,a)*yl(0,a);
    ux(2,0) = ux(2,0) + Nwx(2,a)*yl(0,a);
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);
    ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a);
    ux(0,2) = ux(0,2) + Nwx(0,a)*yl(2,a);
    ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a);
    ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a);

    Nxn(a)  = Nwx(0,a)*nV(0) + Nwx(1,a)*nV(1) + Nwx(2,a)*nV(2);
  }

  double p = 0.0;
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(3,a);
  }

  Vector<double> uh(3); 

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      uh(0) = uh(0) + Nw(a)*yl(4,a);
      uh(1) = uh(1) + Nw(a)*yl(5,a);
      uh(2) = uh(2) + Nw(a)*yl(6,a);
    }
  }

  double un = (u(0)-uh(0))*nV(0) + (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2);
  un = (fabs(un) - un) * 0.50;

  u = u - ub;
  double ubn  = u(0)*nV(0) + u(1)*nV(1) + u(2)*nV(2);

  //  Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  Array<double> es(3,3);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,1) = ux(1,1) + ux(1,1);
  es(2,2) = ux(2,2) + ux(2,2);
  es(1,0) = ux(1,0) + ux(0,1);
  es(2,1) = ux(2,1) + ux(1,2);
  es(0,2) = ux(0,2) + ux(2,0);
  es(0,1) = es(1,0);
  es(1,2) = es(2,1);
  es(2,0) = es(0,2);

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + 
               es(0,1)*es(0,1) + es(1,1)*es(1,1);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  // sigma.n (deviatoric)
  Vector<double> sgmn(2);
  sgmn(0) = mu*(es(0,0)*nV(0) + es(1,0)*nV(1));
  sgmn(1) = mu*(es(0,1)*nV(0) + es(1,1)*nV(1));

  Vector<double> rV(2);
  rV(0) = p*nV(0) - sgmn(0) + (tauT + rho*un)*u(0) + (tauN-tauT)*ubn*nV(0);
  rV(1) = p*nV(1) - sgmn(1) + (tauT + rho*un)*u(1) + (tauN-tauT)*ubn*nV(1);

  Array<double> rM(2,2);
  rM(0,0) = -mu*( u(0)*nV(0) + u(0)*nV(0) );
  rM(1,0) = -mu*( u(0)*nV(1) + u(1)*nV(0) );
  rM(0,1) = -mu*( u(1)*nV(0) + u(0)*nV(1) );
  rM(1,1) = -mu*( u(1)*nV(1) + u(1)*nV(1) );

  // Local residue (momentum)
  //
  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + w*(Nw(a)*rV(0) + Nwx(0,a)*rM(0,0) + Nwx(1,a)*rM(1,0));
    lR(1,a) = lR(1,a) + w*(Nw(a)*rV(1) + Nwx(0,a)*rM(0,1) + Nwx(1,a)*rM(1,1));
  }

  // Local residue (continuity)
  for (int a = 0; a < eNoNq; a++) {
    lR(2,a) = lR(2,a) - w*Nq(a)*ubn;
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T3 = Nw(a)*Nw(b);
      double T1 = (tauT + rho*un)*T3 - mu*(Nw(a)*Nxn(b) + Nxn(a)*Nw(b));
      double T2{0.0};

      // dRm_a1/du_b1
      T2 = (tauN - tauT)*T3*nV(0)*nV(0) - mu*(Nw(a)*Nwx(0,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(0));
      lK(0,a,b) = lK(0,a,b) + wl*(T2 + T1);

      // dRm_a1/du_b2
      T2 = (tauN - tauT)*T3*nV(0)*nV(1) - mu*(Nw(a)*Nwx(0,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(0));
      lK(1,a,b) = lK(1,a,b) + wl*T2;

      // dRm_a2/du_b1
      T2 = (tauN - tauT)*T3*nV(1)*nV(0) - mu*(Nw(a)*Nwx(1,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(1));
      lK(3,a,b) = lK(3,a,b) + wl*T2;

      // dRm_a2/du_b2
      T2 = (tauN - tauT)*T3*nV(1)*nV(1) - mu*(Nw(a)*Nwx(1,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(1));
      lK(4,a,b) = lK(4,a,b)  + wl*(T2 + T1);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T1 = wl*Nw(a)*Nq(b);

      // dRm_a1/dp_b
      lK(2,a,b)  = lK(2,a,b)  + T1*nV(0);

      // dRm_a2/dp_b
      lK(5,a,b) = lK(5,a,b)  + T1*nV(1);
    }
  }

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      double T1 = wl*Nq(a)*Nw(b);

      // dRc_a/du_b1
      lK(6,a,b) = lK(6,a,b) - T1*nV(0);

      // dRc_a/du_b2
      lK(7,a,b) = lK(7,a,b) - T1*nV(1);
    }
  }
}

//-----------------
// construct_fluid
//-----------------
// This is for solving fluid transport equation solving Navier-Stokes
// equations. Dirichlet boundary conditions are either treated
// strongly or weakly.
//
// Modifies:
//
void construct_fluid(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg)
{
  #define n_debug_construct_fluid
  #ifdef debug_construct_fluid
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  double start_time = utils::cput();
  #endif

  using namespace consts;

  const int eNoN = lM.eNoN;
  bool vmsStab = false;

  if (lM.nFs == 1) {
     vmsStab = true;
  } else {
     vmsStab = false;
  }

  // l = 3, if nsd==2 ; else 6;
  const int l = com_mod.nsymd;
  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  #ifdef debug_construct_fluid
  dmsg << "cEq: " << cEq;
  dmsg << "eq.sym: " << eq.sym;
  dmsg << "eq.dof: " << eq.dof;
  dmsg << "eNoN: " << eNoN;
  dmsg << "vmsStab: " << vmsStab;
  dmsg << "l: " << l;
  dmsg << "tDof: " << tDof;
  dmsg << "dof: " << dof;
  dmsg << "lM.nEl: " <<  lM.nEl;
  dmsg << "nsd: " <<  nsd;
  #endif

  // FLUID: dof = nsd+1
  Vector<int> ptr(eNoN); 
  Array<double> xl(nsd,eNoN); 
  Array<double> al(tDof,eNoN); 
  Array<double> yl(tDof,eNoN); 
  Array<double> bfl(nsd,eNoN); 
  Array<double> lR(dof,eNoN); 
  Array3<double> lK(dof*dof,eNoN,eNoN);

  // Loop over all elements of mesh
  //
  int num_c = lM.nEl / 10;

  for (int e = 0; e < lM.nEl; e++) {
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;

    if (cPhys != EquationType::phys_fluid) {
      continue;
    }

    //  Update shape functions for NURBS
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    // Create local copies
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
      }
    }

    // Initialize residue and tangents
    lR = 0.0;
    lK = 0.0;
    std::array<fsType,2> fs;

    // Set function spaces for velocity and pressure.
    fs::get_thood_fs(com_mod, fs, lM, vmsStab, 1);

    // Define element coordinates appropriate for function spaces
    Array<double> xwl(nsd,fs[0].eNoN); 
    Array<double> Nwx(nsd,fs[0].eNoN); 
    Array<double> Nwxx(l,fs[0].eNoN);

    Array<double> xql(nsd,fs[1].eNoN); 
    Array<double> Nqx(nsd,fs[1].eNoN);

    #ifdef debug_construct_fluid
    dmsg;
    dmsg << "l: " << l;
    dmsg << "fs[0].eNoN: " << fs[0].eNoN;
    dmsg << "fs[1].eNoN: " << fs[1].eNoN;
    #endif

    xwl = xl;

    for (int i = 0; i < xql.nrows(); i++) { 
      for (int j = 0; j < fs[1].eNoN; j++) { 
        xql(i,j) = xl(i,j);
      }
    }

    // Gauss integration 1
    //
    #ifdef debug_construct_fluid
    dmsg;
    dmsg << "Gauss integration 1 ... ";
    dmsg << "fs[1].nG: " << fs[0].nG;
    dmsg << "fs[1].lShpF: " << fs[0].lShpF;
    dmsg << "fs[2].nG: " << fs[1].nG;
    dmsg << "fs[2].lShpF: " << fs[1].lShpF;
    #endif

    double Jac{0.0};
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < fs[0].nG; g++) {
      if (g == 0 || !fs[1].lShpF) {
        auto Nx = fs[1].Nx.rslice(g);
        nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs[0].lShpF) {
        auto Nx = fs[0].Nx.rslice(g);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }

        if (!vmsStab) {
          auto Nx = fs[0].Nx.rslice(g);
          auto Nxx = fs[0].Nxx.rslice(g);
          nn::gn_nxx(l, fs[0].eNoN, nsd, nsd, Nx, Nxx, xwl, Nwx, Nwxx); 
        }
      }

      double w = fs[0].w(g) * Jac;

      // Compute momentum residual and tangent matrix.
      //
      if (nsd == 3) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_3d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, 
            Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_2d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, 
            Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);
      }
    } // g: loop

    // Set function spaces for velocity and pressure.
    //
    fs::get_thood_fs(com_mod, fs, lM, vmsStab, 2);

    // Gauss integration 2
    //
    for (int g = 0; g < fs[1].nG; g++) {
      if (g == 0 || !fs[0].lShpF) {
        auto Nx = fs[0].Nx.rslice(g);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs[1].lShpF) {
        auto Nx = fs[1].Nx.rslice(g);
        nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }
      double w = fs[1].w(g) * Jac;

      // Compute continuity residual and tangent matrix.
      //
      if (nsd == 3) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_3d_c(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_2d_c(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);
      }

    } // g: loop

    // Assembly
    //dmsg << "Assembly  ... ";

#ifdef WITH_TRILINOS
    if (eq.assmTLS) {
      trilinos_doassem_(const_cast<int&>(eNoN), const_cast<int*>(ptr.data()), lK.data(), lR.data());
    } else {
#endif
      lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
#ifdef WITH_TRILINOS
    }

#endif 

  } // e: loop

  #ifdef debug_construct_fluid
  double end_time = utils::cput();
  double etime = end_time - start_time;
  #endif
}

//------------
// fluid_2d_c
//------------
// Reproduces Fortran 'FLUID2D_C()'.
//
void fluid_2d_c(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_fluid_2d_c 
  #ifdef debug_fluid_2d_c 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  dmsg << "Nw: " << Nw;
  dmsg << "Nq: " << Nq;
  dmsg << "Nwx: " << Nwx;
  dmsg << "Nqx: " << Nqx;
  dmsg << "w: " << w;
  dmsg << "bfl: " << bfl;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  const double ctM = 1.0;
  const double ctC = 36.0;

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  Vector<double> f(2);
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  f[1] = dmn.prop[PhysicalProperyType::f_y];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  #ifdef debug_fluid_2d_c 
  dmsg << "T1: " << T1;
  dmsg << "amd: " << amd;
  dmsg << "wl: " << wl;
  #endif

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  Vector<double> ud{-f[0], -f[1]};
  Vector<double> u(2);
  Array<double> ux(2,2);
  Array3<double> uxx(2,2,2);

  for (int a = 0; a < eNoNw; a++) {
    ud(0) = ud(0) + Nw(a)*(al(0,a)-bfl(0,a));
    ud(1) = ud(1) + Nw(a)*(al(1,a)-bfl(1,a));

    u(0) = u(0) + Nw(a)*yl(0,a);
    u(1) = u(1) + Nw(a)*yl(1,a);

    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);

    uxx(0,0,0) = uxx(0,0,0) + Nwxx(0,a)*yl(0,a);
    uxx(1,0,1) = uxx(1,0,1) + Nwxx(1,a)*yl(0,a);
    uxx(1,0,0) = uxx(1,0,0) + Nwxx(2,a)*yl(0,a);

    uxx(0,1,0) = uxx(0,1,0) + Nwxx(0,a)*yl(1,a);
    uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a);
    uxx(1,1,0) = uxx(1,1,0) + Nwxx(2,a)*yl(1,a);
  }

  double divU = ux(0,0) + ux(1,1);
  uxx(0,0,1) = uxx(1,0,0);
  uxx(0,1,1) = uxx(1,1,0);

  Vector<double> d2u2(2);
  d2u2(0) = uxx(0,0,0) + uxx(1,0,1);
  d2u2(1) = uxx(0,1,0) + uxx(1,1,1);


  // Pressure and its gradient
  Vector<double> px(2);
  for (int a = 0; a < eNoNq; a++) {
    px(0) = px(0) + Nqx(0,a)*yl(2,a);
    px(1) = px(1) + Nqx(1,a)*yl(2,a);
  }

  //  Update convection velocity relative to mesh velocity
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u(0) = u(0) - Nw(a)*yl(3,a);
      u(1) = u(1) - Nw(a)*yl(4,a);
     }
  }

  // Strain rate tensor 2*e_ij := (u_ij + u_ji)
  Array<double> es(2,2);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,0) = ux(1,0) + ux(0,1);
  es(1,1) = ux(1,1) + ux(1,1);
  es(0,1) = es(1,0);

  #ifdef debug_fluid_2d_c 
  dmsg << "ud: " << ud;
  dmsg << "u: " << u;
  dmsg << "ux: " << ux;
  dmsg << "uxx: " << uxx;
  dmsg << "px: " << px;
  dmsg << "es: " << es;
  #endif

  Array<double> esNx(2,eNoNw);
  for (int a = 0; a < eNoNw; a++) {
    esNx(0,a) = es(0,0)*Nwx(0,a) + es(1,0)*Nwx(1,a);
    esNx(1,a) = es(0,1)*Nwx(0,a) + es(1,1)*Nwx(1,a);
  }

  Array3<double> es_x(2,2,2);
  for (int k = 0; k < 2; k++) { 
    es_x(0,0,k) = uxx(0,0,k) + uxx(0,0,k);
    es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k);
    es_x(1,0,k) = uxx(1,0,k) + uxx(0,1,k);
    es_x(0,1,k) = es_x(1,0,k);
  }

  Vector<double> mu_x(2);
  mu_x(0) = (es_x(0,0,0)*es(0,0) + es_x(1,1,0)*es(1,1))*0.50 +  es_x(1,0,0)*es(1,0);
  mu_x(1) = (es_x(0,0,1)*es(0,0) + es_x(1,1,1)*es(1,1))*0.50 +  es_x(1,0,1)*es(1,0);

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(0,1)*es(0,1) + es(1,1)*es(1,1);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     mu_g = mu_g / gam;
  }
  std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
  //mu_x(:) = mu_g * mu_x(:)

  #ifdef debug_fluid_2d_c 
  dmsg;
  dmsg << "esNx: " << esNx;
  dmsg << "es_x: " << es_x;
  dmsg << "mu_x: " << mu_x;
  dmsg << "gam: " << gam;
  #endif

  // Stabilization parameters
  //
  Array3<double> updu(2,2,eNoNw);
  Vector<double> up(2);
  double tauM{0.0};

  if (vmsFlag) {
    double kT = 4.0 * pow(ctM/dt,2.0);
    double kU = u(0)*u(0)*Kxi(0,0) + u(1)*u(0)*Kxi(1,0) + u(0)*u(1)*Kxi(0,1) + u(1)*u(1)*Kxi(1,1);
    double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1);

    kS = ctC * kS * pow(mu/rho,2.0);
    tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

    Vector<double> rV(2);
    rV(0) = ud(0) + u(0)*ux(0,0) + u(1)*ux(1,0);
    rV(1) = ud(1) + u(0)*ux(0,1) + u(1)*ux(1,1);

    Vector<double> rS(2);
    rS(0) = mu_x(0)*es(0,0) + mu_x(1)*es(1,0) + mu*d2u2(0);
    rS(1) = mu_x(0)*es(0,1) + mu_x(1)*es(1,1) + mu*d2u2(1);

    up(0) = -tauM*(rho*rV(0) + px(0) - rS(0));
    up(1) = -tauM*(rho*rV(1) + px(1) - rS(1));

    for (int a = 0; a < eNoNw; a++) {
      double uNx = u(0)*Nwx(0,a) + u(1)*Nwx(1,a);
      T1 = -rho*uNx + mu*(Nwxx(0,a) + Nwxx(1,a)) + mu_x(0)*Nwx(0,a) + mu_x(1)*Nwx(1,a);

      updu(0,0,a) = mu_x(0)*Nwx(0,a) + d2u2(0)*mu_g*esNx(0,a) + T1;
      updu(1,0,a) = mu_x(1)*Nwx(0,a) + d2u2(0)*mu_g*esNx(1,a);

      updu(0,1,a) = mu_x(0)*Nwx(1,a) + d2u2(1)*mu_g*esNx(0,a);
      updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1;
    } 

    #ifdef debug_fluid_2d_c 
    dmsg;
    dmsg << "kT: " << kT;
    dmsg << "kS: " << kS;
    dmsg << "tauM: " << tauM;
    dmsg << "rV: " << rV;
    dmsg << "rS: " << rS;
    dmsg << "up: " << up;
    dmsg << "updu: " << updu;
    #endif

  } else {
    tauM = 0.0;
    up = 0.0;
    updu = 0.0;
  }

  // Local residue
  //
  for (int a = 0; a < eNoNq; a++) {
    double upNx = up(0)*Nqx(0,a) + up(1)*Nqx(1,a);
    lR(2,a) = lR(2,a) + w*(Nq(a)*divU - upNx);
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    T1 = rho*amd*Nw(b);

    for (int a = 0; a < eNoNq; a++) {
      // dRc_a/dU_b1
      double T2 = Nqx(0,a)*(updu(0,0,b) - T1) + Nqx(1,a)*updu(0,1,b);
      lK(6,a,b) = lK(6,a,b) + wl*(Nq(a)*Nwx(0,b) - tauM*T2);

      // dRc_a/dU_b2
      T2 = Nqx(0,a)*updu(1,0,b) + Nqx(1,a)*(updu(1,1,b) - T1);
      lK(7,a,b) = lK(7,a,b) + wl*(Nq(a)*Nwx(1,b) - tauM*T2);
    }
  }

  if (vmsFlag) {
    for (int b = 0; b < eNoNq; b++) {
      for (int a = 0; a < eNoNq; a++) {
        // dC/dP
        double NxNx = Nqx(0,a)*Nqx(0,b) + Nqx(1,a)*Nqx(1,b);
        lK(8,a,b) = lK(8,a,b) + wl*tauM*NxNx;
      }
    }
  }
}

//------------
// fluid_2d_m
//------------
//
void fluid_2d_m(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_fluid_2d_m
  #ifdef debug_fluid_2d_m
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double ctM = 1.0;
  double ctC = 36.0;

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  Vector<double> f(2);
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  f[1] = dmn.prop[PhysicalProperyType::f_y];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  double wr = w*rho;

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  Vector<double> ud{-f[0], -f[1]};
  Vector<double> u(2);
  Array<double> ux(2,2);
  Array3<double> uxx(2,2,2);

  for (int a = 0; a < eNoNw; a++) {
    ud(0) = ud(0) + Nw(a)*(al(0,a)-bfl(0,a));
    ud(1) = ud(1) + Nw(a)*(al(1,a)-bfl(1,a));

    u(0) = u(0) + Nw(a)*yl(0,a);
    u(1) = u(1) + Nw(a)*yl(1,a);

    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);

    uxx(0,0,0) = uxx(0,0,0) + Nwxx(0,a)*yl(0,a);
    uxx(1,0,1) = uxx(1,0,1) + Nwxx(1,a)*yl(0,a);
    uxx(1,0,0) = uxx(1,0,0) + Nwxx(2,a)*yl(0,a);

    uxx(0,1,0) = uxx(0,1,0) + Nwxx(0,a)*yl(1,a);
    uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a);
    uxx(1,1,0) = uxx(1,1,0) + Nwxx(2,a)*yl(1,a);
  }

  double divU = ux(0,0) + ux(1,1);
  uxx(0,0,1) = uxx(1,0,0);
  uxx(0,1,1) = uxx(1,1,0);

  Vector<double> d2u2(2);
  d2u2(0) = uxx(0,0,0) + uxx(1,0,1);
  d2u2(1) = uxx(0,1,0) + uxx(1,1,1);

  // Pressure and its gradient
  double p = 0.0;
  Vector<double> px(2);
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(2,a);
    px(0) = px(0) + Nqx(0,a)*yl(2,a);
    px(1) = px(1) + Nqx(1,a)*yl(2,a);
  }

  //  Update convection velocity relative to mesh velocity
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u(0) = u(0) - Nw(a)*yl(3,a);
      u(1) = u(1) - Nw(a)*yl(4,a);
     }
  }

  // Strain rate tensor 2*e_ij := (u_ij + u_ji)
  Array<double> es(2,2);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,1) = ux(1,1) + ux(1,1);
  es(1,0) = ux(1,0) + ux(0,1);
  es(0,1) = es(1,0);

  Array<double> esNx(2,eNoNw);
  for (int a = 0; a < eNoNw; a++) {
    esNx(0,a) = es(0,0)*Nwx(0,a) + es(1,0)*Nwx(1,a);
    esNx(1,a) = es(0,1)*Nwx(0,a) + es(1,1)*Nwx(1,a);
  }

  Array3<double> es_x(2,2,2);
  for (int k = 0; k < 2; k++) { 
    es_x(0,0,k) = uxx(0,0,k) + uxx(0,0,k);
    es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k);
    es_x(1,0,k) = uxx(1,0,k) + uxx(0,1,k);
    es_x(0,1,k) = es_x(1,0,k);
  }

  Vector<double> mu_x(2);
  mu_x(0) = (es_x(0,0,0)*es(0,0) + es_x(1,1,0)*es(1,1))*0.50 +  es_x(1,0,0)*es(1,0);
  mu_x(1) = (es_x(0,0,1)*es(0,0) + es_x(1,1,1)*es(1,1))*0.50 +  es_x(1,0,1)*es(1,0);

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(0,1)*es(0,1) + es(1,1)*es(1,1);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     mu_g = mu_g / gam;
  }
  std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
  //mu_x(:) = mu_g * mu_x(:)

  //  Stabilization parameters
  double kT = 4.0 * pow(ctM/dt,2.0);

  double kU = u(0)*u(0)*Kxi(0,0) + u(1)*u(0)*Kxi(1,0) + u(0)*u(1)*Kxi(0,1) + u(1)*u(1)*Kxi(1,1);
  double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1);
  kS = ctC * kS * pow(mu/rho,2.0);
  double tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

  Vector<double> rV(2);
  rV(0) = ud(0) + u(0)*ux(0,0) + u(1)*ux(1,0);
  rV(1) = ud(1) + u(0)*ux(0,1) + u(1)*ux(1,1);

  Vector<double> rS(2);
  rS(0) = mu_x(0)*es(0,0) + mu_x(1)*es(1,0) + mu*d2u2(0);
  rS(1) = mu_x(0)*es(0,1) + mu_x(1)*es(1,1) + mu*d2u2(1);

  Vector<double> up(2);
  up(0) = -tauM*(rho*rV(0) + px(0) - rS(0));
  up(1) = -tauM*(rho*rV(1) + px(1) - rS(1));

  double tauC, tauB, pa;
  double eps = std::numeric_limits<double>::epsilon();
  Vector<double> ua(2);

  if (vmsFlag) {
     tauC = 1.0 / (tauM * (Kxi(0,0) + Kxi(1,1)));
     tauB = up(0)*up(0)*Kxi(0,0) + up(1)*up(0)*Kxi(1,0) + up(0)*up(1)*Kxi(0,1) + up(1)*up(1)*Kxi(1,1);

    if (utils::is_zero(tauB)) {
      tauB = eps;
    }
     tauB = rho / sqrt(tauB);

     ua(0) = u(0) + up(0);
     ua(1) = u(1) + up(1);
     pa = p - tauC*divU;
  } else {
     tauC = 0.0;
     tauB = 0.0;
     ua = u;
     pa = p;
  }

  rV(0) = tauB*(up(0)*ux(0,0) + up(1)*ux(1,0));
  rV(1) = tauB*(up(0)*ux(0,1) + up(1)*ux(1,1));

  Array<double> rM(2,2);
  rM(0,0) = mu*es(0,0) - rho*up(0)*ua(0) + rV(0)*up(0) - pa;
  rM(1,0) = mu*es(1,0) - rho*up(0)*ua(1) + rV(0)*up(1);

  rM(0,1) = mu*es(0,1) - rho*up(1)*ua(0) + rV(1)*up(0);
  rM(1,1) = mu*es(1,1) - rho*up(1)*ua(1) + rV(1)*up(1) - pa;

  rV(0) = ud(0) + ua(0)*ux(0,0) + ua(1)*ux(1,0);
  rV(1) = ud(1) + ua(0)*ux(0,1) + ua(1)*ux(1,1);

  // Local residue
  Array3<double> updu(2,2,eNoNw);
  Vector<double> uNx(eNoNw), upNx(eNoNw), uaNx(eNoNw); 

  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + wr*Nw(a)*rV(0) + w*(Nwx(0,a)*rM(0,0) + Nwx(1,a)*rM(1,0));
    lR(1,a) = lR(1,a) + wr*Nw(a)*rV(1) + w*(Nwx(0,a)*rM(0,1) + Nwx(1,a)*rM(1,1));

    // Quantities used for stiffness matrix
    uNx(a) = u(0)*Nwx(0,a) + u(1)*Nwx(1,a);
    upNx(a) = up(0)*Nwx(0,a) + up(1)*Nwx(1,a);

    if (vmsFlag) {
      uaNx(a) = uNx(a) + upNx(a);
    } else {
      uaNx(a) = uNx(a);
    }

    T1 = -rho*uNx(a) + mu*(Nwxx(0,a) + Nwxx(1,a)) + mu_x(0)*Nwx(0,a) + mu_x(1)*Nwx(1,a);

    updu(0,0,a) = mu_x(0)*Nwx(0,a) + d2u2(0)*mu_g*esNx(0,a) + T1;
    updu(1,0,a) = mu_x(1)*Nwx(0,a) + d2u2(0)*mu_g*esNx(1,a);

    updu(0,1,a) = mu_x(0)*Nwx(1,a) + d2u2(1)*mu_g*esNx(0,a);
    updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1;
  }

  // Tangent (stiffness) matrices
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      rM(0,0) = Nwx(0,a)*Nwx(0,b);
      rM(1,0) = Nwx(1,a)*Nwx(0,b);
      rM(0,1) = Nwx(0,a)*Nwx(1,b);
      rM(1,1) = Nwx(1,a)*Nwx(1,b);

      double NxNx = Nwx(0,a)*Nwx(0,b) + Nwx(1,a)*Nwx(1,b);
      T1 = mu*NxNx + rho*amd*Nw(b)*(Nw(a) + rho*tauM*uaNx(a)) + rho*Nw(a)*(uNx(b)+upNx(b)) + tauB*upNx(a)*upNx(b);

      // dRm_a1/du_b1
      double T2 = (mu + tauC)*rM(0,0) + esNx(0,a)*mu_g*esNx(0,b) - rho*tauM*uaNx(a)*updu(0,0,b);
      lK(0,a,b) = lK(0,a,b)  + wl*(T2 + T1);

      // dRm_a1/du_b2
      T2 = mu*rM(1,0) + tauC*rM(0,1) + esNx(0,a)*mu_g*esNx(1,b) - rho*tauM*uaNx(a)*updu(1,0,b);
      lK(1,a,b) = lK(1,a,b) + wl*(T2);

      // dRm_a2/du_b1
      T2 = mu*rM(0,1) + tauC*rM(1,0) + esNx(1,a)*mu_g*esNx(0,b) - rho*tauM*uaNx(a)*updu(0,1,b);
      lK(3,a,b) = lK(3,a,b) + wl*(T2);

      // dRm_a2/du_b2
      T2 = (mu + tauC)*rM(1,1) + esNx(1,a)*mu_g*esNx(1,b) - rho*tauM*uaNx(a)*updu(1,1,b);
      lK(4,a,b) = lK(4,a,b) + wl*(T2 + T1);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      T1 = rho*tauM*uaNx(a);

      // dRm_a1/dp_b
      lK(2,a,b)  = lK(2,a,b)  - wl*(Nwx(0,a)*Nq(b) - Nqx(0,b)*T1);

      // dRm_a2/dp_b
      lK(5,a,b) = lK(5,a,b)  - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T1);
    }
  }
}

//-----------
// fluid_3d_c
//-----------
// Element continuity residue.
//
void fluid_3d_c(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_fluid3d_c
  #ifdef debug_fluid3d_c
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  double start_time = utils::cput();
  #endif
  
  // Maximum size of arrays sized by (3,eNoNw) -> (3,MAX_SIZE).
  const int MAX_SIZE = 8;

  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  const double ctM  = 1.0;
  const double ctC  = 36.0;

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  double f[3];
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  f[1] = dmn.prop[PhysicalProperyType::f_y];
  f[2] = dmn.prop[PhysicalProperyType::f_z];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  double wr = w*rho;

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  double ud[3] = {-f[0], -f[1], -f[2]};
  double u[3] = {0.0};
  double ux[3][3] = {0.0};
  double uxx[3][3][3] = {0.0};

  for (int a = 0; a < eNoNw; a++) {
    ud[0] = ud[0] + Nw(a)*(al(0,a)-bfl(0,a));
    ud[1] = ud[1] + Nw(a)*(al(1,a)-bfl(1,a));
    ud[2] = ud[2] + Nw(a)*(al(2,a)-bfl(2,a));

    u[0] = u[0] + Nw(a)*yl(0,a);
    u[1] = u[1] + Nw(a)*yl(1,a);
    u[2] = u[2] + Nw(a)*yl(2,a);

    ux[0][0] = ux[0][0] + Nwx(0,a)*yl(0,a);
    ux[1][0] = ux[1][0] + Nwx(1,a)*yl(0,a);
    ux[2][0] = ux[2][0] + Nwx(2,a)*yl(0,a);
    ux[0][1] = ux[0][1] + Nwx(0,a)*yl(1,a);
    ux[1][1] = ux[1][1] + Nwx(1,a)*yl(1,a);
    ux[2][1] = ux[2][1] + Nwx(2,a)*yl(1,a);
    ux[0][2] = ux[0][2] + Nwx(0,a)*yl(2,a);
    ux[1][2] = ux[1][2] + Nwx(1,a)*yl(2,a);
    ux[2][2] = ux[2][2] + Nwx(2,a)*yl(2,a);

    uxx[0][0][1] += Nwxx(0,a)*yl(0,a);
    uxx[1][0][1] += Nwxx(1,a)*yl(0,a);
    uxx[2][0][2] += Nwxx(2,a)*yl(0,a);
    uxx[1][0][0] += Nwxx(3,a)*yl(0,a);
    uxx[2][0][1] += Nwxx(4,a)*yl(0,a);
    uxx[0][0][2] += Nwxx(5,a)*yl(0,a);

    uxx[0][1][0] += Nwxx(0,a)*yl(1,a);
    uxx[1][1][1] += Nwxx(1,a)*yl(1,a);
    uxx[2][1][2] += Nwxx(2,a)*yl(1,a);
    uxx[1][1][0] += Nwxx(3,a)*yl(1,a);
    uxx[2][1][1] += Nwxx(4,a)*yl(1,a);
    uxx[0][1][2] += Nwxx(5,a)*yl(1,a);

    uxx[0][2][0] += Nwxx(0,a)*yl(2,a);
    uxx[1][2][1] += Nwxx(1,a)*yl(2,a);
    uxx[2][2][2] += Nwxx(2,a)*yl(2,a);
    uxx[1][2][0] += Nwxx(3,a)*yl(2,a);
    uxx[2][2][1] += Nwxx(4,a)*yl(2,a);
    uxx[0][2][2] += Nwxx(5,a)*yl(2,a);
  }

  double divU = ux[0][0] + ux[1][1] + ux[2][2];

  uxx[0][0][1] = uxx[1][0][0];
  uxx[1][0][2] = uxx[2][0][1];
  uxx[2][0][0] = uxx[0][0][2];

  uxx[0][1][1] = uxx[1][1][0];
  uxx[1][1][2] = uxx[2][1][1];
  uxx[2][1][0] = uxx[0][1][2];

  uxx[0][2][1] = uxx[1][2][0];
  uxx[1][2][2] = uxx[2][2][1];
  uxx[2][2][0] = uxx[0][2][2];

  double d2u2[3] = {0.0};
  d2u2[0] = uxx[0][0][0] + uxx[1][0][1] + uxx[2][0][2];
  d2u2[1] = uxx[0][1][0] + uxx[1][1][1] + uxx[2][1][2];
  d2u2[2] = uxx[0][2][0] + uxx[1][2][1] + uxx[2][2][2];

  // Pressure and its gradient
  //
  double px[3] = {0.0};

  for (int a = 0; a < eNoNq; a++) {
    px[0] = px[0] + Nqx(0,a)*yl(3,a);
    px[1] = px[1] + Nqx(1,a)*yl(3,a);
    px[2] = px[2] + Nqx(2,a)*yl(3,a);
  }

  // Update convection velocity relative to mesh velocity
  //
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u[0] = u[0] - Nw(a)*yl(4,a);
      u[1] = u[1] - Nw(a)*yl(5,a);
      u[2] = u[2] - Nw(a)*yl(6,a);
    }
  }

  // Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  //
  double es[3][3] = {0.0};
  es[0][0] = ux[0][0] + ux[0][0];
  es[1][1] = ux[1][1] + ux[1][1];
  es[2][2] = ux[2][2] + ux[2][2];
  es[1][0] = ux[1][0] + ux[0][1];
  es[2][1] = ux[2][1] + ux[1][2];
  es[0][2] = ux[0][2] + ux[2][0];
  es[0][1] = es[1][0];
  es[1][2] = es[2][1];
  es[2][0] = es[0][2];

  double esNx[3][MAX_SIZE];

  for (int a = 0; a < eNoNw; a++) {
    esNx[0][a] = es[0][0]*Nwx(0,a) + es[1][0]*Nwx(1,a) + es[2][0]*Nwx(2,a);
    esNx[1][a] = es[0][1]*Nwx(0,a) + es[1][1]*Nwx(1,a) + es[2][1]*Nwx(2,a);
    esNx[2][a] = es[0][2]*Nwx(0,a) + es[1][2]*Nwx(1,a) + es[2][2]*Nwx(2,a);
  }

  double es_x[3][3][3] = {0.0};

  for (int k = 0; k < 3; k++) { 
     es_x[0][0][k] = uxx[0][0][k] + uxx[0][0][k];
     es_x[1][1][k] = uxx[1][1][k] + uxx[1][1][k];
     es_x[2][2][k] = uxx[2][2][k] + uxx[2][2][k];
     es_x[1][0][k] = uxx[1][0][k] + uxx[0][1][k];
     es_x[2][1][k] = uxx[2][1][k] + uxx[1][2][k];
     es_x[0][2][k] = uxx[0][2][k] + uxx[2][0][k];

     es_x[0][1][k] = es_x[1][0][k];
     es_x[1][2][k] = es_x[2][1][k];
     es_x[2][0][k] = es_x[0][2][k];
  }

  double mu_x[3];

  mu_x[0] = (es_x[0][0][0]*es[0][0] + es_x[1][1][0]*es[1][1]
          +  es_x[2][2][0]*es[2][2])*0.5
          +  es_x[1][0][0]*es[1][0] + es_x[2][1][0]*es[2][1]
          +  es_x[0][2][0]*es[0][2];

  mu_x[1] = (es_x[0][0][1]*es[0][0] + es_x[1][1][1]*es[1][1]
          +  es_x[2][2][1]*es[2][2])*0.5
          +  es_x[1][0][1]*es[1][0] + es_x[2][1][1]*es[2][1]
          +  es_x[0][2][1]*es[0][2];

  mu_x[2] = (es_x[0][0][2]*es[0][0] + es_x[1][1][2]*es[1][1]
          +  es_x[2][2][2]*es[2][2])*0.5
          +  es_x[1][0][2]*es[1][0] + es_x[2][1][2]*es[2][1]
          +  es_x[0][2][2]*es[0][2];

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es[0][0]*es[0][0] + es[1][0]*es[1][0] + es[2][0]*es[2][0]
             + es[0][1]*es[0][1] + es[1][1]*es[1][1] + es[2][1]*es[2][1]
             + es[0][2]*es[0][2] + es[1][2]*es[1][2] + es[2][2]*es[2][2];
  gam = sqrt(0.5*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     mu_g = mu_g / gam;
  }

  for (int i = 0; i < 3; i++) {
    mu_x[i] = mu_g * mu_x[i];
  }
  //std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
  //mu_x(:) = mu_g * mu_x(:)

  // Stabilization parameters
  //
  double up[3] = {0.0};
  double updu[3][3][MAX_SIZE] = {0.0};
  double tauM = 0.0;

  if (vmsFlag) {
    double kT = 4.0 * pow(ctM/dt,2.0);

    double kU = u[0]*u[0]*Kxi(0,0) + u[1]*u[0]*Kxi(1,0) + u[2]*u[0]*Kxi(2,0)
              + u[0]*u[1]*Kxi(0,1) + u[1]*u[1]*Kxi(1,1) + u[2]*u[1]*Kxi(2,1)
              + u[0]*u[2]*Kxi(0,2) + u[1]*u[2]*Kxi(1,2) + u[2]*u[2]*Kxi(2,2);

    double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(2,0)*Kxi(2,0)
              + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
              + Kxi(0,2)*Kxi(0,2) + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2);

    kS = ctC * kS * pow(mu/rho,2.0);
    tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

    double rV[3];
    rV[0] = ud[0] + u[0]*ux[0][0] + u[1]*ux[1][0] + u[2]*ux[2][0];
    rV[1] = ud[1] + u[0]*ux[0][1] + u[1]*ux[1][1] + u[2]*ux[2][1];
    rV[2] = ud[2] + u[0]*ux[0][2] + u[1]*ux[1][2] + u[2]*ux[2][2];

    double rS[3];
    rS[0] = mu_x[0]*es[0][0] + mu_x[1]*es[1][0] + mu_x[2]*es[2][0] + mu*d2u2[0];
    rS[1] = mu_x[0]*es[0][1] + mu_x[1]*es[1][1] + mu_x[2]*es[2][1] + mu*d2u2[1];
    rS[2] = mu_x[0]*es[0][2] + mu_x[1]*es[1][2] + mu_x[2]*es[2][2] + mu*d2u2[2];

    up[0] = -tauM*(rho*rV[0] + px[0] - rS[0]);
    up[1] = -tauM*(rho*rV[1] + px[1] - rS[1]);
    up[2] = -tauM*(rho*rV[2] + px[2] - rS[2]);

    for (int a = 0; a < eNoNw; a++) {
      double uNx = u[0]*Nwx(0,a) + u[1]*Nwx(1,a) + u[2]*Nwx(2,a);
      T1 = -rho*uNx + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a);

      updu[0][0][a] = mu_x[0]*Nwx(0,a) + d2u2[0]*mu_g*esNx[0][a] + T1;
      updu[1][0][a] = mu_x[1]*Nwx(0,a) + d2u2[0]*mu_g*esNx[1][a];
      updu[2][0][a] = mu_x[2]*Nwx(0,a) + d2u2[0]*mu_g*esNx[2][a];
  
      updu[0][1][a] = mu_x[0]*Nwx(1,a) + d2u2[1]*mu_g*esNx[0][a];
      updu[1][1][a] = mu_x[1]*Nwx(1,a) + d2u2[1]*mu_g*esNx[1][a] + T1;
      updu[2][1][a] = mu_x[2]*Nwx(1,a) + d2u2[1]*mu_g*esNx[2][a];
  
      updu[0][2][a] = mu_x[0]*Nwx(2,a) + d2u2[2]*mu_g*esNx[0][a];
      updu[1][2][a] = mu_x[1]*Nwx(2,a) + d2u2[2]*mu_g*esNx[1][a];
      updu[2][2][a] = mu_x[2]*Nwx(2,a) + d2u2[2]*mu_g*esNx[2][a] + T1;
    }

  } else {
    tauM = 0.0;
    std::memset(up, 0, sizeof up);
    std::memset(updu, 0, sizeof updu);
  }

  //  Local residue
  //
  for (int a = 0; a < eNoNq; a++) {
    double upNx = up[0]*Nqx(0,a) + up[1]*Nqx(1,a) + up[2]*Nqx(2,a);
    lR(3,a) = lR(3,a) + w*(Nq(a)*divU - upNx);
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    T1 = rho*amd*Nw(b);

    for (int a = 0; a < eNoNq; a++) {
      // dRc_a/dU_b1
      double T2 = Nqx(0,a)*(updu[0][0][b] - T1) + Nqx(1,a)*updu[0][1][b] + Nqx(2,a)*updu[0][2][b];
      lK(12,a,b) = lK(12,a,b) + wl*(Nq(a)*Nwx(0,b) - tauM*T2);

      // dRc_a/dU_b2
      T2 = Nqx(0,a)*updu[1][0][b] + Nqx(1,a)*(updu[1][1][b] - T1) + Nqx(2,a)*updu[1][2][b];
      lK(13,a,b) = lK(13,a,b) + wl*(Nq(a)*Nwx(1,b) - tauM*T2);

      // dRc_a/dU_b3
      T2 = Nqx(0,a)*updu[2][0][b] + Nqx(1,a)*updu[2][1][b] + Nqx(2,a)*(updu[2][2][b] - T1);
      lK(14,a,b) = lK(14,a,b) + wl*(Nq(a)*Nwx(2,b) - tauM*T2);
    }
  }

  if (vmsFlag) {
    for (int b = 0; b < eNoNq; b++) {
      for (int a = 0; a < eNoNq; a++) {
        // dC/dP
        double NxNx = Nqx(0,a)*Nqx(0,b) + Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b);
        lK(15,a,b) = lK(15,a,b) + wl*tauM*NxNx;
      }
    }
  }
}

//------------
// fluid_3d_m
//------------
// Element momentum residue.
//
//  Modifies:
//    lR(dof,eNoN)  - Residue
//    lK(dof*dof,eNoN,eNoN) - Stiffness matrix
//
void fluid_3d_m(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w,
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx,
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_fluid_3d_m
  #ifdef debug_fluid_3d_m
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  double start_time = utils::cput();
  #endif

  // Maximum size of arrays sized by (3,eNoNw) -> (3,MAX_SIZE).
  const int MAX_SIZE = 8;

  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double ctM  = 1.0;
  double ctC  = 36.0;

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  std::array<double,3> f;
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  f[1] = dmn.prop[PhysicalProperyType::f_y];
  f[2] = dmn.prop[PhysicalProperyType::f_z];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  double wr = w*rho;

  #ifdef debug_fluid_3d_m
  dmsg << "rho: " << rho;
  dmsg << "T1: " << T1;
  dmsg << "wl: " << wl;
  dmsg << "wr: " << wr;
  #endif

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  std::array<double,3> ud{-f[0], -f[1], -f[2]};
  //ud  = -f
  double u[3] = {0.0};
  double ux[3][3] = {0.0};
  double uxx[3][3][3] = {0.0};

  for (int a = 0; a < eNoNw; a++) {
    ud[0] = ud[0] + Nw(a)*(al(0,a)-bfl(0,a));
    ud[1] = ud[1] + Nw(a)*(al(1,a)-bfl(1,a));
    ud[2] = ud[2] + Nw(a)*(al(2,a)-bfl(2,a));

    u[0] = u[0] + Nw(a)*yl(0,a);
    u[1] = u[1] + Nw(a)*yl(1,a);
    u[2] = u[2] + Nw(a)*yl(2,a);
    
    ux[0][0] += Nwx(0,a)*yl(0,a);
    ux[1][0] += Nwx(1,a)*yl(0,a);
    ux[2][0] += Nwx(2,a)*yl(0,a);
    ux[0][1] += Nwx(0,a)*yl(1,a);
    ux[1][1] += Nwx(1,a)*yl(1,a);
    ux[2][1] += Nwx(2,a)*yl(1,a);
    ux[0][2] += Nwx(0,a)*yl(2,a);
    ux[1][2] += Nwx(1,a)*yl(2,a);
    ux[2][2] += Nwx(2,a)*yl(2,a);

    uxx[0][0][0] += Nwxx(0,a)*yl(0,a);
    uxx[1][0][1] += Nwxx(1,a)*yl(0,a);
    uxx[2][0][2] += Nwxx(2,a)*yl(0,a);
    uxx[1][0][0] += Nwxx(3,a)*yl(0,a);
    uxx[2][0][1] += Nwxx(4,a)*yl(0,a);
    uxx[0][0][2] += Nwxx(5,a)*yl(0,a);

    uxx[0][1][0] += Nwxx(0,a)*yl(1,a);
    uxx[1][1][1] += Nwxx(1,a)*yl(1,a);
    uxx[2][1][2] += Nwxx(2,a)*yl(1,a);
    uxx[1][1][0] += Nwxx(3,a)*yl(1,a);
    uxx[2][1][1] += Nwxx(4,a)*yl(1,a);
    uxx[0][1][2] += Nwxx(5,a)*yl(1,a);

    uxx[0][2][0] += Nwxx(0,a)*yl(2,a);
    uxx[1][2][1] += Nwxx(1,a)*yl(2,a);
    uxx[2][2][2] += Nwxx(2,a)*yl(2,a);
    uxx[1][2][0] += Nwxx(3,a)*yl(2,a);
    uxx[2][2][1] += Nwxx(4,a)*yl(2,a);
    uxx[0][2][2] += Nwxx(5,a)*yl(2,a);
  }

  double divU = ux[0][0] + ux[1][1] + ux[2][2];
  #ifdef debug_fluid_3d_m
  dmsg << "divU: " << divU;
  #endif

  uxx[0][0][1] = uxx[1][0][0];
  uxx[1][0][2] = uxx[2][0][1];
  uxx[2][0][0] = uxx[0][0][2];

  uxx[0][1][1] = uxx[1][1][0];
  uxx[1][1][2] = uxx[2][1][1];
  uxx[2][1][0] = uxx[0][1][2];

  uxx[0][2][1] = uxx[1][2][0];
  uxx[1][2][2] = uxx[2][2][1];
  uxx[2][2][0] = uxx[0][2][2];

  std::array<double,3> d2u2{0.0};
  d2u2[0] = uxx[0][0][0] + uxx[1][0][1] + uxx[2][0][2];
  d2u2[1] = uxx[0][1][0] + uxx[1][1][1] + uxx[2][1][2];
  d2u2[2] = uxx[0][2][0] + uxx[1][2][1] + uxx[2][2][2];

  // Pressure and its gradient
  //
  double p = 0.0;
  double px[3] = {0.0};

  for (int a = 0; a < eNoNq; a++) {
    p  = p + Nq(a)*yl(3,a);
    px[0] = px[0] + Nqx(0,a)*yl(3,a);
    px[1] = px[1] + Nqx(1,a)*yl(3,a);
    px[2] = px[2] + Nqx(2,a)*yl(3,a);
  }
  #ifdef debug_fluid_3d_m
  dmsg << "p: " << p;
  dmsg << "px: " << px[0] << " " << px[1] << " " << px[2];
  #endif

  // Update convection velocity relative to mesh velocity
  //
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u[0] = u[0] - Nw(a)*yl(4,a);
      u[1] = u[1] - Nw(a)*yl(5,a);
      u[2] = u[2] - Nw(a)*yl(6,a);
    }
  }

  // Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  //
  double es[3][3] = {0.0};
  es[0][0] = ux[0][0] + ux[0][0];
  es[1][1] = ux[1][1] + ux[1][1];
  es[2][2] = ux[2][2] + ux[2][2];
  es[1][0] = ux[1][0] + ux[0][1];
  es[2][1] = ux[2][1] + ux[1][2];
  es[0][2] = ux[0][2] + ux[2][0];
  es[0][1] = es[1][0];
  es[1][2] = es[2][1];
  es[2][0] = es[0][2];

  double esNx[3][MAX_SIZE];

  for (int a = 0; a < eNoNw; a++) {
    esNx[0][a] = es[0][0]*Nwx(0,a) + es[1][0]*Nwx(1,a) + es[2][0]*Nwx(2,a);
    esNx[1][a] = es[0][1]*Nwx(0,a) + es[1][1]*Nwx(1,a) + es[2][1]*Nwx(2,a);
    esNx[2][a] = es[0][2]*Nwx(0,a) + es[1][2]*Nwx(1,a) + es[2][2]*Nwx(2,a);
  }

  double es_x[3][3][3];

  for (int k = 0; k < 3; k++) { 
     es_x[0][0][k] = uxx[0][0][k] + uxx[0][0][k];
     es_x[1][1][k] = uxx[1][1][k] + uxx[1][1][k];
     es_x[2][2][k] = uxx[2][2][k] + uxx[2][2][k];
     es_x[1][0][k] = uxx[1][0][k] + uxx[0][1][k];
     es_x[2][1][k] = uxx[2][1][k] + uxx[1][2][k];
     es_x[0][2][k] = uxx[0][2][k] + uxx[2][0][k];

     es_x[0][1][k] = es_x[1][0][k];
     es_x[1][2][k] = es_x[2][1][k];
     es_x[2][0][k] = es_x[0][2][k];
  }

  std::array<double,3> mu_x{0.0};

  mu_x[0] = (es_x[0][0][0]*es[0][0] + es_x[1][1][0]*es[1][1]
          +  es_x[2][2][0]*es[2][2])*0.5
          +  es_x[1][0][0]*es[1][0] + es_x[2][1][0]*es[2][1]
          +  es_x[0][2][0]*es[0][2];

  mu_x[1] = (es_x[0][0][1]*es[0][0] + es_x[1][1][1]*es[1][1]
          +  es_x[2][2][1]*es[2][2])*0.5
          +  es_x[1][0][1]*es[1][0] + es_x[2][1][1]*es[2][1]
          +  es_x[0][2][1]*es[0][2];

  mu_x[2] = (es_x[0][0][2]*es[0][0] + es_x[1][1][2]*es[1][1]
          +  es_x[2][2][2]*es[2][2])*0.5
          +  es_x[1][0][2]*es[1][0] + es_x[2][1][2]*es[2][1]
          +  es_x[0][2][2]*es[0][2];

  #ifdef debug_fluid_3d_m
  dmsg << "mu_x: " << mu_x[0] << " " << mu_x[1] << " " << mu_x[2];
  #endif

  // Shear-rate := (2*e_ij*e_ij)^.5
  //
  //dmsg << "Compute shear rate ... ";
  double gam = es[0][0]*es[0][0] + es[1][0]*es[1][0] + es[2][0]*es[2][0]
             + es[0][1]*es[0][1] + es[1][1]*es[1][1] + es[2][1]*es[2][1]
             + es[0][2]*es[0][2] + es[1][2]*es[1][2] + es[2][2]*es[2][2];
  gam = sqrt(0.5*gam);
  #ifdef debug_fluid_3d_m
  dmsg << "gam: " << gam;
  #endif

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  //
  double mu, mu_s, mu_g;
  fluid::get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     mu_g = mu_g / gam;
  }
  std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
  //mu_x(:) = mu_g * mu_x(:)

  // Stabilization parameters
  //
  double kT = 4.0 * pow(ctM/dt,2.0);

  double kU = u[0]*u[0]*Kxi(0,0) + u[1]*u[0]*Kxi(1,0) + u[2]*u[0]*Kxi(2,0)
            + u[0]*u[1]*Kxi(0,1) + u[1]*u[1]*Kxi(1,1) + u[2]*u[1]*Kxi(2,1)
            + u[0]*u[2]*Kxi(0,2) + u[1]*u[2]*Kxi(1,2) + u[2]*u[2]*Kxi(2,2);

  double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(2,0)*Kxi(2,0)
            + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
            + Kxi(0,2)*Kxi(0,2) + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2);
  kS = ctC * kS * pow(mu/rho,2.0);
  double tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

  #ifdef debug_fluid_3d_m
  dmsg << "kT: " << kT;
  dmsg << "kU: " << kU;
  dmsg << "kS: " << kS;
  dmsg << "tauM: " << tauM;
  #endif

  double rV[3] = {0.0};
  rV[0] = ud[0] + u[0]*ux[0][0] + u[1]*ux[1][0] + u[2]*ux[2][0];
  rV[1] = ud[1] + u[0]*ux[0][1] + u[1]*ux[1][1] + u[2]*ux[2][1];
  rV[2] = ud[2] + u[0]*ux[0][2] + u[1]*ux[1][2] + u[2]*ux[2][2];

  double rS[3] = {0.0};
  rS[0] = mu_x[0]*es[0][0] + mu_x[1]*es[1][0] + mu_x[2]*es[2][0] + mu*d2u2[0];
  rS[1] = mu_x[0]*es[0][1] + mu_x[1]*es[1][1] + mu_x[2]*es[2][1] + mu*d2u2[1];
  rS[2] = mu_x[0]*es[0][2] + mu_x[1]*es[1][2] + mu_x[2]*es[2][2] + mu*d2u2[2];

  double up[3] = {0.0};
  up[0] = -tauM*(rho*rV[0] + px[0] - rS[0]);
  up[1] = -tauM*(rho*rV[1] + px[1] - rS[1]);
  up[2] = -tauM*(rho*rV[2] + px[2] - rS[2]);

  double tauC, tauB, pa;
  double eps = std::numeric_limits<double>::epsilon();
  double ua[3] = {0.0};

  if (vmsFlag) {
    tauC = 1.0 / (tauM * (Kxi(0,0) + Kxi(1,1) + Kxi(2,2)));
    tauB = up[0]*up[0]*Kxi(0,0) + up[1]*up[0]*Kxi(1,0)
         + up[2]*up[0]*Kxi(2,0) + up[0]*up[1]*Kxi(0,1)
         + up[1]*up[1]*Kxi(1,1) + up[2]*up[1]*Kxi(2,1)
         + up[0]*up[2]*Kxi(0,2) + up[1]*up[2]*Kxi(1,2)
         + up[2]*up[2]*Kxi(2,2);

    if (utils::is_zero(tauB)) {
      tauB = eps;
    }
    tauB = rho / sqrt(tauB);

    ua[0] = u[0] + up[0];
    ua[1] = u[1] + up[1];
    ua[2] = u[2] + up[2];
    pa = p - tauC*divU;

  } else {
    tauC = 0.0;
    tauB = 0.0;
    for (int i = 0; i < 3; i++) {
      ua[i] = u[i];
    }
    pa = p;
  }

  rV[0] = tauB*(up[0]*ux[0][0] + up[1]*ux[1][0] + up[2]*ux[2][0]);
  rV[1] = tauB*(up[0]*ux[0][1] + up[1]*ux[1][1] + up[2]*ux[2][1]);
  rV[2] = tauB*(up[0]*ux[0][2] + up[1]*ux[1][2] + up[2]*ux[2][2]);

  double rM[3][3];
  rM[0][0] = mu*es[0][0] - rho*up[0]*ua[0] + rV[0]*up[0] - pa;
  rM[1][0] = mu*es[1][0] - rho*up[0]*ua[1] + rV[0]*up[1];
  rM[2][0] = mu*es[2][0] - rho*up[0]*ua[2] + rV[0]*up[2];

  rM[0][1] = mu*es[0][1] - rho*up[1]*ua[0] + rV[1]*up[0];
  rM[1][1] = mu*es[1][1] - rho*up[1]*ua[1] + rV[1]*up[1] - pa;
  rM[2][1] = mu*es[2][1] - rho*up[1]*ua[2] + rV[1]*up[2];

  rM[0][2] = mu*es[0][2] - rho*up[2]*ua[0] + rV[2]*up[0];
  rM[1][2] = mu*es[1][2] - rho*up[2]*ua[1] + rV[2]*up[1];
  rM[2][2] = mu*es[2][2] - rho*up[2]*ua[2] + rV[2]*up[2] - pa;

  rV[0] = ud[0] + ua[0]*ux[0][0] + ua[1]*ux[1][0] + ua[2]*ux[2][0];
  rV[1] = ud[1] + ua[0]*ux[0][1] + ua[1]*ux[1][1] + ua[2]*ux[2][1];
  rV[2] = ud[2] + ua[0]*ux[0][2] + ua[1]*ux[1][2] + ua[2]*ux[2][2];

  //  Local residue
  //
  double updu[3][3][MAX_SIZE] = {0.0};
  double uNx[MAX_SIZE] = {0.0};
  double upNx[MAX_SIZE] = {0.0}; 
  double uaNx[MAX_SIZE] = {0.0}; 

  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + wr*Nw(a)*rV[0] + w*(Nwx(0,a)*rM[0][0] + Nwx(1,a)*rM[1][0] + Nwx(2,a)*rM[2][0]);
    lR(1,a) = lR(1,a) + wr*Nw(a)*rV[1] + w*(Nwx(0,a)*rM[0][1] + Nwx(1,a)*rM[1][1] + Nwx(2,a)*rM[2][1]);
    lR(2,a) = lR(2,a) + wr*Nw(a)*rV[2] + w*(Nwx(0,a)*rM[0][2] + Nwx(1,a)*rM[1][2] + Nwx(2,a)*rM[2][2]);

    // Quantities used for stiffness matrix
    uNx[a]  = u[0]*Nwx(0,a)  + u[1]*Nwx(1,a)  + u[2]*Nwx(2,a);
    upNx[a] = up[0]*Nwx(0,a) + up[1]*Nwx(1,a) + up[2]*Nwx(2,a);

    if (vmsFlag) {
       uaNx[a] = uNx[a] + upNx[a];
    } else {
       uaNx[a] = uNx[a];
    }

    T1 = -rho*uNx[a] + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a);

    updu[0][0][a] = mu_x[0]*Nwx(0,a) + d2u2[0]*mu_g*esNx[0][a] + T1;
    updu[1][0][a] = mu_x[1]*Nwx(0,a) + d2u2[0]*mu_g*esNx[1][a];
    updu[2][0][a] = mu_x[2]*Nwx(0,a) + d2u2[0]*mu_g*esNx[2][a];

    updu[0][1][a] = mu_x[0]*Nwx(1,a) + d2u2[1]*mu_g*esNx[0][a];
    updu[1][1][a] = mu_x[1]*Nwx(1,a) + d2u2[1]*mu_g*esNx[1][a] + T1;
    updu[2][1][a] = mu_x[2]*Nwx(1,a) + d2u2[1]*mu_g*esNx[2][a];

    updu[0][2][a] = mu_x[0]*Nwx(2,a) + d2u2[2]*mu_g*esNx[0][a];
    updu[1][2][a] = mu_x[1]*Nwx(2,a) + d2u2[2]*mu_g*esNx[1][a];
    updu[2][2][a] = mu_x[2]*Nwx(2,a) + d2u2[2]*mu_g*esNx[2][a] + T1;
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      rM[0][0] = Nwx(0,a)*Nwx(0,b);
      rM[1][0] = Nwx(1,a)*Nwx(0,b);
      rM[2][0] = Nwx(2,a)*Nwx(0,b);
      rM[0][1] = Nwx(0,a)*Nwx(1,b);
      rM[1][1] = Nwx(1,a)*Nwx(1,b);
      rM[2][1] = Nwx(2,a)*Nwx(1,b);
      rM[0][2] = Nwx(0,a)*Nwx(2,b);
      rM[1][2] = Nwx(1,a)*Nwx(2,b);
      rM[2][2] = Nwx(2,a)*Nwx(2,b);

      double NxNx = Nwx(0,a)*Nwx(0,b) + Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b);
      T1 = mu*NxNx + rho*amd*Nw(b)*(Nw(a) + rho*tauM*uaNx[a]) + rho*Nw(a)*(uNx[b]+upNx[b]) + tauB*upNx[a]*upNx[b];

      // dRm_a1/du_b1
      double T2 = (mu + tauC)*rM[0][0] + esNx[0][a]*mu_g*esNx[0][b] - rho*tauM*uaNx[a]*updu[0][0][b];
      lK(0,a,b)  = lK(0,a,b)  + wl*(T2 + T1);

      // dRm_a1/du_b2
      T2 = mu*rM[1][0] + tauC*rM[0][1] + esNx[0][a]*mu_g*esNx[1][b] - rho*tauM*uaNx[a]*updu[1][0][b];
      lK(1,a,b)  = lK(1,a,b)  + wl*(T2);

      // dRm_a1/du_b3
      T2 = mu*rM[2][0] + tauC*rM[0][2] + esNx[0][a]*mu_g*esNx[2][b] - rho*tauM*uaNx[a]*updu[2][0][b];
      lK(2,a,b)  = lK(2,a,b)  + wl*(T2);

      // dRm_a2/du_b1
      T2 = mu*rM[0][1] + tauC*rM[1][0] + esNx[1][a]*mu_g*esNx[0][b] - rho*tauM*uaNx[a]*updu[0][1][b];
      lK(4,a,b) = lK(4,a,b)  + wl*(T2);

      // dRm_a2/du_b2
      T2 = (mu + tauC)*rM[1][1] + esNx[1][a]*mu_g*esNx[1][b] - rho*tauM*uaNx[a]*updu[1][1][b];
      lK(5,a,b)  = lK(5,a,b)  + wl*(T2 + T1);

      // dRm_a2/du_b3
      T2 = mu*rM[2][1] + tauC*rM[1][2] + esNx[1][a]*mu_g*esNx[2][b] - rho*tauM*uaNx[a]*updu[2][1][b];
      lK(6,a,b)  = lK(6,a,b)  + wl*(T2);

      // dRm_a3/du_b1
      T2 = mu*rM[0][2] + tauC*rM[2][0] + esNx[2][a]*mu_g*esNx[0][b] - rho*tauM*uaNx[a]*updu[0][2][b];
      lK(8,a,b)  = lK(8,a,b)  + wl*(T2);

      // dRm_a3/du_b2
      T2 = mu*rM[1][2] + tauC*rM[2][1] + esNx[2][a]*mu_g*esNx[1][b] - rho*tauM*uaNx[a]*updu[1][2][b];
      lK(9,a,b) = lK(9,a,b) + wl*(T2);

      // dRm_a3/du_b3;
      T2 = (mu + tauC)*rM[2][2] + esNx[2][a]*mu_g*esNx[2][b] - rho*tauM*uaNx[a]*updu[2][2][b];
      lK(10,a,b) = lK(10,a,b) + wl*(T2 + T1);
      //dmsg << "lK(10,a,b): " << lK(10,a,b);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      T1 = rho*tauM*uaNx[a];

      // dRm_a1/dp_b
      lK(3,a,b)  = lK(3,a,b)  - wl*(Nwx(0,a)*Nq(b) - Nqx(0,b)*T1);

      // dRm_a2/dp_b
      lK(7,a,b)  = lK(7,a,b)  - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T1);

      // dRm_a3/dp_b
      lK(11,a,b) = lK(11,a,b) - wl*(Nwx(2,a)*Nq(b) - Nqx(2,b)*T1);
    }
  }
}
//---------------
// get_viscosity
//---------------
//
void get_viscosity(const ComMod& com_mod, const dmnType& lDmn, double& gamma, double& mu, double& mu_s, double& mu_x)
{
  using namespace consts;

  mu = 0.0;
  mu_s = 0.0;
  mu_x = 0.0;

  double mu_i, mu_o, lam, a, n, T1, T2;

  switch (lDmn.visc.viscType) {

    case FluidViscosityModelType::viscType_Const:
      mu = lDmn.visc.mu_i;
      mu_s = mu;
      mu_x = 0.0;
    break;

    case FluidViscosityModelType::viscType_CY:
      mu_i = lDmn.visc.mu_i;
      mu_o = lDmn.visc.mu_o;
      lam = lDmn.visc.lam;
      a = lDmn.visc.a;
      n = lDmn.visc.n;

      T1 = 1.0 + pow(lam*gamma, a);
      T2 = pow(T1,((n-1.0)/a));
      mu = mu_i + (mu_o - mu_i)*T2;
      mu_s = mu_i;

      T1 = T2 / T1;
      T2 = pow(lam,a) * pow(gamma,(a-1.0)) * T1;
      mu_x = (mu_o - mu_i) * (n - 1.0) * T2;
    break;

    case FluidViscosityModelType::viscType_Cass:
      mu_i = lDmn.visc.mu_i;
      mu_o = lDmn.visc.mu_o;
      lam  = lDmn.visc.lam;

      if (gamma < lam) { 
         mu_o = mu_o / sqrt(lam);
         gamma = lam;
      } else { 
         mu_o = mu_o / sqrt(gamma);
      }

      mu  = (mu_i + mu_o) * (mu_i + mu_o);
      mu_s = mu_i * mu_i;
      mu_x = 2.0 * mu_o * (mu_o + mu_i) / gamma;
    break;
  } 
}

};
