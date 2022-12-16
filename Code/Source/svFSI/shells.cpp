
// This routines is for solving nonlinear shell mechanics problem
// using linear triangle finite elements and IGA.

#include "shells.h"

#include "nn.h"

namespace shells {

//----------
// shell_fp
//----------
//
void shell_fp(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& dl, const Array<double>& xl, const Vector<double>& tfl, Array<double>& lR, Array3<double>& lK)
{
  int nsd = com_mod.nsd;
  int dof = com_mod.dof;
  double dt = com_mod.dt;
  auto& cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];

  double afl = eq.af * eq.beta*dt*dt;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  // Get the current configuration and traction vector
  //
  Array<double> xc(3,eNoN);
  double tfn = 0.0;

  for (int a = 0; a < eNoN; a++) {
    xc(0,a) = xl(0,a) + dl(i,a);
    xc(1,a) = xl(1,a) + dl(j,a);
    xc(2,a) = xl(2,a) + dl(k,a);
    tfn = tfn + N(a)*tfl(a);
  }

  double wl = w * tfn;

  // Covariant and contravariant bases in current config
  Vector<double> nV(3); 
  Array<double> gCov(3,2);
  Array<double> gCnv(3,2);
  nn::gnns(nsd, eNoN, Nx, xc, nV, gCov, gCnv);
  //CALL GNNS(eNoN, Nx, xc, nV, gCov, gCnv)

  // Local residue
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) - wl*N(a)*nV(0);
    lR(1,a) = lR(1,a) - wl*N(a)*nV(1);
    lR(2,a) = lR(2,a) - wl*N(a)*nV(2);
  }

  // Local stiffness: mass matrix and stiffness contribution due to
  // follower traction load
  //
  double T1 = afl*wl*0.5;

  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      auto lKp = gCov.col(0) * (N(b)*Nx(1,a) - N(a)*Nx(2,b)) - gCov.col(1)*(N(b)*Nx(0,a) - N(a)*Nx(0,b));
      //lKp(:) = gCov(:,1)*(N(b)*Nx(2,a) - N(a)*Nx(2,b)) - gCov(:,2)*(N(b)*Nx(1,a) - N(a)*Nx(1,b))

      lK(1,a,b) = lK(1,a,b) - T1*lKp(2);
      lK(2,a,b) = lK(2,a,b) + T1*lKp(1);

      lK(dof+1,a,b) = lK(dof+1,a,b) + T1*lKp(3);
      lK(dof+3,a,b) = lK(dof+3,a,b) - T1*lKp(1);

      lK(2*dof,a,b) = lK(2*dof,a,b) - T1*lKp(1);
      //lK(2*dof+1,a,b) = lK(2*dof+1,a,b) - T1*lKp(2)

      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + T1*lKp(0);
      //lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + T1*lKp(1)
    }
  }
}


};

