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
