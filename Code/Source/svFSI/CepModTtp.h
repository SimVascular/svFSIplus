
// The classes defined here duplicate the data structures in the Fortran TPPMOD module defined 
// in CEPMOD_TTP.f and PARAMS_TPP.f files. 

// This module defines data structures for ten Tusscher-Panfilov
// epicardial cellular activation model for cardiac electrophysiology

#ifndef CEP_MOD_TTP_H 
#define CEP_MOD_TTP_H 

#include "Array.h"
#include "Vector.h"
#include <array>

#include <optional>
#include <functional>

template <class T>
T& make_ref(T&& x) { return x; }

class CepModTtp
{
  public:
    CepModTtp();
    ~CepModTtp();


//--------------------------------------------------------------------
//
//     Constants for TenTusscher-Panfilov Ventricular Myocyte Model.
//
//--------------------------------------------------------------------

//     Default model parameters
//     R: Gas constant
      double Rc = 8314.472;     // units: J/mol/K

//     T: Temperature
      double Tc = 310.0;          // units: K

//     F: Faraday constant
      double Fc = 96485.3415;    // units: C/mmol

//     Cm: Cell capacitance per unit surface area
      double Cm = 0.185;         // units: uF/cm^{2}

//     sV: Surface to volume ratio
      double sV = 0.2;           // units: um^{-1}

//     rho: Cellular resistivity
      double rho = 162.0;         // units: \Omega-cm

//     V_c: Cytoplasmic volume
      double V_c = 16.404E-3;    // units: um^{3}

//     V_sr: Sacroplasmic reticulum volume
      double V_sr = 1.094E-3;    // units: um^{3}

//     V_ss: Subspace volume
      double V_ss = 5.468E-5;    // units: um^{3}

//     K_o: Extracellular K concentration
      double K_o = 5.4;          // units: mM

//     Na_o: Extracellular Na concentration
      double Na_o = 140.0;        // units: mM

//     Ca_o: Extracellular Ca concentration
      double Ca_o = 2.0;          // units: mM

//     G_Na: Maximal I_Na conductance
      double G_Na = 14.838;      // units: nS/pF

//     G_K1: Maximal I_K1 conductance
      double G_K1 = 5.405;       // units: nS/pF

//     G_to: Maximal epicardial I_to conductance, units: nS/pF
      Vector<double> G_to = {0.294, 0.073, 0.294};

//     G_Kr: Maximal I_Kr conductance
      double G_Kr = 0.153;      // units: nS/pF

//     G_Kr for spiral wave breakup
//      double G_Kr = 0.172;     // units: nS/pF

//     G_Ks: Maximal epicardial I_Ks conductance, units: nS/pF
      Vector<double> G_Ks = {0.392, 0.392, 0.098};

//     G_Ks for spiral wave breakup (epi)
//      double G_Ks(3) = (/0.441, 0.392_RKIND, 0.098_RKIND/)

//     p_KNa: Relative I_Ks permeability to Na
      double p_KNa = 3.E-2;     // dimensionless

//     G_CaL: Maximal I_CaL conductance
      double G_CaL = 3.98E-5;   // units: cm^{3}/uF/ms

//     K_NaCa: Maximal I_NaCa
      double K_NaCa = 1000.;    // units: pA/pF

//     gamma: Voltage dependent parameter of I_NaCa
      double gamma = 0.35;      // dimensionless

//     K_mCa: Ca_i half-saturation constant for I_NaCa
      double K_mCa = 1.38;      // units: mM

//     K_mNai: Na_i half-saturation constant for I_NaCa
      double K_mNai = 87.5;     // units: mM

//     K_sat: Saturation factor for I_NaCa
      double K_sat = 0.1;       // dimensionless

//     alpha: Factor enhancing outward nature of I_NaCa
      double alpha = 2.5;       // dimensionless

//     p_NaK: Maximal I_NaK
      double p_NaK = 2.724;     // units: pA/pF

//     K_mK: K_o half-saturation constant of I_NaK
      double K_mK = 1.;         // units: mM

//     K_mNa: Na_i half-saturation constant of I_NaK
      double K_mNa = 40.;       // units: mM

//     G_pK: Maximal I_pK conductance
      double G_pK = 1.46E-2;    // units: nS/pF

//     G_pK for spiral wave breakup
//      double G_pK = 2.19E-3;    // units: nS/pF

//     G_pCa: Maximal I_pCa conductance
      double G_pCa = 0.1238;    // units: pA/pF

//     G_pCa for spiral wave breakup
//      double G_pCa = 0.8666;    // units: pA/pF

//     K_pCa: Half-saturation constant of I_pCa
      double K_pCa = 5.E-4;     // units: mM

//     G_bNa: Maximal I_bNa conductance
      double G_bNa = 2.9E-4;    // units: nS/pF

//     G_bCa: Maximal I_bCa conductance
      double G_bCa = 5.92E-4;   // units: nS/pF

//     Vmax_up: Maximal I_up conductance
      double Vmax_up = 6.375E-3;// units: mM/ms

//     K_up: Half-saturation constant of I_up
      double K_up = 2.5E-4;     // units: mM

//     V_rel: Maximal I_rel conductance
      double V_rel = 0.102;     // units: mM/ms

//     k1p: R to O and RI to I, I_rel transition rate
      double k1p = 0.15;         // units: mM^{-2}/ms

//     k2p: O to I and R to RI, I_rel transition rate
      double k2p = 4.5E-2;       // units: mM^{-1}/ms

//     k3: O to R and I to RI, I_rel transition rate
      double k3 = 6.E-2;         // units: ms^{-1}

//     k4: I to O and Ri to I, I_rel transition rate
      double k4 = 5.E-3;         // units: ms^{-1}

//     EC: Ca_sr half-saturation constant of k_casr
      double EC = 1.5;           // units: mM

//     max_sr: Maximum value of k_casr
      double max_sr = 2.5;       // dimensionless

//     min_sr: Minimum value of k_casr
      double min_sr = 1.;        // dimensionless

//     V_leak: Maximal I_leak conductance
      double V_leak = 3.6E-4;    // units: mM/ms

//     V_xfer: Maximal I_xfer conductance
      double V_xfer = 3.8E-3;    // units: mM/ms

//     Buf_c: Total cytoplasmic buffer concentration
      double Buf_c = 0.2;        // units: mM

//     K_bufc: Ca_i half-saturation constant for cytplasmic buffer
      double K_bufc = 1.E-3;     // units: mM

//     Buf_sr: Total sacroplasmic buffer concentration
      double Buf_sr = 10.;       // units: mM

//     K_bufsr: Ca_sr half-saturation constant for subspace buffer
      double K_bufsr = 0.3;      // units: mM

//     Buf_ss: Total subspace buffer concentration
      double Buf_ss = 0.4;       // units: mM

//     K_bufss: Ca_ss half-saturation constant for subspace buffer
      double K_bufss = 2.5E-4;   // units: mM

//     Resting potential
      double Vrest = -85.23;     // units: mV

//     Electromechanics coupling parameters: active stress model
//     Ca_rest: Resting Ca concentration
      double Ca_rest = 5.E-5;    // units: mM

//     Ca_crit: Critical Ca concentration
      double Ca_crit = 8.E-4;    // units: mM

//     eta_T: Saturation of concentration
      double eta_T = 12.5;       // units: MPa/mM

//     eps_0: Minimum activation
      double eps_0 = 0.1;        // units: ms^{-1}

//     eps_i: Maximum activation
      double eps_i = 1. ;        // units: ms^{-1}

//     Transition rate
      double xi_T = 4.E3;        // units: mM^{-1}

//     Electromechanics coupling parameters: active strain model
//

//     Active force of sacromere (-mM^{-2})
      double alFa = -4.E6;

//     Resting Ca concentration (mM)
      double c_Ca0 = 2.155E-4;

//     Viscous-type constant (ms-mM^{-2})
      double mu_Ca = 5.E6;

//     Force-length relationship parameters
//     Initial length of sacromeres (um)
      double SL0 = 1.95;

//     Min. length of sacromeres (um)
      double SLmin = 1.7;

//     Max. length of sacromeres (um)
      double SLmax = 2.6;

//     Fourier coefficients
      double f0  = -4333.618335582119;
      double fc1 =  2570.395355352195;
      double fs1 = -2051.827278991976;
      double fc2 =  1329.53611689133;
      double fs2 =  302.216784558222;
      double fc3 =  104.943770305116;
      double fs3 =  218.375174229422;

//-----------------------------------------------------------------------
//     Scaling factors
//     Voltage scaling
      double Vscale  = 1.;

//     Time scaling
      double Tscale  = 1.;

//     Voltage offset parameter
      double Voffset = 0.;

//-----------------------------------------------------------------------
//     Variables
//     Reverse potentials for Na, K, Ca
      double E_Na;
      double E_K;
      double E_Ca;
      double E_Ks;
//     Cellular transmembrane currents
//     I_Na: Fast sodium current
      double I_Na;

//     I_K1: inward rectifier outward current
      double I_K1;

//     I_to: transient outward current
      double I_to;

//     I_Kr: rapid delayed rectifier current
      double I_Kr;

//     I_Ks: slow delayed rectifier current
      double I_Ks;

//     I_CaL: L-type Ca current
      double I_CaL;

//     I_NaCa: Na-Ca exchanger current
      double I_NaCa;

//     I_NaK: Na-K pump current
      double I_NaK;

//     I_pCa: plateau Ca current
      double I_pCa;

//     I_pK: plateau K current
      double I_pK;

//     I_bCa: background Ca current
      double I_bCa;

//     I_lean: background Na current
      double I_bNa;

//     I_leak: sacroplasmic reticulum Ca leak current
      double I_leak;

//     I_up: sacroplasmic reticulum Ca pump current
      double I_up;

//     I_rel: Ca induced Ca release current
      double I_rel;

//     I_xfer: diffusive Ca current
      double I_xfer;
//-----------------------------------------------------------------------
//     State variables
      double V;
      double K_i;
      double Na_i;
      double Ca_i;
      double Ca_ss;
      double Ca_sr;
      double R_bar;

//     Gating variables (runtime, steady state)
      double xr1, xr1i;
      double xr2, xr2i;
      double xs, xsi;
      double m, mi;
      double h, hi;
      double j, ji;
      double d, di;
      double f, fi;
      double f2, f2i;
      double fcass, fcassi;
      double s, si;
      double r, ri;

//     Other variables
      double k1;
      double k2;
      double k_casr;
      double O;

//     Jacobian variables
      double E_Na_Nai, E_K_Ki, E_Ca_Cai, E_Ks_Ki, E_Ks_Nai;
      double I_Na_V, I_Na_Nai;
      double I_to_V, I_to_Ki;
      double I_K1_V, I_K1_Ki;
      double I_Kr_V, I_Kr_Ki;
      double I_Ks_V, I_Ks_Ki, I_Ks_Nai;
      double I_CaL_V, I_CaL_Cass;
      double I_NaCa_V, I_NaCa_Nai, I_NaCa_Cai;
      double I_NaK_V, I_NaK_Nai;
      double I_pCa_Cai;
      double I_pK_V, I_pK_Ki;
      double I_bCa_V, I_bCa_Cai;
      double I_bNa_V, I_bNa_Nai;
      double I_leak_Cai, I_leak_Casr;
      double I_up_Cai;
      double I_rel_Cass, I_rel_Casr, I_rel_Rbar;
      double I_xfer_Cai, I_xfer_Cass;
      double k_casr_sr, k1_casr, O_Casr, O_Cass, O_Rbar;

    void actv_strn(const double c_Ca, const double I4f, const double dt, double& gf);
    void actv_strs(const double c_Ca, const double dt, double& Tact, double& epsX);

    void getf(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, 
        Vector<double>& dX, const double I_stim, const double K_sac, Vector<double>& RPAR);

    void getj(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, 
        Array<double>& JAC, const double Ksac);

    void init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg);

    void init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg,
        Vector<double>& X0, Vector<double>& Xg0);

    void integ_cn2(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg,
        const double Ts, const double dt, const double Istim, const double Ksac, 
        Vector<int>& IPAR, Vector<double>& RPAR);

    void integ_fe(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
        const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR);

    void integ_rk(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
        const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR);

    void update_g(const int i, const double dt, const int n, const int nG, const Vector<double>& X, 
        Vector<double>& Xg);

};

#endif

