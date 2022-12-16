#ifndef MAT_MODELS_H 
#define MAT_MODELS_H 

#include "Array.h"
#include "CepMod.h"
#include "ComMod.h"
#include "Tensor4.h"

#include "mat_fun.h"

namespace mat_models {

void actv_strain(const ComMod& com_mod, const CepMod& cep_mod, const double gf, 
    const int nfd, const Array<double>& fl, Array<double>& Fa);

void cc_to_voigt(const int nsd, const Tensor4<double>& CC, Array<double>& Dm);

void get_fib_stress(const ComMod& com_mod, const CepMod& cep_mod, const fibStrsType& Tfl, double& g);

void get_pk2cc(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const Array<double>& F, const int nfd,
    const Array<double>& fl, const double ya, Array<double>& S, Array<double>& Dm);

void get_pk2cc_dev(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const Array<double>& F, const int nfd,
    const Array<double>& fl, const double ya, Array<double>& S, Array<double>& Dm, double& Ja);

void get_tau(const ComMod& com_mod, const dmnType& lDmn, const double detF, const double Je, double& tauM, double& tauC);

void get_svol_p(const ComMod& com_mod, const CepMod& cep_mod, const stModelType& stM, const double J, 
    double& p, double& pl);

void g_vol_pen(const ComMod& com_mod, const dmnType& lDmn, const double p, 
    double& ro, double& bt, double& dro, double& dbt, const double Ja);

};

//#include "mat_models_fixed.h"

#endif

