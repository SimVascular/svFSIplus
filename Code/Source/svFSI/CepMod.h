
// The classes defined here duplicate the data structures in the Fortran CEPMOD module
// defined in CEPMOD.f. 

// This module defines data structures for cardiac electrophysiology
// model equation. It also interfaces with individual modules for
// the cellular activation model.


#ifndef CEP_MOD_H 
#define CEP_MOD_H 

#include "CepModAp.h"
#include "CepModBo.h"
#include "CepModFn.h"
#include "CepModTtp.h"
#include "consts.h"

#include "Array.h"
#include "Vector.h"
#include <map>

//----------------------------
// ElectrophysiologyModelType
//----------------------------
//
// Type of cardiac electrophysiology models.
//
enum class ElectrophysiologyModelType {
  NA = 100, 
  AP = 101,
  BO = 102, 
  FN = 103, 
  TTP = 104
};

extern const std::map<std::string,ElectrophysiologyModelType> cep_model_name_to_type;

// Print ElectrophysiologyModelType as a string.
static std::ostream &operator << ( std::ostream& strm, ElectrophysiologyModelType type)
{
  const std::map<ElectrophysiologyModelType, std::string> names = { 
    {ElectrophysiologyModelType::NA, "NA"}, 
    {ElectrophysiologyModelType::AP,"AP"}, 
    {ElectrophysiologyModelType::BO, "BO"}, 
    {ElectrophysiologyModelType::FN, "FN"}, 
    {ElectrophysiologyModelType::TTP, "TTP"}, 
  };
  return strm << names.at(type);
}

//--------------------
// TimeIntegratioType
//--------------------
// Time integration scheme.
//
enum class TimeIntegratioType {
  NA = 200, 
  FE = 201,
  RK4 = 202, 
  CN2 = 203
};

extern const std::map<std::string,TimeIntegratioType> cep_time_int_to_type;

static std::ostream &operator << ( std::ostream& strm, TimeIntegratioType type)
{
  const std::map<TimeIntegratioType, std::string> names = { 
    {TimeIntegratioType::NA, "NA"}, 
    {TimeIntegratioType::FE, "FE"}, 
    {TimeIntegratioType::RK4, "RK4"}, 
    {TimeIntegratioType::CN2, "CN2"}, 
  };
  return strm << names.at(type);
}

//--------------------
// TimeIntegratioType

//--------
// cmType
//--------
// Time integration scheme and related parameters
//
class odeType {
  public:
    odeType() {};

    // Time integration method type
    TimeIntegratioType tIntType = TimeIntegratioType::NA;
    //int tIntType = tIntType_NA;

    // Max. iterations for Newton-Raphson method
    int maxItr = 5;

    // Absolute tolerance
    double absTol = 1.E-8;

    // Relative tolerance
    double relTol = 1.E-4;
};

// External stimulus type
class stimType
{
  public:
    // start time
    double Ts = 0.0;

    // duration of stimulus
    double Td = 0.0;

    // cycle length
    double CL = 0.0;

    // stimulus amplitude
    double A = 0.0;
};

// ECG leads type
class ecgLeadsType
{
  public:
    // Number of leads
    int num_leads = 0;

    // x coordinates
    Vector<double> x_coords;

    // y coordinates
    Vector<double> y_coords;

    // z coordinates
    Vector<double> z_coords;

    // Pseudo ECG over each lead
    Vector<double> pseudo_ECG;

    // Output files
    std::vector<std::string> out_files;
};

//!     Cardiac electrophysiology model type

class cepModelType
{
  public:
    cepModelType();
    ~cepModelType();

    // Type of cardiac electrophysiology model
    ElectrophysiologyModelType cepType = ElectrophysiologyModelType::NA;

    // Number of state variables
    int nX = 0;

    // Number of gating variables
    int nG = 0;

    //  Number of fiber directions
    int nFn = 0;

    //  Myocardium zone id, default to epicardium.
    int imyo = 1;

    //  Time step for integration
    double dt = 0.0;

    //  Constant for stretch-activated-currents
    double Ksac = 0.0;

    //  Isotropic conductivity
    double Diso = 0.0;

    //  Anisotropic conductivity
    Vector<double> Dani;

    //  External stimulus
    stimType Istim;

    //  Time integration options
    odeType odes;
};

  //     Cardiac electromechanics model type
class cemModelType
{
  public:
    //  Whether electrophysiology and mechanics are coupled
    bool cpld = false;
    //bool cpld = .FALSE.

    //  Whether active stress formulation is employed
    bool aStress = false;
    //bool aStress = .FALSE.

    //  Whether active strain formulation is employed
    bool aStrain = false;
    //bool aStrain = .FALSE.

    //  Local variable integrated in time
    //    := activation force for active stress model
    //    := fiber stretch for active strain model
    Vector<double> Ya;
};

//--------
// CepMod
//--------
//
class CepMod 
{
  public:

    // Whether cardiac electrophysiology is solved
    bool cepEq;

    // Max. dof in cellular activation model
    int nXion = 0;

    // Unknowns stored at all nodes
    Array<double> Xion;

    // Cardiac electromechanics type
    cemModelType cem;

    // Interface for Aliev-Panfilov cellular activation model.
    CepModAp ap;

    // Interface for ABueno-Orovio cellular activation model.
    CepModBo bo;

    // Interface for Fitzhugh-Nagumo cellular activation model.
    CepModFn fn;

    // Interface for Tusscher-Panfilov cellular activation model.
    CepModTtp ttp;

    // ECG leads
    ecgLeadsType ecgleads;
};

#endif

