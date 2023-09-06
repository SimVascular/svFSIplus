#ifndef CONSTS_H 
#define CONSTS_H 

#include <iostream>
#include <limits>
#include <map>
#include <type_traits>

// The enums here replicate the PARAMETERs defined
// in CONSTS.f.

namespace consts {

const double pi = 3.1415926535897932384626;

const int maxNSD = 3;

const int maxNProp = 20;

const int maxOutput = 5;

// Use inf numeric values to represent a value that is not set.
const int int_inf = std::numeric_limits<int>::infinity();
const double double_inf = std::numeric_limits<double>::infinity();

//----------
// enum_int
//-----------
//
template<typename T>
int enum_int(T value)
{
  return static_cast<int>(value);
}

//---------
// present
//---------
// Check if a value is set to infinity.
//
template<typename T>
bool present(T value)
{
  return (value != std::numeric_limits<T>::infinity()); 
}

//---------------
// BodyForceType
//---------------
// Body force types: volumetric (default), traction, Neumann
// (pressure based), time dependence (steady, unsteady, spatially
// varying, general)
//
enum class BodyForceType 
{
  bfType_vol = 0, 
  bfType_trac = 1,
  bfType_Neu = 2, 
  bfType_std = 3, 
  bfType_ustd = 4,
  bfType_spl = 5, 
  bfType_gen = 6
};

//-----------------------
// BoundaryConditionType
//-----------------------
//Boundary conditions type.
//
// BC types are stored as bitwise values.
//
// boundary conditions types. Items of this list can be combined
//
//   BCs from imposing perspective can be Neu/Dir/per
//
//   BCs time dependence can be std/ustd/cpl/gen/res
//
//   BCs spatial distribution can be flat/para/ud
//
//   Beside these nodes at the boundary perimeter can be set to
//   zero and flux through surface can be assigned instead of nodal
//   values.
//
//   Dirichlet, Neumann, Traction, CMM, Robin, steady, unsteady,
//   coupled, general (combination of ud/ustd), resistance, imposed
//   flux, zero out perimeter, impose BC on the integral of state
//   variable or D (instead of Y), flat profile, parabolic profile,
//   user defined profile, backflow stabilization, BCs for shells
//   (fixed, hinged, free, symmetric), undeforming Neu, RCR-Neu

enum class BoundaryConditionType 
{
  bType_Dir = 0,       // Dirichlet
  bType_Neu = 1,       // Neumann
  bType_trac = 2,      // Traction
  bType_CMM = 3,       // CMM
  bType_Robin = 4,     // RObin
  bType_std = 5,       // steady
  bType_ustd = 6,      // unsteady
  bType_cpl = 7,       // coupled
  bType_gen = 8,       // general
  bType_res = 9,       // resistance
  bType_flx = 10,      // imposed flux 
  bType_zp = 11,       // zero out perimeter
  bType_impD = 12,     // impose BC on the integral of state variable or D (instead of Y)
  bType_flat =13,      // flat profile
  bType_para = 14,     // parabolic profile
  bType_ud = 15,       // user defined profile
  bType_bfs = 16,      // backflow stabilization
  bType_fix = 17,      // shell fixed
  bType_hing = 18,     // shell hinged
  bType_free = 19,     // shell free
  bType_symm = 20,     // shell symmetric 
  bType_undefNeu = 21, // undeforming Neu
  bType_RCR = 22       // RCR-Neu
};

// Define constants using smaller name and integer value (needed for bitwise operations).
//
constexpr auto BC_CMM = BoundaryConditionType::bType_CMM;
constexpr auto iBC_CMM = static_cast<int>(BoundaryConditionType::bType_CMM);

constexpr auto BC_cpl = BoundaryConditionType::bType_cpl;
constexpr auto iBC_cpl = static_cast<int>(BoundaryConditionType::bType_cpl);

constexpr auto BC_Dir = BoundaryConditionType::bType_Dir;
constexpr auto iBC_Dir = static_cast<int>(BoundaryConditionType::bType_Dir);

constexpr auto BC_fix = BoundaryConditionType::bType_fix;
constexpr auto iBC_fix = static_cast<int>(BoundaryConditionType::bType_fix);

constexpr auto BC_flat = BoundaryConditionType::bType_flat;
constexpr auto iBC_flat = static_cast<int>(BoundaryConditionType::bType_flat);

constexpr auto BC_free = BoundaryConditionType::bType_free;
constexpr auto iBC_free = static_cast<int>(BoundaryConditionType::bType_free);

constexpr auto BC_gen = BoundaryConditionType::bType_gen;
constexpr auto iBC_gen = static_cast<int>(BoundaryConditionType::bType_gen);

constexpr auto BC_hing = BoundaryConditionType::bType_hing;
constexpr auto iBC_hing = static_cast<int>(BoundaryConditionType::bType_hing);

constexpr auto BC_impD = BoundaryConditionType::bType_impD;
constexpr auto iBC_impD = static_cast<int>(BoundaryConditionType::bType_impD);

constexpr auto BC_Neu = BoundaryConditionType::bType_Neu;
constexpr auto iBC_Neu = static_cast<int>(BoundaryConditionType::bType_Neu);

constexpr auto BC_para = BoundaryConditionType::bType_para;
constexpr auto iBC_para = static_cast<int>(BoundaryConditionType::bType_para);

constexpr auto BC_RCR = BoundaryConditionType::bType_RCR;
constexpr auto iBC_RCR = static_cast<int>(BoundaryConditionType::bType_RCR);

constexpr auto BC_res = BoundaryConditionType::bType_res;
constexpr auto iBC_res = static_cast<int>(BoundaryConditionType::bType_res);

constexpr auto BC_Robin = BoundaryConditionType::bType_Robin;
constexpr auto iBC_Robin = static_cast<int>(BoundaryConditionType::bType_Robin);

constexpr auto BC_std = BoundaryConditionType::bType_std;
constexpr auto iBC_std = static_cast<int>(BoundaryConditionType::bType_std);

constexpr auto BC_symm = BoundaryConditionType::bType_symm;
constexpr auto iBC_symm = static_cast<int>(BoundaryConditionType::bType_symm);

constexpr auto BC_trac = BoundaryConditionType::bType_trac;
constexpr auto iBC_trac = static_cast<int>(BoundaryConditionType::bType_trac);

constexpr auto BC_undefNeu = BoundaryConditionType::bType_undefNeu;
constexpr auto iBC_undefNeu = static_cast<int>(BoundaryConditionType::bType_undefNeu);

constexpr auto BC_ustd = BoundaryConditionType::bType_ustd;
constexpr auto iBC_ustd = static_cast<int>(BoundaryConditionType::bType_ustd);

//-----------------------
// ConstitutiveModelType
//-----------------------
// Constitutive model (isochoric) type for structure equation:
//
enum class ConstitutiveModelType 
{
  stIso_NA = 600,
  stIso_StVK = 601, 
  stIso_mStVK = 602, 
  stIso_nHook = 603,
  stIso_MR = 604, 
  stIso_HGO = 605, 
  stIso_lin = 606,
  stIso_Gucci = 607, 
  stIso_HO = 608,
  stIso_HO_ma = 610,
  stIso_LS = 611,
  stVol_NA = 650,
  stVol_Quad = 651, 
  stVol_ST91 = 652, 
  stVol_M94 = 653
};

// Map for constitutive_model string to ConstitutiveModelType. 
extern const std::map<std::string,ConstitutiveModelType> constitutive_model_name_to_type;

//------------------
// ContactModelType
//------------------
//
enum class ContactModelType
{
  cntctM_NA = 800,
  cntctM_penalty = 801,
  cntctM_potential = 802
};

// Map for model type string to ContactModelType. 
extern const std::map<std::string,ContactModelType> contact_model_name_to_type;

// Differenty type of coupling for cplBC. 
//
// INTEGER(KIND=IKIND), PARAMETER :: cplBC_NA = 400, cplBC_I = 401,
// cplBC_SI = 402, cplBC_E = 403
//
// INTEGER(KIND=IKIND), PARAMETER :: cplBC_Dir = 66112, cplBC_Neu = 66113
//
enum class CplBCType 
{
  cplBC_NA = 400,
  cplBC_Dir = 66112,   // Dirichlet type coupling
  cplBC_E = 403,       // explicit
  cplBC_I = 401,       // implicit
  cplBC_Neu = 66113,   // Neumann type coupling
  cplBC_SI = 402,      // semi-implicit
};

// Map for cplBC type to CplBCType. 
extern const std::map<std::string,CplBCType> cplbc_name_to_type;

//-------------
// ElementType
//-------------
// Element type replicating eType_NA, eType_PNT, etc. 
//
enum class ElementType 
{
  NA = 100, 
  PNT = 101,
  LIN1 = 102, 
  LIN2 = 103, 
  TRI3 = 104,
  TRI6 = 105, 
  QUD4 = 106, 
  QUD8 = 107,
  QUD9 = 108, 
  TET4 = 109, 
  TET10 = 110,
  HEX8 = 111, 
  HEX20 = 112, 
  HEX27 = 113,
  WDG = 114, 
  NRB = 115
};

extern const std::map<ElementType,int> element_type_to_elem_nonb;

// Template for printing ElementType.
/*
template<typename T>
std::ostream& operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e)
{
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}
*/

//--------------
// EquationType
//--------------
// Types of equations that are included in this solver.
//
//  Fluid equation (Navier-Stokes), nonlinear structure (pure d), heat
//  equation, linear elasticity, heat in fluid (advection-diffusion),
//  fluid-structure-interaction, mesh motion, Shell mechanics,
//  Coupled-Momentum-Method, Cardiac Electro-Physiology,
//  Nonlinear structure (v-p), Stokes equations
//
enum class EquationType 
{
  phys_NA = 200, 
  phys_fluid = 201,
  phys_struct = 202,  // nonlinear structure (pure d)
  phys_heatS = 203, 
  phys_lElas = 204,
  phys_heatF = 205, 
  phys_FSI = 206, 
  phys_mesh = 207,    // solves a modified lElas for mesh motion; should be used with FSI
  phys_shell = 208,   // solves nonlinear thin shell mechanics (Kirchhoff-Love theory)
  phys_CMM = 209, 
  phys_CEP = 210,
  phys_ustruct = 211,  // Nonlinear elastodynamics using mixed VMS-stabilized formulation 
  phys_stokes = 212
};

constexpr auto Equation_CMM = EquationType::phys_CMM;
constexpr auto Equation_CEP = EquationType::phys_CEP;
constexpr auto Equation_fluid = EquationType::phys_fluid;
constexpr auto Equation_FSI = EquationType::phys_FSI;
constexpr auto Equation_heatF = EquationType::phys_heatF;
constexpr auto Equation_heatS = EquationType::phys_heatS;
constexpr auto Equation_lElas = EquationType::phys_lElas;
constexpr auto Equation_mesh = EquationType::phys_mesh;
constexpr auto Equation_shell = EquationType::phys_shell;
constexpr auto Equation_stokes = EquationType::phys_stokes;
constexpr auto Equation_struct = EquationType::phys_struct;
constexpr auto Equation_ustruct = EquationType::phys_ustruct;

extern const std::map<std::string,EquationType> equation_name_to_type;

//-------------------
// MeshGeneratorType
//-------------------
//
enum class MeshGeneratorType
{
  RMSH_TETGEN = 1,
  RMSH_MESHSIM = 2
};

// Map for string to MeshGeneratorType. 
extern const std::map<std::string,MeshGeneratorType> mesh_generator_name_to_type;

//------------
// OutputType
//------------
//
enum class OutputType 
{
  outGrp_NA = 500, 
  outGrp_A = 501,
  outGrp_Y = 502, 
  outGrp_D = 503, 
  outGrp_I = 504, 
  outGrp_WSS = 505, 
  outGrp_trac = 506, 
  outGrp_vort = 507, 
  outGrp_vortex = 508,
  outGrp_stInv = 509, 
  outGrp_eFlx = 510, 
  outGrp_hFlx = 511,
  outGrp_absV = 512, 
  outGrp_fN = 513, 
  outGrp_fA = 514,
  outGrp_stress = 515, 
  outGrp_cauchy = 516, 
  outGrp_mises = 517,
  outGrp_J = 518, 
  outGrp_F = 519, 
  outGrp_strain = 520,
  outGrp_divV = 521, 
  outGrp_Visc = 522,
  outGrp_fS = 523,
  outGrp_C = 524, 
  outGrp_I1 = 525,

  out_velocity = 599,
  out_pressure = 598, 
  out_temperature = 597, 
  out_voltage = 596,
  out_acceleration = 595, 
  out_displacement = 594, 
  out_integ =593,
  out_WSS = 592, 
  out_traction = 591, 
  out_vorticity = 590,
  out_vortex = 589, 
  out_strainInv = 588, 
  out_energyFlux = 587,
  out_heatFlux = 586, 
  out_absVelocity = 585, 
  out_fibDir = 584,
  out_fibAlign = 583, 
  out_stress = 582, 
  out_cauchy = 581,
  out_mises = 580, 
  out_jacobian = 579, 
  out_defGrad = 578,
  out_strain = 577, 
  out_divergence = 576, 
  out_viscosity = 575,
  out_fibStrn = 574,
  out_CGstrain = 573,
  out_CGInv1 = 572
};

//---------------------
// PhysicalProperyType
//---------------------
// Possible physical properties. Current maxNPror is 20.
//
enum class PhysicalProperyType 
{
  NA = 0, 
  fluid_density = 1, 
  solid_density = 2, 
  solid_viscosity = 3, 
  elasticity_modulus = 4,
  poisson_ratio = 5, 
  conductivity = 6, 
  f_x = 7,                     // internal force x
  f_y = 8,                     // internal force y
  f_z = 9,                     // internal force z
  backflow_stab = 10,          // stabilization coeff. for backflow divergence
  source_term = 11,            // external source
  damping = 12,
  shell_thickness = 13, 
  ctau_M = 14,                 // stabilization coeffs. for USTRUCT (momentum, continuity)
  ctau_C = 15
};

//--------------------
// PreconditionerType
//--------------------
//
enum class PreconditionerType 
{
  PREC_NONE = 700,
  PREC_FSILS = 701, 
  PREC_TRILINOS_DIAGONAL = 702,
  PREC_TRILINOS_BLOCK_JACOBI = 703, 
  PREC_TRILINOS_ILU = 704,
  PREC_TRILINOS_ILUT = 705, 
  PREC_TRILINOS_IC = 706,
  PREC_TRILINOS_ICT = 707, 
  PREC_TRILINOS_ML = 708,
  PREC_RCS = 709
};

// Map for preconditioner type string to pair (PreconditionerType enum, bool(true if Trilinos precondition)). 
//  
using PreconditionerMapType = std::pair<PreconditionerType,bool>;
//extern const std::map<std::string,PreconditionerType> preconditioner_name_to_type;
extern const std::map<std::string,PreconditionerMapType> preconditioner_name_to_type;

//------------
// SolverType
//------------
//
enum class SolverType
{
  lSolver_NA = 799,
  lSolver_CG = 798, 
  lSolver_GMRES = 797, 
  lSolver_NS = 796,
  lSolver_BICGS = 795
};

// Map for solver type string to SolverType enum. 
extern const std::map<std::string,SolverType> solver_name_to_type;

//-------------------------
// FluidViscosityModelType
//-------------------------
//
enum class FluidViscosityModelType 
{
  viscType_CY = 697, 
  viscType_Cass = 696,
  viscType_Const = 698, 
  viscType_NA = 699
};

// Map for fluid viscosity model string to FluidViscosityModelType. 
extern const std::map<std::string,FluidViscosityModelType> fluid_viscosity_model_name_to_type;

// Template for printing enum class types as an int.
template<typename T>
std::ostream& operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e)
{
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}

};

#endif

