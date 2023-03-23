
// The classes defined here duplicate the data structures in the Fortran COMMOD module
// defined in MOD.f. 
//
// All of the data structures used for the mesh, boundarsy conditions and solver parameters, etc. 
// are defined here.

#ifndef COMMOD_H 
#define COMMOD_H 

#include "Array.h"
#include "Array3.h"
#include "CepMod.h"
#include "ChnlMod.h"
#include "CmMod.h"
#include "Timer.h"
#include "Vector.h"

#include "DebugMsg.h"

#include "consts.h"

#include "fils_struct.hpp"

#include <array>
#include <iostream>
#include <string>
#include <vector>

//--------
// fcType
//--------
// Fourier coefficients that are used to specify unsteady BCs
//
class fcType
{
  public:

    bool defined() { return n != 0; };

    // If this is a ramp function
    bool lrmp = false;

    // Number of Fourier coefficient
    int n = 0;

    // No. of dimensions (scalar or vector)
    int d = 0;
   
    // Initial value
    Vector<double> qi;

    // Time derivative of linear part
    Vector<double> qs;

    // Period
    double T = 0.0;

    // Initial time
    double ti = 0.0;

    // Imaginary part of coefficint
    Array<double> i;

    // Real part of coefficint
    Array<double> r;
};

//--------
// MBType
//--------
// Moving boundary data structure (used for general BC)
//
class MBType
{
  public:

    bool defined() { return dof != 0; };

    // Degrees of freedom of d(:,.,.)
    int dof = 0;

    // Number of time points to be read
    int nTP = 0;

    // The period of data
    double period = 0.0;

    // Time points
    Vector<double> t;

    // Displacements at each direction, location, and time point
    Array3<double> d;
};

//---------
// rcrType
//---------
//
class rcrType
{
  public:

    // Proximal resistance
    double Rp = 0.0;

    // Capacitance
    double C = 0.0;

    // Distance resistance
    double Rd = 0.0;

    // Distal pressure
    double Pd = 0.0;

    // Initial value
    double Xo = 0.0;
};



//--------
// bcType
//--------
// Boundary condition data type
//
class bcType
{
  public:

    // Strong/Weak application of Dirichlet BC
    bool weakDir = false;

    // Whether load vector changes with deformation
    // (Neu - struct/ustruct only)
    bool flwP = false;

    // Robin: apply only in normal direction
    bool rbnN = false;

    // Pre/Res/Flat/Para... boundary types
    //
    // This stores differnt BCs as bitwise values. 
    //
    int bType = 0;

    // Pointer to coupledBC%face
    int cplBCptr = -1;

    // The face index that corresponds to this BC
    int iFa = -1;

    // The mesh index that corresponds to this BC
    int iM = -1;

    // Pointer to FSILS%bc
    int lsPtr = -1;

    // Undeforming Neu BC master-slave node parameters.
    int masN = 0;

    // Defined steady value
    double g = 0.0;

    // Neu: defined resistance
    double r = 0.0;

    // Robin: stiffness
    double k = 0.0;

    // Robin: damping
    double c = 0.0;

    // Penalty parameters for weakly applied Dir BC
    Vector<double> tauB{0.0, 0.0};
    //double tauB[2];

    // Direction vector for imposing the BC
    Vector<int> eDrn;

    // Defined steady vector (traction)
    Vector<double> h;

    // Spatial dependant BC (profile data)
    Vector<double> gx;

    // General BC (unsteady and UD combination)
    //
    // This is declare ALLOCATABLE in MOD.f. 
    //
    MBType gm;

    // Time dependant BC (Unsteady imposed value);
    //
    // This is declare ALLOCATABLE in MOD.f. 
    //
    fcType gt;

    // Neu: RCR
    rcrType RCR;
};

//--------
// bsType
//--------
// Class storing data for B-Splines.
//
class bsType 
{
  public:

    // Number of knots (p + nNo + 1)
    int n = 0;

    // Number of Gauss points for integration
    int nG = 0;

    // Number of knot spans (element)
    int nEl = 0;

    // Number of control points (nodes)
    int nNo = 0;

    // Number of sample points in each element (for output)
    int nSl = 0;

    // The order
    int p = 0;

    // Knot vector.
    Vector<double> xi;
};

//--------
// fsType
//--------
// Function spaces (basis) type.
//
class fsType {

  public:

    fsType();

    void destroy();

    // Whether the basis function is linear
    bool lShpF = false;

    // Element type
    consts::ElementType eType = consts::ElementType::NA;

    // Number of basis functions, typically equals msh%eNoN
    int eNoN = 0;

    // Number of Gauss points for integration
    int nG = 0;

    // Gauss weights
    Vector<double> w;

    // Gauss integration points in parametric space
    Array<double> xi;

    // Bounds on Gauss integration points in parametric space
    Array<double> xib;

    // Parent shape function
    Array<double> N;

    // Bounds on shape functions
    Array<double> Nb;

    // Parent shape functions gradient
    Array3<double> Nx;

    // Second derivatives of shape functions - used for shells & IGA
    Array3<double> Nxx;
};

//--------
// bfType
//--------
// Body force data structure type
//
class bfType
{
  public:

    std::string file_name;
    std::string mesh_name;

    // Type of body force applied
    int bType = 0;

    // No. of dimensions (1 or nsd)
    int dof = 0;

    // Mesh index corresponding to this body force
    int iM = -1;

    // Steady value
    Vector<double> b;

    // Steady but spatially dependant
    Array<double> bx;

    // Time dependant (unsteady imposed value)
    fcType bt;

    // General (unsteady and spatially dependent combination)
    MBType bm;
};

// Fiber stress type
class fibStrsType
{
  public:

    // Type of fiber stress
    int fType = 0;

    // Constant steady value
    double g = 0.0;

    // Unsteady time-dependent values
    fcType gt;
};

//-------------
// stModelType
//-------------
// Structural domain type
//
class stModelType
{
  public:
    // Type of constitutive model (volumetric) for struct/FSI
    consts::ConstitutiveModelType volType = consts::ConstitutiveModelType::stIso_NA;

    // Penalty parameter
    double Kpen = 0.0;

    // Type of constitutive model (isochoric) for struct/FSI
    consts::ConstitutiveModelType isoType = consts::ConstitutiveModelType::stIso_NA;

    // Parameters specific to the constitutive model (isochoric)
    // NeoHookean model (C10 = mu/2)
    double C10 = 0.0;

    // Mooney-Rivlin model (C10, C01)
    double C01 = 0.0;

    // Holzapfel model(a, b, aff, bff, ass, bss, afs, bfs, kap)
    double a = 0.0;
    double b = 0.0;
    double aff = 0.0;
    double bff = 0.0;
    double ass = 0.0;
    double bss = 0.0;
    double afs = 0.0;
    double bfs = 0.0;

    // Collagen fiber dispersion parameter (Holzapfel model)
    double kap = 0.0;

    // Fiber reinforcement stress
    fibStrsType Tf;
};

//---------------
// viscModelType
//---------------
// Fluid viscosity model type
//
class viscModelType
{
  public:

    // Type of constitutive model for fluid viscosity
    consts::FluidViscosityModelType viscType = consts::FluidViscosityModelType::viscType_NA;

    // Limiting zero shear-rate viscosity value
    double mu_o = 0.0;

    // Limiting high shear-rate viscosity (asymptotic) value
    double mu_i = 0.0;

    // Strain-rate tensor multiplier
    double lam = 0.0;

    // Strain-rate tensor exponent
    double a = 0.0;

    // Power-law exponent
    double n = 0.0;
};

//     Domain type is to keep track with element belong to which domain
//     and also different physical quantities
class dmnType
{
  public:
    dmnType();
    ~dmnType();

    // The domain ID. Default includes entire domain
    int Id = -1;

    // Which physics must be solved in this domain
    consts::EquationType phys = consts::EquationType::phys_NA;

    // The volume of this domain
    double v = 0.0;

    // General physical properties such as density, elastic modulus...
    // FIX davep double prop[maxNProp] ;
    std::map<consts::PhysicalProperyType,double> prop;
    //double prop[consts::maxNProp];

    // Electrophysiology model
    cepModelType cep;

    // Structure material model
    stModelType stM;

    // Viscosity model for fluids
    viscModelType visc;
};

//---------
// adjType
//---------
// Mesh adjacency (neighboring element for each element)
//
class adjType
{
  public:
    void destroy();

    // No of non-zeros
    int nnz = 0;

    // Column pointer
    Vector<int> pcol;

    // Row pointer
    Vector<int> prow;

};

// Tracer type used for immersed boundaries. Identifies traces of
// nodes and integration points on background mesh elements
class traceType
{
  public:
    // No. of non-zero nodal traces
    int n = 0;

    // No. of non-zero integration point traces
    int nG = 0;

    // Local to global nodes maping nNo --> tnNo
    Vector<int> gN;

    // Self pointer of each trace to the IB integration point and
    // element ID
    Array<int> gE;

    // Nodal trace pointer array stores two values for each trace.
    // (1) background mesh element to which the trace points to,
    // (2) mesh ID
    Array<int> nptr;

    // Integration point tracer array stores two values for each trace
    // (1) background mesh element to which the trace points to,
    // (2) mesh ID
    Array<int> gptr;

    // Parametric coordinate for each nodal trace
    Array<double> xi;

    // Parametric coordinate for each Gauss point trace
    Array<double> xiG;
};

//----------
// faceType
//----------
// The face type containing mesh at boundary
//
class faceType
{
  public:
    faceType();
    ~faceType();

    void destroy();

    //faceType& operator=(const faceType& rhs);

    // Parametric direction normal to this face (NURBS)
    int d = 0;

    // Number of nodes (control points) in a single element
    int eNoN = 0;

    // Element type
    consts::ElementType eType = consts::ElementType::NA;

    // The mesh index that this face belongs to
    int iM = 0;

    // Number of elements
    int nEl = 0;

    // Global number of elements
    int gnEl = 0;

    // Number of function spaces
    int nFs = 0;

    // Number of Gauss points for integration
    int nG = 0;

    // Number of nodes
    int nNo = 0;

    // Global element Ids
    Vector<int> gE;

    // Global node Ids
    Vector<int> gN;

    // Global to local maping tnNo --> nNo
    Vector<int> lN;

    // Connectivity array
    Array<int> IEN;

    // EBC array (gE + gIEN)
    Array<int> gebc;

    // Surface area
    double area = 0.0;

    // Gauss point weights
    Vector<double> w;

    // Position coordinates
    Array<double> x;

    // Gauss points in parametric space
    Array<double> xi;

    // Shape functions at Gauss points
    Array<double> N;

    // Normal vector to each nodal point
    Array<double> nV;

    // Shape functions derivative at Gauss points
    // double Nx(:,:,:);
    Array3<double> Nx;

    // Second derivatives of shape functions - for shells & IGA
    // double Nxx(:,:,:);
    Array3<double> Nxx;

    // Face name for flux files
    std::string name;

    // Face nodal adjacency
    adjType nAdj;

    // Face element adjacency
    adjType eAdj;

    // Function spaces (basis)
    std::vector<fsType> fs;
};

// Declared type for outputed variables
class outputType
{
  public:

    // Is this output suppose to be written into VTK, boundary, vol
    std::vector<bool> wtn{false, false, false};

    // The group that this belong to (one of outType_*)
    consts::OutputType grp = consts::OutputType::outGrp_NA;
    //int grp;

    // Length of the outputed variable
    int l = 0;

    // Offset from the first index
    int o = 0;

    // The name to be used for the output and also in input file
    std::string name;
};


// Linear system of equations solver type
class lsType
{
  public:

    // LS solver                     (IN)
    consts::SolverType LS_type = consts::SolverType::lSolver_NA;

    // Preconditioner                (IN)
    consts::PreconditionerType PREC_Type = consts::PreconditionerType::PREC_NONE;

    // Successful solving            (OUT)
    bool suc = false;

    // Maximum iterations            (IN)
    int mItr = 1000;

    // Space dimension               (IN)
    int sD = 0;

    // Number of iteration           (OUT)
    int itr = 0;

    // Number of Ax multiple         (OUT)
    int cM = 0;

    // Number of |x| norms           (OUT)
    int cN = 0;

    // Number of <x.y> dot products  (OUT)
    int cD = 0;

    // Only for data alignment       (-)
    int reserve = 0;

    // Absolute tolerance            (IN)
    double absTol = 1e-08;

    // Relative tolerance            (IN)
    double relTol = 0.0;

    // Initial norm of residual      (OUT)
    double iNorm = 0.0;

    // Final norm of residual        (OUT)
    double fNorm = 0.0;

    // Res. rduction in last itr.    (OUT)
    double dB = 0.0;

    // Calling duration              (OUT)
    double callD = 0.0;
};


// Contact model type
class cntctModelType
{
  public:
    // Contact model
    int cType = 0;

    // Penalty parameter
    double k = 0.0;

    // Min depth of penetration
    double h = 0.0;

    // Max depth of penetration
    double c = 0.0;

    // Min norm of face normals in contact
    double al = 0.0;

    // Tolerance
    double tol = 0.0; 
};

class cplFaceType
{
  public:
    // GenBC_Dir/GenBC_Neu
    consts::CplBCType bGrp = consts::CplBCType::cplBC_NA;

    // Pointer to X
    int Xptr = -1;

    // Internal genBC use
    int eqv = 0;

    // Flow rates at t
    double Qo = 0.0;

    // Flow rates at t+dt
    double Qn = 0.0;

    // Pressures at t
    double Po = 0.0;

    // Pressures at t+dt
    double Pn = 0.0;

    // Imposed flow/pressure
    double y = 0.0;

    // Name of the face
    std::string name;

    // RCR type BC
    rcrType RCR;
};


//-----------
// cplBCType
//-----------
// For coupled 0D-3D problems
//
class cplBCType
{
  public:
    cplBCType();
    // Is multi-domain active
    bool coupled = false;

    //  Whether to use genBC
    bool useGenBC = false;

    //  Whether to initialize RCR from flow data
    bool initRCR = false;

    //  Number of coupled faces
    int nFa = 0;

    //  Number of unknowns in the 0D domain
    int nX = 0;

    //  Number of output variables addition to nX
    int nXp = 0;

    //  Implicit/Explicit/Semi-implicit schemes
    consts::CplBCType schm = consts::CplBCType::cplBC_NA;
    //int schm = cplBC_NA;

    //  Path to the 0D code binary file
    std::string binPath;

    //  File name for communication between 0D and 3D
    std::string commuName;
    //std::string commuName = ".CPLBC_0D_3D.tmp";

    //  The name of history file containing "X"
    std::string saveName;
    //std::string(LEN=stdL) :: saveName = "LPN.dat";

    //  New time step unknowns in the 0D domain
    Vector<double> xn;

    //  Old time step unknowns in the 0D domain
    Vector<double> xo;

    //  Output variables to be printed
    Vector<double> xp;

    //  Data structure used for communicating with 0D code
    std::vector<cplFaceType> fa;
};

//---------
// mshType
//---------
// This is the container for a mesh or NURBS patch, those specific
// to NURBS are noted
//
class mshType
{
  public:
    mshType();
    std::string dname = "";

/*
    mshType(const mshType &other) 
    { 
      std::cout << "c c c c c mshType copy c c c c c" << std::endl;
    } 
      
    mshType& operator = (const mshType &other)
    {
      std::cout << "= = = = = mshType assignment = = = = =" << std::endl;
      return *this;
    } 
*/

    ~mshType() 
    { 
      std::cout << "- - - - -  mshType dtor - - - - -   dname: " << dname << std::endl;
    };

    // Whether the shape function is linear
    bool lShpF = false;

    // Whether the mesh is shell
    bool lShl = false;

    // Whether the mesh is fibers (Purkinje)
    bool lFib = false;

    // Element type
    consts::ElementType eType = consts::ElementType::NA;
    //int eType = eType_NA

    // Number of nodes (control points) in a single element
    int eNoN = 0;

    // Global number of elements (knot spans)
    int gnEl = 0;

    // Global number of nodes (control points)
    int gnNo = 0;

    // Number of element face. Used for reading Gambit mesh files
    int nEf = 0;

    // Number of elements (knot spans)
    int nEl = 0;

    // Number of faces
    int nFa = 0;

    // Number of function spaces
    int nFs = 0;

    // Number of Gauss points for integration
    int nG = 0;

    // Number of nodes (control points) for 2D elements?
    int nNo = 0;

    // Number of elements sample points to be outputs (NURBS)
    int nSl = 0;

    // The element type recognized by VTK format
    int vtkType = 0;

    // Number of fiber directions
    int nFn = 0;

    // Mesh scale factor
    double scF = 0.0;

    // IB: Mesh size parameter
    double dx = 0.0;

    // Element distribution between processors
    Vector<int> eDist;

    // Element domain ID number
    Vector<int> eId;

    // Global nodes maping nNo --> tnNo
    Vector<int> gN;

    // GLobal projected nodes mapping projected -> unprojected mapping
    Vector<int> gpN;

    // Global connectivity array mappig eNoN,nEl --> gnNo
    Array<int> gIEN;

    // The connectivity array mapping eNoN,nEl --> nNo
    Array<int> IEN;

    // gIEN mapper from old to new
    Vector<int> otnIEN;

    // Local knot pointer (NURBS)
    Array<int> INN;

    // Global to local maping tnNo --> nNo
    Vector<int> lN;

    // Shells: extended IEN array with neighboring nodes
    Array<int> eIEN;

    // Shells: boundary condition variable
    Array<int> sbc;

    // IB: Whether a cell is a ghost cell or not
    Vector<int> iGC;

    // Control points weights (NURBS)
    Vector<double> nW;

    // Gauss weights
    Vector<double> w;

    // Gauss integration points in parametric space
    Array<double> xi;

    // Bounds on parameteric coordinates
    Array<double> xib;

    // Position coordinates
    Array<double> x;

    // Parent shape function
    Array<double> N;

    // Shape function bounds
    Array<double> Nb;

    // Normal vector to each nodal point (for Shells)
    Array<double> nV;

    // Fiber orientations stored at the element level - used for
    // electrophysiology and solid mechanics
    Array<double> fN;

    // Parent shape functions gradient
    // double Nx(:,:,:)
    Array3<double> Nx;

    // Second derivatives of shape functions - used for shells & IGA
    // davep double Nxx(:,:,:)
    Array3<double> Nxx;

    // Mesh Name
    std::string name;

    // Mesh nodal adjacency
    adjType nAdj;

    // Mesh element adjacency
    adjType eAdj;

    // Function spaces (basis)
    std::vector<fsType> fs;

    // BSpline in different directions (NURBS)
    std::vector<bsType> bs;

    // Faces are stored in this variable
    std::vector<faceType> fa;

    // IB: tracers
    traceType trc;

  private:
    //mshType(const mshType&);
    //mshType& operator=(const mshType&);

};

//--------
// eqType
//--------
// Equation type
//
class eqType
{
  public:
    eqType();
    ~eqType();

    // Should be satisfied in a coupled/uncoupled fashion
    bool coupled = false;
    //bool coupled = .TRUE.

    // Satisfied/not satisfied
    bool ok = false;

    //  Use C++ Trilinos framework for the linear solvers
    bool useTLS = false;

    //  Use C++ Trilinos framework for assembly and for linear solvers
    bool assmTLS = false;

    //  Degrees of freedom
    int dof = 0;

    //  Pointer to end of unknown Yo(:,s:e)
    int e = -1;

    //  Pointer to start of unknown Yo(:,s:e)
    int s = -1;

    //  Number of performed iterations
    int itr = 0;

    //  Maximum iteration for this eq.
    int maxItr = 5;

    //  Minimum iteration for this eq.
    int minItr = 1;

    //  Number of possible outputs
    int nOutput = 0;

    //  IB: Number of possible outputs
    int nOutIB = 0;

    //  Number of domains
    int nDmn = 0;

    //  IB: Number of immersed domains
    int nDmnIB = 0;

    //  Number of BCs
    int nBc = 0;

    //  IB: Number of BCs on immersed surfaces
    int nBcIB = 0;

    //  Number of BFs
    int nBf = 0;

    //  Type of equation fluid/heatF/heatS/lElas/FSI
    consts::EquationType phys = consts::EquationType::phys_NA;

    // Parameters used for the Generalized α− Method.
    //
    //  \alpha_f
    double af = 0.0;

    //  \alpha_m
    //
    // For second order equations: am = (2.0 - roInf) / (1.0 + roInf)
    // First order equations: am = 0.5 * (3.0 - roInf) / (1.0 + roInf)
    //
    double am = 0.0;

    //  \beta
    double beta = 0.0;

    //  \gamma
    double gam = 0.0;


    //  Initial norm of residual
    double iNorm = 0.0;

    //  First iteration norm
    double pNorm = 0.0;

    //  \rho_{infinity}
    double roInf = 0.0;

    //  Accepted relative tolerance
    double tol = 0.0;

    //  Equation symbol
    std::string sym;
    //std::string(LEN=2) :: sym = "NA";

    //  type of linear solver
    lsType ls;

    //  FSILS type of linear solver
    fsi_linear_solver::FSILS_lsType FSILS;

    //  BCs associated with this equation;
    std::vector<bcType> bc;

    //  IB: BCs associated with this equation on immersed surfaces
    std::vector<bcType> bcIB;

    //  domains that this equation must be solved
    std::vector<dmnType> dmn;

    //  IB: immersed domains that this equation must be solved
    std::vector<dmnType> dmnIB;

    //  Outputs
    std::vector<outputType> output;

    //  IB: Outputs
    std::vector<outputType> outIB;

    //  Body force associated with this equation
    std::vector<bfType> bf;
};

//----------
// dataType
//----------
// This type will be used to write data in the VTK files.
//
class dataType
{
  public:

    //  Element number of nodes
    int eNoN = 0;

    //  Number of elements
    int nEl = 0;

    //  Number of nodes
    int nNo = 0;

    //  vtk type
    int vtkType = 0;

    //  Connectivity array
    Array<int> IEN;

    //  Element based variables to be written
    Array<double> xe;

    //  All the variables after transformation to global format
    Array<double> gx;

    //  All the variables to be written (including position)
    Array<double> x;
};


class rmshType
{
  public:

    rmshType();

    // Whether remesh is required for problem or not
    bool isReqd = false;

    // Method for remeshing: 1-TetGen, 2-MeshSim
    consts::MeshGeneratorType method = consts::MeshGeneratorType::RMSH_TETGEN;

    // Counter to track number of remesh done
    int cntr = 0;

    // Time step from which remeshing is done
    int rTS = 0;

    // Time step freq for saving data
    int cpVar = 0;

    // Time step at which forced remeshing is done
    int fTS = 1000;

    // Time step frequency for forced remeshing
    int freq = 1000;

    // Time where remeshing starts
    double time = 0.0;

    // Mesh quality parameters
    double minDihedAng = 0.0;
    double maxRadRatio = 0.0;

    // Edge size of mesh
    Vector<double> maxEdgeSize;

    // Initial norm of an equation
    Vector<double> iNorm;

    // Copy of solution variables where remeshing starts
    Array<double> A0;
    Array<double> Y0;
    Array<double> D0;

    // Flag is set if remeshing is required for each mesh
    std::vector<bool> flag;
};

//------------
// ibCommType
//------------
//
class ibCommType
{
  public:
    //  Num traces (nodes) local to each process
    Vector<int> n;

    //  Pointer to global trace (node num) stacked contiguously
    Vector<int> gN;

    //  Num traces (Gauss points) local to each process
    Vector<int> nG;

    //  Pointer to global trace (Gauss point) stacked contiguously
    Vector<int> gE;
};


// Immersed Boundary (IB) data type
//
class ibType
{
  public:

    //  Whether any file being saved
    bool savedOnce = false;
    //bool savedOnce = .FALSE.

    //  IB method
    int mthd = 0;
    //int mthd = ibMthd_NA;

    //  IB coupling
    int cpld = 0;
    //int cpld = ibCpld_NA;

    //  IB interpolation method
    int intrp = 0;
    //int intrp = ibIntrp_NA;

    //  Current IB domain ID
    int cDmn = 0;

    //  Current equation
    int cEq = 0;

    //  Total number of IB nodes
    int tnNo = 0;

    //  Number of IB meshes
    int nMsh = 0;

    //  IB call duration (1: total time; 2: update; 3,4: communication)
    double callD[4] = {0.0, 0.0, 0.0, 0.0};

    //  IB Domain ID
    Vector<int> dmnID;

    //  Row pointer (for sparse LHS matrix storage)
    Vector<int> rowPtr;

    //  Column pointer (for sparse LHS matrix storage)
    Vector<int> colPtr;

    //  IB position coordinates
    Array<double> x;

    //  Velocity (new)
    Array<double> Yb;

    //  Time derivative of displacement (old)
    Array<double> Auo;

    //  Time derivative of displacement (new)
    Array<double> Aun;

    //  Time derivative of displacement (n+am)
    Array<double> Auk;

    //  Displacement (old)
    Array<double> Ubo;

    //  Displacement (new)
    Array<double> Ubn;

    //  Displacement (n+af)
    Array<double> Ubk;

    //  Displacement (projected on background mesh, old)
    Array<double> Uo;

    //  Displacement (projected on background mesh, new, n+af)
    Array<double> Un;

    //  Residue (FSI force)
    Array<double> R;

    //  Residue (displacement, background mesh)
    Array<double> Ru;

    //  Residue (displacement, IB mesh)
    Array<double> Rub;

    //  LHS tangent matrix for displacement
    Array<double> Ku;

    //  DERIVED class VARIABLES IB meshes;
    std::vector<mshType> msh;

    //  IB communicator
    ibCommType cm;
};

// Data type for Trilinos Linear Solver related arrays
//
class tlsType
{
  public:

    //  Local to global mapping
    Vector<int> ltg;

    //  Factor for Dirichlet BCs
    Array<double> W;

    //  Residue
    Array<double> R;
};


//--------
// ComMod 
//--------
// The ComMod class duplicates the data structures in the Fortran COMMOD module
// defined in MOD.f. 
//
// The data members here are the global variables exposed by the COMMOD module.
//
class ComMod {

  public:
    ComMod();
    ~ComMod();

    //----- bool members -----//

    // Whether there is a requirement to update mesh and Dn-Do variables
    bool dFlag = false;

    // Whether mesh is moving
    bool mvMsh = false;

    // Whether to averaged results
    bool saveAve = false;

    // Whether to save to VTK files
    bool saveVTK = false;

    // Whether any file being saved
    bool savedOnce = false;

    // Whether to use separator in output
    bool sepOutput = false;

    // Whether start from beginning or from simulations
    bool stFileFlag = false;

    // Whether to overwrite restart file or not
    bool stFileRepl = false;

    // Restart simulation after remeshing
    bool resetSim = false;

    // Check IEN array for initial mesh
    bool ichckIEN = false;

    // Reset averaging variables from zero
    bool zeroAve = false;

    // Whether CMM equation is initialized
    bool cmmInit = false;

    // Whether variable wall properties are used for CMM
    bool cmmVarWall = false;

    // Whether shell equation is being solved
    bool shlEq = false;

    // Whether PRESTRESS is being solved
    bool pstEq = false;

    // Whether velocity-pressure based structural dynamics solver is used
    bool sstEq = false;

    // Whether to detect and apply any contact model
    bool iCntct = false;

    // Whether any Immersed Boundary (IB) treatment is required
    bool ibFlag = false;

    // Postprocess step - convert bin to vtk
    bool bin2VTK = false;


    //----- int members -----//

    // Current domain
    int cDmn = 0;

    // Current equation
    int cEq = 0;

    // Current time step
    int cTS = 0;

    std::array<double,3> timeP;

    // Starting time step
    int startTS = 0;

    // Current equation degrees of freedom
    int dof = 0;

    // Global total number of nodes
    int gtnNo = 0;

    // Number of equations
    int nEq = 0;

    // Number of faces in the LHS passed to FSILS
    int nFacesLS = 0;

    // Number of meshes
    int nMsh = 0;

    // Number of spatial dimensions
    int nsd = 0;

    // Number of time steps
    int nTS = 0;

    // Number of initialization time steps
    int nITs = 0;

    // stFiles record length
    int recLn = 0;

    // Start saving after this number of time step
    int saveATS = 0;

    // Increment in saving solutions
    int saveIncr = 0;

    // Stamp ID to make sure simulation is compatible with stFiles
    std::array<int,7> stamp;

    // Increment in saving restart file
    int stFileIncr = 0;

    // Total number of degrees of freedom per node
    int tDof = 0;

    // Total number of nodes
    int tnNo = 0;

    // Restart Time Step
    int rsTS = 0;

    // Number of stress values to be stored
    int nsymd = 0;


    //----- double members -----//

    // Time step size
    double dt = 0.0;

    // Time
    double time = 0.0;


    //----- string members -----//

    // Initialization file path
    std::string iniFilePath;

    // Saved output file name
    std::string saveName;

    // Restart file name
    std::string stFileName;

    // Stop_trigger file name
    std::string stopTrigName;

    // ALLOCATABLE DATA

    // Column pointer (for sparse LHS matrix structure)
    // Modified in: lhsa()
    Vector<int> colPtr;

    // Domain ID
    Vector<int>  dmnId;

    // Local to global pointer tnNo --> gtnNo
    Vector<int> ltg;

    // Row pointer (for sparse LHS matrix structure)
    // Modified in: lhsa()
    Vector<int> rowPtr;

    // Array that maps global node id to rowN in the matrix
    // Modified in: lhsa()
    Vector<int> idMap;

    // Boundary nodes set for CMM initialization and for zeroing-out
    // non-wall nodal displacements
    Vector<int> cmmBdry;

    // IB: iblank used for immersed boundaries (1 => solid, 0 => fluid)
    Vector<int> iblank;

    // Old time derivative of variables (acceleration)
    Array<double>  Ao;

    // New time derivative of variables
    Array<double>  An;

    // Old integrated variables (dissplacement)
    Array<double>  Do;

    // New integrated variables
    Array<double>  Dn;

    // Residual vector
    Array<double>  R;

    // LHS matrix
    Array<double>  Val;

    // Position vector
    Array<double>  x;

    // Old variables (velocity)
    Array<double>  Yo;

    // New variables
    Array<double>  Yn;

    // Body force
    Array<double>  Bf;

    //-----------------------------------------------------
    // Additional arrays for velocity-based formulation of 
    // nonlinear solid mechanics.
    //-----------------------------------------------------

    // Time derivative of displacement
    Array<double>  Ad;

    // Residue of the displacement equation
    Array<double>  Rd;

    // LHS matrix for displacement equation
    Array<double>  Kd;

    // Variables for prestress calculations
    Array<double>  pS0;
    Array<double>  pSn;
    Vector<double>  pSa;

    // Temporary storage for initializing state variables
    Vector<double> Pinit;
    Array<double>  Vinit;
    Array<double>  Dinit;

    // CMM-variable wall properties: 1-thickness, 2-Elasticity modulus
    Array<double>  varWallProps;

    //------------------------
    // DERIVED TYPE VARIABLES
    //------------------------

    // Coupled BCs structures used for multidomain simulations
    cplBCType cplBC;

    // All data related to equations are stored in this container
    std::vector<eqType> eq;

    // FSILS data structure to produce LHS sparse matrix
    fsi_linear_solver::FSILS_lhsType lhs;

    // All the meshes are stored in this variable
    std::vector<mshType> msh;

    // Input/output to the screen is handled by this structure
    chnlType std, err, wrn, dbg;

    // To group above channels
    ioType io;

    // The general communicator
    cmType cm;

    // Remesher type
    rmshType rmsh;

    // Contact model type
    cntctModelType cntctM;

    // IB: Immersed boundary data structure
    ibType ib;

    // Trilinos Linear Solver data type
    tlsType  tls;

    bool debug_active = false;

    Timer timer;
};

#endif

