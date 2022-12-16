
#ifndef FSI_LINEAR_SOLVER_H 
#define FSI_LINEAR_SOLVER_H 

#include "CmMod.h"

#include "mpi.h"

#include <map>

/*
SELECTED_REAL_KIND(P,R) returns the kind value of a real data type with 

  1) decimal precision of at least P digits, 

  2) exponent range of at least R

INTEGER, PARAMETER :: LSRP4  = SELECTED_REAL_KIND(6, 37)
INTEGER, PARAMETER :: LSRP8  = SELECTED_REAL_KIND(15, 307)      // C++ double
INTEGER, PARAMETER :: LSRP16 = SELECTED_REAL_KIND(33, 4931)     // C++ long double
INTEGER, PARAMETER :: LSRP = LSRP8
*/

// The enums here replicate the PARAMETERs defined
// in FSILS_STRUCT.h.

namespace fsi_linear_solver {

enum class BcType
{
  BC_TYPE_Dir = 0,
  BC_TYPE_Neu = 1
};

enum class BcopType
{
  BCOP_TYPE_ADD = 0,
  BCOP_TYPE_PRE = 1
};

//------------------
// LinearSolverType 
//------------------
//
enum LinearSolverType
{
  LS_TYPE_CG = 798,
  LS_TYPE_GMRES = 797, 
  LS_TYPE_NS = 796, 
  LS_TYPE_BICGS = 795
};

//-----------------
// FSILS_commuType
//-----------------
//
class FSILS_commuType 
{
  public:
    int foC;             // Free of created          (USE)
    int masF;            // If this the master       (USE)
    int master;          // Master ID                (USE)
    int task;            // ID of this proc.         (USE)
    int tF;              // Task in FORTRAN indexing (USE)
    int nTasks;          // Total number of tasks    (USE)
    MPI_Comm comm;       // MPI communicator         (IN)
};

//--------------
// FSILS_cSType
//--------------
//
class FSILS_cSType
{
  public:
    // The processor to communicate with
    int iP = -1;

    // Number of data to be commu
    int n = 0;

    // Pointer to the data for commu
    Vector<int> ptr;
};

//----------------
// FSILS_faceType
//----------------
//
class FSILS_faceType
{
  public:
    // Free or created                (USE)
    bool foC = false;

    // Neu: P/Q coupling              (USE)
    bool coupledFlag = false;

    // Neu: shared between proces     (USE)
    bool sharedFlag = false;

    // Included in the computations   (IN)
    bool incFlag = false;

    // Number of nodes                (IN)
    int nNo = 0;

    // Degrees of freedom for val     (IN)
    int dof = 0;

    // Dir/Neu                        (IN)
    BcType bGrp = BcType::BC_TYPE_Dir;

    // Only for data alignment
    int reserved = 0;

    // Global node number             (IN)
    Vector<int> glob;

    // ||Sai||**2._LSRP                   (USE)
    double nS = 0.0;

    // Neu: P = res*Q                 (IN)
    double res = 0.0;

    // nodal Sai for Neu              (IN)
    Array<double> val;

    // Neu W*Sai                      (TMP)
    Array<double> valM;
};

//---------------
// FSILS_lhsType
//---------------
//
// Modified in:
//
//  fsils_lhs_create()
//
class FSILS_lhsType 
{
  public:
    // Free of created                     (USE)
    bool foC = false; 

    bool debug_active = false; 

    // Global number of nodes              (IN)
    int gnNo = 0;

    // Number of nodes                     (IN)
    int nNo = 0;

    // Number of non-zero in lhs           (IN)
    int nnz = 0;

    //  Number of faces                     (IN)
    int nFaces = 0;

    // nNo of this proc                    (USE)
    int mynNo = 0;

    // nNo of shared with lower proc       (USE)
    int shnNo = 0;

    // Number of communication requests    (USE)
    int nReq = 0;

    // Column pointer                      (USE)
    Vector<int> colPtr;

    // Row pointer                         (USE)
    Array<int> rowPtr;

    // Diagonal pointer                    (USE)
    Vector<int> diagPtr;

    // Mapping of nodes                    (USE)
    Vector<int> map;

    FSILS_commuType commu;

    std::vector<FSILS_cSType> cS;

    std::vector<FSILS_faceType> face;
};

//-----------------
// FSILS_subLsType
//-----------------
//
class FSILS_subLsType 
{
  public:
    bool suc;       // Successful solving          (OUT)
    //int suc;       // Successful solving          (OUT)
    int mItr;      // Maximum iteration           (IN)
    int sD;        // Space dimension             (IN)
    int itr;       // Number of iteration         (OUT)
    int cM;        // Number of Ax multiply       (OUT)
    int cN;        // Number of |x| norms         (OUT)
    int cD;        // Number of <x.y> dot products(OUT)
    int reserve;   // Only for data alignment     (-)
    double absTol; // Absolute tolerance          (IN)
    double relTol; // Relative tolerance          (IN)
    double iNorm;  // Initial norm of residual    (OUT)
    double fNorm;  // Final norm of residual      (OUT)
    double dB;     // Res. rduction in last itr.  (OUT)
    double callD;  // Calling duration            (OUT)
};

class FSILS_lsType 
{
  public:
    int foC;                     // Free of created             (USE)
    LinearSolverType LS_type;    // Which one of LS             (IN)
    int Resm;                    // Contribution of mom. res.   (OUT)
    int Resc;                    // Contribution of cont. res.  (OUT)
    FSILS_subLsType GM;
    FSILS_subLsType CG;
    FSILS_subLsType RI;
};


};

#endif

