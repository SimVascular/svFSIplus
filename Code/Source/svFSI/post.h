#ifndef POST_H 
#define POST_H 

#include "Simulation.h"
#include "consts.h"

namespace post {

void all_post(Simulation* simulation, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputType outGrp, const int iEq);

void bpost(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputType outGrp);

void div_post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, const int iEq);

void fib_algn_post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lD, const int iEq);

void fib_dir_post(Simulation* simulation, const mshType& lM, const int nFn, Array<double>& res, const Array<double>& lD, const int iEq);

void fib_strech(Simulation* simulation, const int iEq, const mshType& lM, const Array<double>& lD, Vector<double>& res);

void post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputType outGrp, const int iEq);

void ppbin2vtk(Simulation* simulation);

void shl_post(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, 
    Vector<double>& resE, const Array<double>& lD, const int iEq, consts::OutputType outGrp);

void tpost(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, Vector<double>& resE, const Array<double>& lD, 
    const Array<double>& lY, const int iEq, consts::OutputType outGrp);

};

#endif

