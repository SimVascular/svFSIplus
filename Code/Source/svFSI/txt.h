#ifndef TXT_H 
#define TXT_H 

#include "Simulation.h"
#include "Array.h"
#include "ComMod.h"

#include "consts.h"

namespace txt_ns {

void cctxt(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const std::array<std::string,2>& fName, const std::vector<bool>& wtn);

void txt(Simulation* simulation, const bool flag);

void wtxt(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const int m, const std::array<std::string,2>& fName, 
    const Array<double>& tmpV, const std::vector<bool>& wtn, const bool div, const bool pflag);

};

#endif

