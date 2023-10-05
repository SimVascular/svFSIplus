
#include "Simulation.h"

#ifndef PIC_H
#define PIC_H

namespace pic {

void picc(Simulation* simulation);

void pic_eth(Simulation* simulation);

void pici(Simulation* simulation, Array<double>& Ag, Array<double>& Yg,
          Array<double>& Dg);

void picp(Simulation* simulation);

};  // namespace pic

#endif
