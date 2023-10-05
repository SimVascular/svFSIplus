#ifndef OUTPUT__H
#define OUTPUT__H

#include <fstream>
#include <iostream>

#include "Simulation.h"

namespace output {

void output_result(Simulation* simulation, std::array<double, 3>& timeP,
                   const int co, const int iEq);

void read_restart_header(ComMod& com_mod, std::array<int, 7>& tStamp,
                         double& timeP, std::ifstream& restart_file);

void write_restart(Simulation* simulation, std::array<double, 3>& timeP);

void write_restart_header(ComMod& com_mod, std::array<double, 3>& timeP,
                          std::ofstream& restart_file);

void write_results(ComMod& com_mod, const std::array<double, 3>& timeP,
                   const std::string& fName, const bool sstEq);

};  // namespace output

#endif
