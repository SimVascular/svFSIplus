/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SIMULATION_H 
#define SIMULATION_H 

#include "ComMod.h"
#include "Parameters.h"
#include "SimulationLogger.h"

#include <string>

class Simulation {

  public:
    Simulation();
    ~Simulation();

    const mshType& get_msh(const std::string& name);

    CepMod& get_cep_mod() { return cep_mod; };
    ChnlMod& get_chnl_mod() { return chnl_mod; };
    ComMod& get_com_mod() { return com_mod; };

    // Read a solver paramerer input XML file.
    void read_parameters(const std::string& fileName);

    // Set simulation and module member data from Parameters.
    void set_module_parameters();

    //----- Fortran subroutines -----//
    //void read_msh();

    //----- Fortran modules -----//
    CepMod cep_mod;
    ChnlMod chnl_mod;
    CmMod cm_mod;
    ComMod com_mod;

    // Solver parameters read in from solver input XML file.
    Parameters parameters;

    // Log solution information.
    SimulationLogger logger;

    // Number of time steps
    int nTs;

    // Simulation initialization file path
    std::string fTmp;

    // Spectral radius of infinite time step; this is later used in equations.
    double roInf;

    // Simulation requires remeshing

    bool isReqd;

    // Name of the history file.
    std::string history_file_name;
};

#endif

