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

// The functions here reproduce the subroutines defined in svFSILS/COMMU.f.

#include "fsils.hpp"
#include "CmMod.h"

#include "mpi.h"

namespace fsi_linear_solver {

/// @brief Modifies:
///   commu.task
///   commu.nTasks
//   
void fsils_commu_create(FSILS_commuType& commu, cm_mod::MpiCommWorldType commi)
{
  // Some of these parameters are set for sequential version
  //
  commu.foC = true;
  commu.comm = commi;
  commu.nTasks = 1;
  commu.task   = 0;
  commu.master = 0;

  MPI_Comm_rank(commi, &commu.task);
  MPI_Comm_size(commi, &commu.nTasks);

  MPI_Allreduce(&commu.task, &commu.master, 1, cm_mod::mpint, MPI_MIN, commi);

  commu.masF = false;
  commu.tF = commu.task;
  if (commu.task == commu.master) {
    commu.masF = true; 
  }
}


void fsils_commu_free(FSILS_commuType& commu)
{
  //IF (.NOT.commu%foC) STOP 'COMMU is not created yet to be freed'

  commu.foC = false;
}

};


