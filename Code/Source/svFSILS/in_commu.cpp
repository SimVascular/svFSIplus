/**
 * Copyright (c) Stanford University, The Regents of the University of California, and others.
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

// The functions here reproduce the subroutines defined in svFSILS/INCOMMU.f.

// To syncronize the data on the boundaries between the processors.
//
// This a both way communication with three main part:
// 1 - rTmp {in master} = R          {from slave}
// 2 - R    {in master} = R + rTmp   {both from master}
// 3 - rTmp {in master} = R          {from master}
// 4 - R    {in slave}  = rTmp       {from master}

#include "fsils.hpp"
#include "CmMod.h"
#include "Array3.h"

#include "fsils_std.h"

namespace fsi_linear_solver {

void fsils_commus(const FSILS_lhsType& lhs, Vector<double>& R)
{
  if (lhs.commu.nTasks == 1) {
    return;
  }

  int nReq = lhs.nReq;
  if (lhs.cS.size() == 0) {
    return;
  }

  int nmax = std::max_element(lhs.cS.begin(), lhs.cS.end(), 
    [](const FSILS_cSType& a, const FSILS_cSType& b){return a.n < b.n;})->n; 

  Array<double> sB(nmax,nReq); 
  Array<double> rB(nmax,nReq); 
  std::vector<MPI_Request> rReq(nReq); 
  std::vector<MPI_Request> sReq(nReq);

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      sB(j,i) = R(k);
    }
  }

  int mpi_tag = 1;

  for (int i = 0; i < nReq; i++) {
    auto rec_err = MPI_Irecv(rB.col_data(i), lhs.cS[i].n, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &rReq[i]);
    auto send_err = MPI_Isend(sB.col_data(i), lhs.cS[i].n, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &sReq[i]);
  }

  // Wait for the MPI receive to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&rReq[i], &stat);
  }

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      R(k) = R(k) + rB(j,i);
    }
  }

  // Wait for the MPI send to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&sReq[i], &stat);
  }
}

/// @brief This a both way communication with three main part:
///
/// 1 - rTmp {in master} = R          {from slave}
/// 2 - R    {in master} = R + rTmp   {both from master}
/// 3 - rTmp {in master} = R          {from master}
/// 4 - R    {in slave}  = rTmp       {from master}
//
void fsils_commuv(const FSILS_lhsType& lhs, int dof, Array<double>& R)
{
  if (lhs.commu.nTasks == 1) {
    return;
  }

  if (lhs.cS.size() == 0) {
    return;
  }

  int nReq = lhs.nReq;
  int nmax = std::max_element(lhs.cS.begin(), lhs.cS.end(), 
      [](const FSILS_cSType& a, const FSILS_cSType& b){return a.n < b.n;})->n; 

  Array3<double> sB(dof,nmax,nReq); 
  Array3<double> rB(dof,nmax,nReq); 
  std::vector<MPI_Request> rReq(nReq); 
  std::vector<MPI_Request> sReq(nReq);

  int n = 1;

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      for (int l = 0; l < dof; l++) { 
        sB(l,j,i) = R(l,k);
      }
    }
  }

  int mpi_tag = 1;

  for (int i = 0; i < nReq; i++) {
    auto rec_err = MPI_Irecv(rB.slice_data(i), lhs.cS[i].n*dof, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &rReq[i]);
    auto send_err = MPI_Isend(sB.slice_data(i), lhs.cS[i].n*dof, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &sReq[i]);
  }

  // Wait for the MPI receive to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&rReq[i], &stat);
  }

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      for (int l = 0; l < dof; l++) { 
        R(l,k) = R(l,k) + rB(l,j,i);
      }
    }
  }

  // Wait for the MPI send to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&sReq[i], &stat);
  }
}

};


