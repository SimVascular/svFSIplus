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

// Subroutines related to initializing linear solver arrays and
// function calls to svFSILS and Trilinos solver library

#include "ls.h"

#include "fsils_api.hpp"
#include "consts.h"

#include <math.h>

namespace ls_ns {

/// @brief Allocate com_mod.R and com_mod.Val arrays.
///
/// Modifies:
///    com_mod.R - Residual vector
///    com_mod.Val - LHS matrix 
///
/// Reproduces 'SUBROUTINE LSALLOC(lEq)'.
//
void ls_alloc(ComMod& com_mod, eqType& lEq)
{
  int dof = com_mod.dof;
  int tnNo = com_mod.tnNo;
  int gtnNo = com_mod.gtnNo;
  auto& lhs = com_mod.lhs;

  com_mod.R.resize(dof,tnNo);

  lEq.linear_algebra->alloc(com_mod, lEq);
}

/// @brief Modifies:    
///  com_mod.R      // Residual vector
///  com_mod.Val    // LHS matrix
///
/// Reproduces ' SUBROUTINE LSSOLVE(lEq, incL, res)'.
//
void ls_solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) 
{
  #define n_debug_ls_solve
  #ifdef debug_ls_solve 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lEq.sym: " << lEq.sym;
  dmsg << "lEq.useTLS: " << lEq.useTLS;
  dmsg << "lEq.assmTLS: " << lEq.assmTLS;
  #endif

  lEq.linear_algebra->solve(com_mod, lEq, incL, res);
}

};


