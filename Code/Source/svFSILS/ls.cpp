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

#include "fsils_api.hpp"
#include "fils_struct.hpp"

namespace fsi_linear_solver {

/// @brief Set solver parameters.
///
/// Reproduces 'SUBROUTINE FSILS_LS_CREATE(ls, LS_type, relTol, absTol, maxItr, dimKry, relTolIn, absTolIn, maxItrIn)'.
//
void fsils_ls_create(FSILS_lsType& ls, LinearSolverType LS_type, double relTol, double absTol, int maxItr, int dimKry, 
  std::array<double,2> relTolIn, std::array<double,2> absTolIn, std::array<int,2> maxItrIn)
{
  using namespace consts;
  ls.foC = true;
  ls.LS_type = LS_type;

  // Set default parameters for each solver type.
  //
  switch (LS_type) {

    case LinearSolverType::LS_TYPE_NS:
      ls.RI.relTol = 0.4;
      ls.GM.relTol = 1.E-2;
      ls.CG.relTol = 0.2;
      ls.RI.mItr = 10;
      ls.GM.mItr = 2;
      ls.CG.mItr = 500;
      ls.GM.sD   = 100;
      ls.RI.sD   = 100;
    break;

    case LinearSolverType::LS_TYPE_GMRES:
      ls.RI.relTol = 0.1;
      ls.RI.mItr   = 4;
      ls.RI.sD     = 250;
    break;

    case LinearSolverType::LS_TYPE_CG:
      ls.RI.relTol = 1.E-2;
      ls.RI.mItr   = 1000;
    break;

    case LinearSolverType::LS_TYPE_BICGS:
      ls.RI.relTol = 1.E-2;
      ls.RI.mItr   = 500;
    break;

    default:
    break;
  }

  ls.RI.absTol = 1.E-10;
  ls.GM.absTol = 1.E-10;
  ls.CG.absTol = 1.E-10;

  // Set values optionally passed in.
  //
  if (present(relTol)) {
    ls.RI.relTol = relTol;
  }

  if (present(absTol)) {
    ls.RI.absTol = absTol;
  }

  if (present(maxItr)) {
    ls.RI.mItr = maxItr;
  }

  if (present(dimKry)) {
    ls.RI.sD = dimKry;
    ls.GM.sD = dimKry;
  }

  if (relTolIn.empty() != 0) {
    ls.GM.relTol = relTolIn[0];
    ls.CG.relTol = relTolIn[1];
  }

  if (absTolIn.empty() != 0) {
    ls.GM.absTol = absTolIn[0];
    ls.CG.absTol = absTolIn[1];
  }

  if (maxItrIn.empty() != 0) {
    ls.GM.mItr = maxItrIn[0];
    ls.CG.mItr = maxItrIn[1];
  }

}

};

