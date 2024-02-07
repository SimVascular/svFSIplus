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

#include "PetscLinearSolver.h"

#include <iostream>

/////////////////////////////////////////////////////////////////
//                    P e t s c I m p l                        //
/////////////////////////////////////////////////////////////////

class PetscLinearSolver::PetscImpl {
  public:
    PetscImpl();
    #ifdef USE_PETSC
    void initialize(ComMod& com_mod);
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);
    #else
    void initialize(ComMod& com_mod) {};
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
    #endif

    void init_dir_and_coupneu_bc_petsc(ComMod& com_mod, const Vector<int>& incL, const Vector<double>& res);

    // Local to global mapping
    Vector<int> ltg_;

    // Factor for Dirichlet BCs
    Array<double> W_;

    // Residue
    Array<double> R_;

    // Factor for Lumped Parameter BCs
    Array<double> V_;
};

#ifdef USE_PETSC
#include "petsc_interface.cpp"
#endif

PetscLinearSolver::PetscImpl::PetscImpl() 
{
  std::cout << "[PetscImpl] ---------- PetscImpl() ---------- " << std::endl;
}

/////////////////////////////////////////////////////////////////
//             P e t s c L i n e a r S o l v e r               //
/////////////////////////////////////////////////////////////////

PetscLinearSolver::PetscLinearSolver()
{
  std::cout << "[PetscLinearSolver] ---------- PetscLinearSolver ---------- " << std::endl;
  #ifndef USE_PETSC
  throw std::runtime_error("[PetscLinearSolver] There is no PETSc interface.");
  #else
  impl = new PetscLinearSolver::PetscImpl();
  solver_type = LinearSolverInterface::petsc; 
  #endif
}

void PetscLinearSolver::initialize(ComMod& com_mod)
{
  std::cout << "[PetscLinearSolver] ---------- initialize ---------- " << std::endl;
  impl->initialize(com_mod);
}

void PetscLinearSolver::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  std::cout << "[PetscLinearSolver] ---------- solve ---------- " << std::endl;
  impl->solve(com_mod, lEq, incL, res);
}

