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

#include "PetscLinearAlgebra.h"

#include <iostream>

// Include PETSc-dependent data structures and functions. 
//
#ifdef USE_PETSC
#include "petsc_impl.cpp"

// If PETSc is not used then define PetscImpl with noop methods.
//
#else
class PetscLinearAlgebra::PetscImpl {
  public:
    PetscImpl(){};
    void initialize(ComMod& com_mod, eqType& lEq) {};
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
};
#endif

/////////////////////////////////////////////////////////////////
//             P e t s c L i n e a r A l g e b r a             //
/////////////////////////////////////////////////////////////////
// The following methods implement the LinearAlgebra interface.

PetscLinearAlgebra::PetscLinearAlgebra()
{
  #ifndef USE_PETSC
  throw std::runtime_error("[PetscLinearAlgebra] There is no PETSc interface.");
  #else
  impl = new PetscLinearAlgebra::PetscImpl();
  interface_type = consts::LinearAlgebraType::petsc; 
  assembly_type = consts::LinearAlgebraType::none;
  preconditioner_type = consts::PreconditionerType::PREC_NONE;
  #endif
}

void PetscLinearAlgebra::alloc(ComMod& com_mod, eqType& lEq)
{
  initialize_fsils(com_mod, lEq);
}

//----------
// assemble
//----------
//
void PetscLinearAlgebra::assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN, 
    const Array3<double>& lK, const Array<double>& lR)
{
  fsils_solver->assemble(com_mod, num_elem_nodes, eqN, lK, lR);
}

//------------
// initialize
//------------
//
void PetscLinearAlgebra::initialize(ComMod& com_mod, eqType& lEq)
{
  impl->initialize(com_mod, lEq);
}

//------------------
// initialize_fsils
//------------------
//
void PetscLinearAlgebra::initialize_fsils(ComMod& com_mod, eqType& lEq)
{
  fsils_solver = LinearAlgebraFactory::create_interface(consts::LinearAlgebraType::fsils);
  fsils_solver->initialize(com_mod, lEq);
  fsils_solver->alloc(com_mod, lEq);
  fsils_solver->set_preconditioner(preconditioner_type);
}

//--------------
// set_assembly
//--------------
//
void PetscLinearAlgebra::set_assembly(consts::LinearAlgebraType atype)
{
  assembly_type = atype;
}

void PetscLinearAlgebra::set_preconditioner(consts::PreconditionerType prec_type)
{
  preconditioner_type = prec_type;
}

void PetscLinearAlgebra::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  std::cout << "[PetscLinearAlgebra] solve" << std::endl;
  impl->solve(com_mod, lEq, incL, res);
}

void PetscLinearAlgebra::solve_assembled(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  impl->solve(com_mod, lEq, incL, res);
}

