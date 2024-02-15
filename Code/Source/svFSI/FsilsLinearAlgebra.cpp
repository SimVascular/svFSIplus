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

#include "FsilsLinearAlgebra.h"
#include "ComMod.h"
#include "fsils_api.hpp"
#include "lhsa.h"
#include <iostream>

// Include FSILS-dependent data structures and functions. 
//
//#ifdef USE_PETSC
//#include "petsc_impl.cpp"

// If PETSc is not used then define FsilsImpl with noop methods.
//
//#else
//class FsilsLinearAlgebra::FsilsImpl {
  //public:
    //FsilsImpl(){};
    //void initialize(ComMod& com_mod) {};
    //void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
//};
//#endif

class FsilsLinearAlgebra::FsilsImpl {
  public:
    FsilsImpl(){};
    void initialize(ComMod& com_mod) {};
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
};

/////////////////////////////////////////////////////////////////
//             F s i l s L i n e a r A l g e b r a             //
/////////////////////////////////////////////////////////////////
// The following methods implement the LinearAlgebra interface.

FsilsLinearAlgebra::FsilsLinearAlgebra()
{
  //std::cout << "[FsilsLinearAlgebra] ---------- FsilsLinearAlgebra ---------- " << std::endl;
  impl = new FsilsLinearAlgebra::FsilsImpl();
  interface_type = consts::LinearAlgebraType::fsils; 
  assembly_type = consts::LinearAlgebraType::fsils; 
  preconditioner_type = consts::PreconditionerType::PREC_FSILS;
}

void FsilsLinearAlgebra::alloc(ComMod& com_mod, eqType& lEq)
{
  //std::cout << "[FsilsLinearAlgebra::alloc] ---------- FsilsLinearAlgebra ---------- " << std::endl;
  int dof = com_mod.dof;
  int tnNo = com_mod.tnNo;
  int gtnNo = com_mod.gtnNo;
  auto& lhs = com_mod.lhs;
  //std::cout << "[FsilsLinearAlgebra::alloc] com_mod.R.size(): " << com_mod.R.size() << std::endl;
  //std::cout << "[FsilsLinearAlgebra::alloc] dof: " << dof << std::endl;
  //std::cout << "[FsilsLinearAlgebra::alloc] com_mod.lhs.nnz: " << com_mod.lhs.nnz << std::endl;

  if (!use_trilinos_assembly) {
    com_mod.Val.resize(dof*dof, com_mod.lhs.nnz);
  }

  if (use_trilinos_preconditioner) { 
    initialize_trilinos(com_mod, lEq);
  }

}

//----------
// assemble
//----------
//
void FsilsLinearAlgebra::assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR)
{
  #if debug_assemble
  std::cout << "[FsilsLinearAlgebra::assemble] ---------- assemble ---------- " << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] num_elem_nodes: " << num_elem_nodes << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] eqN.size(): " << eqN.size() << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] lK.size(): " << lK.size() << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] lR.size(): " << lR.size() << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] assembly_set: " << assembly_set << std::endl;
  #endif

  if (assembly_set) {
    trilinos_solver->assemble(com_mod, num_elem_nodes, eqN, lK, lR);
  } else {
    lhsa_ns::do_assem(com_mod, num_elem_nodes, eqN, lK, lR);
  }
}

void FsilsLinearAlgebra::initialize(ComMod& com_mod)
{
  //std::cout << "[FsilsLinearAlgebra] ---------- initialize ---------- " << std::endl;
  //impl->initialize(com_mod);
}

//---------------------
// initialize_trilinos
//---------------------
//
void FsilsLinearAlgebra::initialize_trilinos(ComMod& com_mod, eqType& lEq)
{
  trilinos_solver = LinearAlgebraFactory::create_interface(consts::LinearAlgebraType::trilinos);
  trilinos_solver->initialize(com_mod);
  trilinos_solver->alloc(com_mod, lEq);
  trilinos_solver->set_preconditioner(preconditioner_type);
}

void FsilsLinearAlgebra::set_assembly(consts::LinearAlgebraType atype)
{
  //std::cout << "[FsilsLinearAlgebra::set_assembly] ---------- set_assembly ---------- " << std::endl;

  if (atype == consts::LinearAlgebraType::none) {
    std::cout << "[FsilsLinearAlgebra::set_assembly] atype == LinearAlgebraType::none" << std::endl;
    return;
  }

  assembly_set = true;
  assembly_type = atype;

  if (assembly_type == consts::LinearAlgebraType::trilinos) { 
    use_trilinos_assembly = true;
  }
}

//--------------------
// set_preconditioner
//--------------------
//
void FsilsLinearAlgebra::set_preconditioner(consts::PreconditionerType prec_type)
{
  //std::cout << "[FsilsLinearAlgebra] ---------- set_preconditioner ---------- " << std::endl;
  preconditioner_type = prec_type;

  if (consts::trilinos_preconditioners.count(prec_type) != 0) {
    std::cout << "[FsilsLinearAlgebra] set Trilinos preconditioner " << std::endl;
    use_trilinos_preconditioner = true;
    //trilinos_solver = LinearAlgebraFactory::create_interface(LinearAlgebraType::trilinos);
  }
}

void FsilsLinearAlgebra::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  //std::cout << "[FsilsLinearAlgebra] ---------- solve ---------- " << std::endl;
  auto& lhs = com_mod.lhs;
  int dof = com_mod.dof;
  auto& R = com_mod.R;      // Residual vector
  auto& Val = com_mod.Val;  // LHS matrix
  //std::cout << "[FsilsLinearAlgebra] dof: " << dof << std::endl;
  //std::cout << "[FsilsLinearAlgebra] R.size(): " << R.size() << std::endl;
  //std::cout << "[FsilsLinearAlgebra] Val.size(): " << Val.size() << std::endl;

  if (use_trilinos_preconditioner) { 

    if (use_trilinos_assembly) {
      trilinos_solver->solve_assembled(com_mod, lEq, incL, res);

    } else {
      trilinos_solver->solve(com_mod, lEq, incL, res);
    }

  } else {

    fsi_linear_solver::fsils_solve(lhs, lEq.FSILS, dof, R, Val, lEq.ls.PREC_Type, incL, res);

  }
}

void FsilsLinearAlgebra::solve_assembled(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  solve(com_mod, lEq, incL, res);
}

