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

#include "TrilinosLinearAlgebra.h"
#include <iostream>

// Include Trilinos-dependent data structures and functions. 
//
#ifdef WITH_TRILINOS
#include "trilinos_impl.cpp"

// If Trilinos is not used then define TrilinosImpl with noop methods.
//
#else
class TrilinosLinearAlgebra::TrilinosImpl {
  public:
    TrilinosImpl(){};
    void alloc(ComMod& com_mod, eqType& lEq){};
    void initialize(ComMod& com_mod) {};
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
};
#endif

/////////////////////////////////////////////////////////////////
//          T r i l i n o s L i n e a r A l g e b r a          //
/////////////////////////////////////////////////////////////////
// The following methods implement the LinearAlgebra interface.

TrilinosLinearAlgebra::TrilinosLinearAlgebra()
{
  //std::cout << "[TrilinosLinearAlgebra] ---------- TrilinosLinearAlgebra ---------- " << std::endl;
  #ifndef WITH_TRILINOS
  throw std::runtime_error("[TrilinosLinearAlgebra] There is no Trilinos interface.");
  #else
  impl = new TrilinosLinearAlgebra::TrilinosImpl();
  interface_type = LinearAlgebraType::petsc; 
  #endif
}

void TrilinosLinearAlgebra::alloc(ComMod& com_mod, eqType& lEq)
{
  //std::cout << "[TrilinosLinearAlgebra] ---------- alloc ---------- " << std::endl;
  impl->alloc(com_mod, lEq);
}

void TrilinosLinearAlgebra::initialize(ComMod& com_mod)
{
  //std::cout << "[TrilinosLinearAlgebra] ---------- initialize ---------- " << std::endl;
  impl->initialize(com_mod);
}

void TrilinosLinearAlgebra::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  std::cout << "[TrilinosLinearAlgebra] ---------- solve ---------- " << std::endl;
  impl->solve(com_mod, lEq, incL, res);
}

