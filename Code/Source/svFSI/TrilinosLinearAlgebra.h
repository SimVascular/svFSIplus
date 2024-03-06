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

#ifndef TRILINOS_LINEAR_ALGEBRA_H 
#define TRILINOS_LINEAR_ALGEBRA_H 

#include "LinearAlgebra.h"

/// @brief The TrilinosLinearAlgebra class implements the LinearAlgebra 
/// interface for the Trilinos numerical linear algebra package.
///
class TrilinosLinearAlgebra : public virtual LinearAlgebra {

  public:
    TrilinosLinearAlgebra();
    ~TrilinosLinearAlgebra();

    virtual void alloc(ComMod& com_mod, eqType& lEq);
    virtual void assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR);
    virtual void check_options(const consts::PreconditionerType prec_cond_type, const consts::LinearAlgebraType assembly_type);
    virtual void initialize(ComMod& com_mod, eqType& lEq);
    virtual void set_assembly(consts::LinearAlgebraType atype);
    virtual void set_preconditioner(consts::PreconditionerType prec_type);
    virtual void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);

  private:
    static std::set<consts::LinearAlgebraType> valid_assemblers;
    void initialize_fsils(ComMod& com_mod, eqType& lEq);
    bool use_fsils_assembly = false;
    /// @brief The FsilsLinearAlgebra object used to assemble local element matrices.
    LinearAlgebra* fsils_solver = nullptr;
    // Private class used to hide Trilinos implementation details.
    class TrilinosImpl;
    TrilinosImpl* impl = nullptr;
};

#endif

