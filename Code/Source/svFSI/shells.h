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

#ifndef SHELLS_H 
#define SHELLS_H 

#include "ComMod.h"
#include "Simulation.h"

namespace shells {

void construct_shell(ComMod& com_mod, const mshType& lM, const Array<double>& Ag,
    const Array<double>& Yg, const Array<double>& Dg);

void shell_3d(ComMod& com_mod, const mshType& lM, const int g, const int eNoN,
    const int nFn, const Array<double>& fN, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& xl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK);

void shell_bend_cst(ComMod& com_mod, const mshType& lM, const int e, const Vector<int>& ptr,
    Array<double>& x0, Array<double>& xc, double bb_0[2][2], double bb_x[2][2],
    Array3<double>& Bb, const bool vflag);

void shell_bf(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& dl, const Array<double>& xl, const Array<double>& tfl, Array<double>& lR, Array3<double>& lK);

void shell_cst(ComMod& com_mod, const mshType& lM, const int e, const int eNoN, const int nFn, const Array<double>& fN,
    const Array<double>& al, const Array<double>& yl, const Array<double>& dl, const Array<double>& xl,
    const Array<double>& bfl, const Vector<int>& ptr);

void shell_fp(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& dl, const Array<double>& xl, const Array<double>& tfl, Array<double>& lR, Array3<double>& lK);

void shl_strs_res(const ComMod& com_mod, const dmnType& lDmn, const int nFn, const Array<double>& fNa0,
    const double aa_0[2][2], const double aa_x[2][2], const double bb_0[2][2], const double bb_x[2][2],
    double& lam3, Array<double>& Sm, Array3<double>& Dm);

};

#endif

