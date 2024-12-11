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

#ifndef POST_H 
#define POST_H 

#include "Simulation.h"
#include "consts.h"

namespace post {

void all_post(Simulation* simulation, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputNameType outGrp, const int iEq);

void bpost(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputNameType outGrp);

void div_post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, const int iEq);

void fib_algn_post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lD, const int iEq);

void fib_dir_post(Simulation* simulation, const mshType& lM, const int nFn, Array<double>& res, const Array<double>& lD, const int iEq);

void fib_strech(Simulation* simulation, const int iEq, const mshType& lM, const Array<double>& lD, Vector<double>& res);

void post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputNameType outGrp, const int iEq);

void ppbin2vtk(Simulation* simulation);

void shl_post(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, 
    Vector<double>& resE, const Array<double>& lD, const int iEq, consts::OutputNameType outGrp);

void tpost(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, Vector<double>& resE, const Array<double>& lD, 
    const Array<double>& lY, const int iEq, consts::OutputNameType outGrp);

};

#endif

