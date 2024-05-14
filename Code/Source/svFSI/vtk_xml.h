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

#include "Simulation.h"
#include "ComMod.h"
#include "Array.h"

#ifndef VTK_XML_H
#define VTK_XML_H

namespace vtk_xml {

void do_test();

void int_msh_data(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, dataType& d, const int outDof, const int nOute);

void read_vtp(const std::string& file_name, faceType& face);

void read_vtp_pdata(const std::string& fName, const std::string& kwrd, const int nsd, const int m, const int idx, faceType& face);

void read_vtu(const std::string& file_name, mshType& mesh);

void read_vtu_pdata(const std::string& fName, const std::string& kwrd, const int nsd, const int m, const int idx, mshType& mesh);

void read_vtus(Simulation* simulation, Array<double>& lA, Array<double>& lY, Array<double>& lD, const std::string& fName);

void write_vtp(ComMod& com_mod, faceType& lFa, const std::string& fName);

void write_vtu(ComMod& com_mod, mshType& lM, const std::string& fName);

void write_vtu_debug(ComMod& com_mod, mshType& lM, const std::string& fName);

void write_vtus(Simulation* simulation, const Array<double>& lA, const Array<double>& lY, const Array<double>& lD, const bool lAve);

};

#endif


