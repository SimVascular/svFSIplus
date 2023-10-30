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

#include "Array3.h"
#include "utils.h"

template<>
double Array3<double>::memory_in_use = 0;

template<>
double Array3<double>::memory_returned = 0;

template<>
int Array3<double>::num_allocated = 0;

template<>
int Array3<double>::active = 0;

template<>
bool Array3<double>::write_enabled = false;

template<>
void Array3<double>::memory(const std::string& prefix)
{
  utils::print_mem("Array3<double>", prefix, memory_in_use, memory_returned);
}

template<>
void Array3<double>::stats(const std::string& prefix)
{
  utils::print_stats("Array3<double>", prefix, num_allocated, active);
}

//--------------------------//
//          int             //
//--------------------------//

template<>
double Array3<int>::memory_in_use = 0;

template<>
double Array3<int>::memory_returned = 0;

template<>
int Array3<int>::num_allocated = 0;

template<>
int Array3<int>::active = 0;

template<>
bool Array3<int>::write_enabled = false;
