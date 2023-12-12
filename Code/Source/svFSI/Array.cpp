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

#include "Array.h"

// b o o l //

template<>
bool Array<bool>::show_index_check_message = true;

template<>
int Array<bool>::id = 0;

template<>
double Array<bool>::memory_in_use = 0;

template<>
double Array<bool>::memory_returned = 0;

template<>
int Array<bool>::num_allocated = 0;

template<>
int Array<bool>::active = 0;

template<>
bool Array<bool>::write_enabled = false;

// d o u b l e //

template<>
bool Array<double>::show_index_check_message = true;

template<>
int Array<double>::id = 0;

template<>
double Array<double>::memory_in_use = 0;

template<>
double Array<double>::memory_returned = 0;

template<>
int Array<double>::num_allocated = 0;

template<>
int Array<double>::active = 0;

template<>
void Array<double>::memory(const std::string& prefix)
{
  utils::print_mem("Array<double>", prefix, memory_in_use, memory_returned);
}

template<>
void Array<double>::stats(const std::string& prefix)
{ 
  utils::print_stats("Array<double>", prefix, num_allocated, active);
}

template<>
bool Array<double>::write_enabled = false;

//  i n t  //

template<>
bool Array<int>::show_index_check_message = true;

template<>
int Array<int>::id = 0;

template<>
double Array<int>::memory_in_use = 0;

template<>
double Array<int>::memory_returned = 0;

template<>
int Array<int>::num_allocated = 0;

template<>
int Array<int>::active = 0;

template<>
void Array<int>::memory(const std::string& prefix)
{
  utils::print_mem("Array<int>", prefix, memory_in_use, memory_returned);
}

template<>
void Array<int>::stats(const std::string& prefix)
{
  utils::print_stats("Array<int>", prefix, num_allocated, active);
}

template<>
bool Array<int>::write_enabled = false;

