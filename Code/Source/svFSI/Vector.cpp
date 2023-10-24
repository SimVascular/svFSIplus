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
#include "utils.h"

template<>
double Vector<double>::memory_in_use = 0;

template<>
double Vector<double>::memory_returned = 0;

template<>
int Vector<double>::num_allocated = 0;

template<>
int Vector<double>::active = 0;

template<>
bool Vector<double>::write_enabled = false;

template<>
void Vector<double>::memory(const std::string& prefix)
{
  utils::print_mem("Vector<double>", prefix, memory_in_use, memory_returned);
}

template<>
void Vector<double>::stats(const std::string& prefix)
{
  utils::print_stats("Vector<double>", prefix, num_allocated, active);
}

//------//
//  int  //
//------//

template<>
double Vector<int>::memory_in_use = 0;

template<>
double Vector<int>::memory_returned = 0;

template<>
int Vector<int>::num_allocated = 0;

template<>
int Vector<int>::active = 0;

template<>
bool Vector<int>::write_enabled = false;

template<>
void Vector<int>::memory(const std::string& prefix)
{
  utils::print_mem("Vector<int>", prefix, memory_in_use, memory_returned);
}

template<>
void Vector<int>::stats(const std::string& prefix)
{
  utils::print_stats("Vector<int>", prefix, num_allocated, active);
}

// Vector<Vector<double>> 
template<>
double Vector<Vector<double>>::memory_in_use = 0;

template<>
double Vector<Vector<double>>::memory_returned = 0;

template<>
int Vector<Vector<double>>::num_allocated = 0;

template<>
int Vector<Vector<double>>::active = 0;

template<>
bool Vector<Vector<double>>::write_enabled = false;

// float //
template<>
double Vector<float>::memory_in_use = 0;

template<>
double Vector<float>::memory_returned = 0;

template<>
int Vector<float>::num_allocated = 0;

template<>
int Vector<float>::active = 0;

template<>
bool Vector<float>::write_enabled = false;

/// @brief Build a prefix for a file name from label which may be from a debugging message.
//
std::string build_file_prefix(const std::string& label)
{
  std::string file_prefix;

  for (auto c : label) {
    if (c == '[') {
      continue;
    }
    if (c == ']') {
      file_prefix.push_back('_');
      continue;
    }
    if (c == ' ') {
      continue;
    }
    if (c == ':') {
      c =  '_';
    }
    file_prefix.push_back(c);
  }

  return file_prefix;
}

