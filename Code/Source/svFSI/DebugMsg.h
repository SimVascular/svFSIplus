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

#ifndef DEBUG_MSG_H 
#define DEBUG_MSG_H 

#include <iostream>

/// @brief The DebugMsg is class is used to print debugging messages. 
class DebugMsg 
{
  public:
    DebugMsg() = delete;
    DebugMsg(const char* function, const int task_id, bool add_endl=true) 
    {
      function_name_ = function;
      task_id_ = task_id;
      prefix_ = "[" + function_name_ + ":" + std::to_string(task_id) + "] ";
      banner_ = prefix_ + "=============== " + function_name_ + " ===============";
      add_endl_ = add_endl;
    }

    void banner() { std::cout << banner_ << std::endl; };
    std::string prefix() { return prefix_; };

    template <class T> DebugMsg& operator<< (const T& x) 
    { 
      if (start_) {
        std::cout << prefix_; 
      }

      std::cout << x;
      count_ += 1;

      if (add_endl_ && count_ == 2) {
        std::cout << std::endl; 
        start_ = true;
        count_ = 0;
      } else {
        start_ = false;
      }

      return *this; 
    };

    DebugMsg& operator<<(std::ostream& (*f)(std::ostream& o)) 
    {
      std::cout << f;
      start_ = true;
      count_ = 0;
      return *this;
    };

  private:
    bool add_endl_ = true;
    std::string banner_;
    int count_ = 0;
    std::string function_name_;
    std::string prefix_;
    bool start_ = true;
    int task_id_;
};

#endif

