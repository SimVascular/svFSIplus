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

#ifndef SIMULATION_LOGGER_H 
#define SIMULATION_LOGGER_H 

#include <fstream>
#include <string>

/// @brief The SimulationLogger class is used to write information 
/// to a text file and optionally to cout.
//
class SimulationLogger {

  public:
    SimulationLogger() { }

    SimulationLogger(const std::string& file_name, bool cout_write=false)
    { 
      this->initialize(file_name, cout_write);
    }

    void initialize(const std::string& file_name, bool cout_write=false) 
    {
      log_file_.open(file_name);
      if (log_file_.fail()) {
        throw std::runtime_error("[SimulationLogger] Unable to open the file '" + file_name + "' for writing.");
      }

      cout_write_ = cout_write;
      file_name_ = file_name;
    }

    ~SimulationLogger() 
    {
      log_file_.close();
    }

  template <class T> SimulationLogger& operator<< (const T& value)
  {
    if (file_name_ == "") { 
      return *this;
    }

    log_file_ << value;

    if (cout_write_) {
      std::cout << value;
    }

    return *this;
  }

  SimulationLogger& operator<<(std::ostream&(*value)(std::ostream& o))
  {
    if (file_name_ == "") { 
      return *this;
    }

    log_file_ << value;

    if (cout_write_) {
      std::cout << value;
    }

    return *this;
  };

  private:
    bool cout_write_ = false;
    std::string file_name_;
    std::ofstream log_file_;
};

#endif


