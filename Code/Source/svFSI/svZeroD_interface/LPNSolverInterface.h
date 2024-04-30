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

#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#ifndef LPNSolverInterface_h
#define LPNSolverInterface_h

//--------------------
// LPNSolverInterface
//--------------------
//
class LPNSolverInterface
{
  public:
    LPNSolverInterface();
    ~LPNSolverInterface();

    void load_library(const std::string& interface_lib);
    void initialize(std::string file_name);
    void increment_time(const double time, std::vector<double>& solution);
    void run_simulation(const double time, std::vector<double>& output_times, 
	                std::vector<double>& output_solutions, int& error_code);
    void update_block_params(std::string block_name, std::vector<double>& new_params);
    void read_block_params(std::string block_name, std::vector<double>& block_params);
    void get_block_node_IDs(std::string block_name, std::vector<int>& IDs);
    void update_state(std::vector<double> state_y, std::vector<double> state_ydot);
    void return_y(std::vector<double>& y);
    void return_ydot(std::vector<double>& ydot);
    void set_external_step_size(double step_size);

    // Interface functions.
    std::string lpn_initialize_name_;
    void (*lpn_initialize_)(std::string, int&, int&, int&, int&, std::vector<std::string>&, 
	                    std::vector<std::string>&);

    std::string lpn_increment_time_name_;
    void (*lpn_increment_time_)(const int, const double, std::vector<double>& solution);

    std::string lpn_run_simulation_name_;
    void (*lpn_run_simulation_)(const int, const double, std::vector<double>& output_times, 
	                        std::vector<double>& output_solutions, int& error_code);
    
    std::string lpn_update_block_params_name_;
    void (*lpn_update_block_params_)(const int, std::string, std::vector<double>& new_params);
    
    std::string lpn_read_block_params_name_;
    void (*lpn_read_block_params_)(const int, std::string, std::vector<double>& block_params);
    
    std::string lpn_get_block_node_IDs_name_;
    void (*lpn_get_block_node_IDs_)(const int, std::string, std::vector<int>& block_params);
    
    std::string lpn_update_state_name_;
    void (*lpn_update_state_)(const int, std::vector<double>, std::vector<double>);
    
    std::string lpn_return_y_name_;
    void (*lpn_return_y_)(const int, std::vector<double>&);
    
    std::string lpn_return_ydot_name_;
    void (*lpn_return_ydot_)(const int, std::vector<double>&);
    
    std::string lpn_set_external_step_size_name_;
    void (*lpn_set_external_step_size_)(const int, double);

    void* library_handle_ = nullptr;
    int problem_id_ = 0;
    int system_size_ = 0;
    int num_cycles_ = 0;
    int pts_per_cycle_ = 0;
    int num_output_steps_ = 0;
    std::vector<std::string> block_names_;
    std::vector<std::string> variable_names_;
};

#endif

