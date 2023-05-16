
#include "utils.h"

#include <bitset>
#include <chrono>
#include <cmath> 
#include <limits>

#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <sys/resource.h>

/* MacOS
#include <mach/task.h>
#include <mach/mach_init.h>
*/

namespace utils {

//-----------
// CountBits
//-----------
// Count the number of bits in an int.
//
int CountBits(int n)
{
  int count = 0;
  while (n) {
    count += n & 1;
    n >>= 1;
  }
  return count;
}

//------
// cput
//------
//
double cput()
{
  auto now = std::chrono::system_clock::now();
  auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);

  auto value = now_ms.time_since_epoch();
  auto duration = value.count() / 1000.0;
  return static_cast<double>(duration);
}

//-------
// cross
//-------
//
Vector<double> 
cross(const Array<double>& V) 
{
  int num_rows = V.nrows();
  Vector<double> U(num_rows);

  if (num_rows == 2) {
    U(0) =  V(1,0);
    U(1) = -V(0,0);
  } else { 
    U(0) = V(1,0)*V(2,1) - V(2,0)*V(1,1);
    U(1) = V(2,0)*V(0,1) - V(0,0)*V(2,1);
    U(2) = V(0,0)*V(1,1) - V(1,0)*V(0,1);
  } 

  return U;
}

//-------
// btest
//-------
//
bool btest(int value, int pos)
{
  return value & (1 << pos); 
}

//---------
// dequeue
//---------
//
bool dequeue(queueType& que, int& iVal) 
{
  int i;
  bool flag;

  if (que.n == 0) { 
    iVal = -1;
    flag = false; 
  } else { 
    iVal = que.v(0);
    for (i = 1; i < que.n; i++) { 
      que.v(i-1) = que.v(i);
    } 
    que.n = que.n - 1;
    flag  = true; 
  }

  return flag;
}

//---------
// enqueue 
//---------
//
void enqueue(queueType& que, int iVal) 
{
  if (que.maxN == 0) {
    que.n = 1;
    que.maxN = 8;
    que.v.resize(que.maxN);
    que.v = -1;
    que.v(0) = iVal;

  } else { 
    if (que.maxN <= que.n) {
      Vector<int> tmp(que.maxN);
      tmp = que.v;
      que.maxN = 4 * que.maxN;
      que.v.resize(que.maxN);
      que.v = -1;

      for (int i = 0; i < que.n; i++) {
        que.v(i) = tmp(i);
      }
    }

    // Check if the new val to be added is already a member of the queue
    bool flag = true;

    for (int i = 0; i < que.n; i++) {
      if (que.v(i) == iVal) {
        flag = false;
        break;
      }
    }

    if (!flag) {
      return;
    }
    que.v(que.n) = iVal;
    que.n = que.n + 1;
  }
}

//-------
// ibclr
//-------
// Clear a bit. 
//
int ibclr(int value, int pos)
{
  return value & ~(1UL << pos);
}

//-------
// ibset 
//-------
// Returns 'value' with the bit at position 'pos' set to one.
//
int ibset(int value, int pos)
{
  return value | 1UL << pos;
}

//---------
// is_zero
//---------
//
bool is_zero(double value1, double value2)
{
  double a = std::abs(value1);
  double b = std::abs(value2);

  if (b > a) { 
    double tmp = a;
    a = b;
    b = tmp;
  } 

  double eps = std::numeric_limits<double>::epsilon();
  double nrm = std::fmax(a,eps);

  if ((a-b)/nrm < 10.0*eps) { 
    return true;
  }

  return false;
}

//-----------
// mem_usage
//-----------
//
double mem_usage(const bool print_usage, const std::string& prefix) 
{
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;
  ret = getrusage(who, &usage);
  double resident_set = usage.ru_maxrss; 

  std::string pr_prefix = "[cpp] ";
  if (prefix != "") {
    pr_prefix += ": " + prefix + " ";
  } 

  if (print_usage) {
    double s2 = (1024.0 * 1024.0);
    int m_rss = static_cast<int>(round(resident_set/s2));
    std::cout << std::endl;
    std::cout << pr_prefix << "#########################################################" << std::endl;
    std::cout << pr_prefix << "###############   m e m o r y   u s a g e  ##############" << std::endl;
    std::cout << pr_prefix << "#########################################################" << std::endl;
    std::cout << pr_prefix << "Current resident set: " << m_rss << std::endl;
    std::cout << pr_prefix << "#########################################################" << std::endl;
    std::cout << std::endl;
  }

  /* MacOS
  task_t task = MACH_PORT_NULL;
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  double vm_usage = 0.0;

  if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) {
    vm_usage = 0.0;
  }
  double s2 = (1024.0 * 1024.0);
  int m_rss = static_cast<int>(round(t_info.resident_size/s2));
  std::cout << "####### Mac Current resident set: " << m_rss << std::endl;
  //vm_usage = t_info.resident_size;
  */

  return resident_set;
}

//------
// norm
//------
// This function will compute second NORM of a vector.
//
// Replicates 'PURE FUNCTION NORMS(U, V)' defined in UTIL.f.
//
double norm(const Vector<double>& U)
{
  double norm = 0.0;
  for (int i = 0; i < U.size(); i++) {
    norm += U(i)*U(i);
  }
  return norm;
}

double norm(const Vector<double>& U, const Vector<double>& V)
{
  double norm = 0.0;
  for (int i = 0; i < U.size(); i++) {
    norm += U(i)*V(i);
  }
  return norm;
}

double norm(const Array<double>& U)
{
  int m = U.nrows();
  int n = U.ncols();
  double norm = 0.0;

  switch (m) { 
    case 1:
      for (int i = 0; i < n; i++) {
        norm = norm + U(0,i)*U(0,i);
      }
    break;

    case 2:
      for (int i = 0; i < n; i++) {
        norm = norm + U(0,i)*U(0,i) + U(1,i)*U(1,i);
      }
    break;

    case 3:
      for (int i = 0; i < n; i++) {
        norm = norm + U(0,i)*U(0,i) + U(1,i)*U(1,i) + U(2,i)*U(2,i);
      }
    break;

    case 4:
      for (int i = 0; i < n; i++) {
        norm = norm + U(0,i)*U(0,i) + U(1,i)*U(1,i) + U(2,i)*U(2,i); + U(3,i)*U(3,i);
      }
    break;

    default: 
      for (int i = 0; i < n; i++) {
        for (int j = 0; i < m; i++) {
          norm = norm + U(j,i)*U(j,i);
        }
      }
    break;
    }

  return norm;
}

//-----------
// print_mem
//-----------
//
void print_mem(const std::string& type, const std::string& prefix, const double memory_in_use, const double memory_returned)
{
  double s = (1024.0);
  double s2 = (1024.0 * 1024.0);
  int muse = static_cast<int>(round(memory_in_use/s));
  int mret = static_cast<int>(round(memory_returned/s));

  double resident_set = mem_usage();
  int m_rss = static_cast<int>(round(resident_set/s2));

  std::string pr_prefix = "[" + type + "] ";
  if (prefix != "") {
    pr_prefix += ": " + prefix + " ";
  } 

  std::cout << std::endl;
  std::cout << pr_prefix << "#########################################################" << std::endl;
  std::cout << pr_prefix << "####################   m e m o r y   ####################" << std::endl;
  std::cout << pr_prefix << "#########################################################" << std::endl;
  std::cout << pr_prefix << "Memory in use: " << muse << std::endl;
  std::cout << pr_prefix << "Memory returned: " << mret << std::endl;
  std::cout << pr_prefix << "Current resident set: " << m_rss << std::endl;
  std::cout << pr_prefix << "#########################################################" << std::endl;
  std::cout << std::endl;
}

//-------------
// print_stats
//-------------
//
void print_stats(const std::string& type, const std::string& prefix, const int allocated, const int active)
{
  std::string pr_prefix = "[" + type + "] ";
  if (prefix != "") {
    pr_prefix += ": " + prefix + " ";
  } 

  std::cout << pr_prefix << "#########################################################" << std::endl;
  std::cout << pr_prefix << "####################   s t a t s   ######################" << std::endl;
  std::cout << pr_prefix << "#########################################################" << std::endl;
  std::cout << pr_prefix << "Number of allocated: " << allocated << std::endl;
  std::cout << pr_prefix << "Number of active: " << active << std::endl;
  std::cout << pr_prefix << "#########################################################" << std::endl;
}

//------------
// pull_stack
//------------
//
bool pull_stack(stackType& stk, int& iVal)
{
  bool flag;

  if (stk.n == 0) {
    iVal = -1;
    flag = false;
  } else { 
    stk.n = stk.n - 1;
    iVal  = stk.v[stk.n];
    flag  = true;
  } 

  return flag;
}

//------------
// push_stack
//------------
// Push a list of values onto the stack.
//
void push_stack(stackType& stk, std::initializer_list<int> values)
{
  for (int value : values) {
    push_stack(stk, value);
  }
}

//------------
// push_stack
//------------
// Push a single value onto the stack.
//
void push_stack(stackType& stk, int iVal)
{
  // This is a new stack
  if (stk.maxN == 0) {
    stk.maxN = 1000;
    stk.v = Vector<int>(stk.maxN);
  }

  stk.v[stk.n] = iVal;
  stk.n += 1;

  if (stk.maxN == stk.n) {
   stk.maxN = 4*stk.maxN;
   stk.v.grow(stk.maxN, -1);
  }
}

//------
// sign
//------
int sign(double value)
{
  int result = -1;

  if (is_zero(value)) { 
    result = 0; 
  } else if (value > 0.0) { 
    result = 1; 
  }

  return result; 
}

//------
// swap
//------
//
void swap(int& value1, int& value2)
{
  std::swap(value1, value2);
}

};
