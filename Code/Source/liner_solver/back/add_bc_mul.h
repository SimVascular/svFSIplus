
#include "fils_struct.hpp"

namespace add_bc_mul {

using namespace fsi_linear_solver;

void add_bc_mul(FSILS_lhsType& lhs, const BcopType op_Type, const int dof, const Array<double>& X, Array<double>& Y);

};
