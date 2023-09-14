
#include "fils_struct.hpp"

namespace bicgs {

void bicgsv(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof,
    const Array<double>& K, Array<double>& R);

void bicgss(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const Vector<double>& K, Vector<double>& R);

};
