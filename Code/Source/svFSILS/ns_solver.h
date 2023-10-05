
#include "fils_struct.hpp"

namespace ns_solver {

void bc_pre(fsi_linear_solver::FSILS_lhsType& lhs, const int nsd, const int dof,
            const int nNo, const int mynNo);

void depart(fsi_linear_solver::FSILS_lhsType& lhs, const int nsd, const int dof,
            const int nNo, const int nnz, const Array<double>& Val,
            Array<double>& Gt, Array<double>& mK, Array<double>& mG,
            Array<double>& mD, Vector<double>& mL);

void ns_solver(fsi_linear_solver::FSILS_lhsType& lhs,
               fsi_linear_solver::FSILS_lsType& ls, const int dof,
               const Array<double>& Val, Array<double>& Ri);

};  // namespace ns_solver
