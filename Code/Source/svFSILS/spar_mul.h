
#include "fils_struct.hpp"

namespace spar_mul {

using namespace fsi_linear_solver;

void fsils_spar_mul_ss(FSILS_lhsType& lhs, const Array<int>& rowPtr,
                       const Vector<int>& colPtr, const Vector<double>& K,
                       const Vector<double>& U, Vector<double>& KU);

void fsils_spar_mul_sv(FSILS_lhsType& lhs, const Array<int>& rowPtr,
                       const Vector<int>& colPtr, const int dof,
                       const Array<double>& K, const Vector<double>& U,
                       Array<double>& KU);

void fsils_spar_mul_vs(FSILS_lhsType& lhs, const Array<int>& rowPtr,
                       const Vector<int>& colPtr, const int dof,
                       const Array<double>& K, const Array<double>& U,
                       Vector<double>& KU);

void fsils_spar_mul_vv(FSILS_lhsType& lhs, const Array<int>& rowPtr,
                       const Vector<int>& colPtr, const int dof,
                       const Array<double>& K, const Array<double>& U,
                       Array<double>& KU);

};  // namespace spar_mul
