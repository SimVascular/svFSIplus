
#include "fils_struct.hpp"

namespace norm {

using namespace fsi_linear_solver;

double fsi_ls_norms(const int nNo, FSILS_commuType& commu, const Vector<double>& U);

double fsi_ls_normv(const int dof, const int nNo, FSILS_commuType& commu, const Array<double>& U);

};
