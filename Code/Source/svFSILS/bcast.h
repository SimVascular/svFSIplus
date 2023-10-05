
#include "fils_struct.hpp"

namespace bcast {

using namespace fsi_linear_solver;

void fsils_bcast(double& u, FSILS_commuType& commu);

void fsils_bcast_v(const int n, Vector<double>& u, FSILS_commuType& commu);

};  // namespace bcast
