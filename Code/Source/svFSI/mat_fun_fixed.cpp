
#include "mat_fun_fixed.h"

namespace mat_fun_fixed {

Array<int> t_ind;

//----------
// ten_init
//----------
// Initialize tensor index pointer
//
void ten_init(const int nd)
{
  if (t_ind.size() != 0) {
    return;
  }

  int nn = pow(nd, 4);
  t_ind.resize(4, nn);

  int ii = 0;
  for (int l = 0; l < nd; l++) {
    for (int k = 0; k < nd; k++) {
      for (int j = 0; j < nd; j++) {
        for (int i = 0; i < nd; i++) {
          t_ind(0,ii) = i;
          t_ind(1,ii) = j;
          t_ind(2,ii) = k;
          t_ind(3,ii) = l;
          ii = ii + 1;
        }
      }
    }
  }
}

};
