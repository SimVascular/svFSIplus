#ifndef MAT_FUN_H
#define MAT_FUN_H

#include "Array.h"
#include "Tensor4.h"
#include "Vector.h"
#include "mat_fun_fixed.h"

/// @brief The classes defined here duplicate the data structures in the
/// Fortran MATFUN module defined in MATFUN.f.
///
/// This module defines data structures for generally performed matrix and
/// tensor operations.
///
/// \todo [TODO:DaveP] this should just be a namespace?
//
namespace mat_fun {
double mat_ddot(const Array<double>& A, const Array<double>& B, const int nd);
double mat_det(const Array<double>& A, const int nd);
Array<double> mat_dev(const Array<double>& A, const int nd);

Array<double> mat_dyad_prod(const Vector<double>& u, const Vector<double>& v,
                            const int nd);

Array<double> mat_id(const int nsd);
Array<double> mat_inv(const Array<double>& A, const int nd, bool debug = false);
Array<double> mat_inv_ge(const Array<double>& A, const int nd,
                         bool debug = false);
Array<double> mat_inv_lp(const Array<double>& A, const int nd);

Vector<double> mat_mul(const Array<double>& A, const Vector<double>& v);
Array<double> mat_mul(const Array<double>& A, const Array<double>& B);
void mat_mul(const Array<double>& A, const Array<double>& B,
             Array<double>& result);
Array<double> mat_symm(const Array<double>& A, const int nd);
Array<double> mat_symm_prod(const Vector<double>& u, const Vector<double>& v,
                            const int nd);
double mat_trace(const Array<double>& A, const int nd);

Tensor4<double> ten_asym_prod12(const Array<double>& A, const Array<double>& B,
                                const int nd);
Tensor4<double> ten_ddot(const Tensor4<double>& A, const Tensor4<double>& B,
                         const int nd);
Tensor4<double> ten_ddot_2412(const Tensor4<double>& A,
                              const Tensor4<double>& B, const int nd);
Tensor4<double> ten_ddot_3424(const Tensor4<double>& A,
                              const Tensor4<double>& B, const int nd);

Tensor4<double> ten_dyad_prod(const Array<double>& A, const Array<double>& B,
                              const int nd);
Tensor4<double> ten_ids(const int nd);
Array<double> ten_mddot(const Tensor4<double>& A, const Array<double>& B,
                        const int nd);

Tensor4<double> ten_symm_prod(const Array<double>& A, const Array<double>& B,
                              const int nd);
Tensor4<double> ten_transpose(const Tensor4<double>& A, const int nd);

Array<double> transpose(const Array<double>& A);

void ten_init(const int nd);

};  // namespace mat_fun

#endif
