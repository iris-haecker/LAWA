#ifndef IRIS_MYRESTRICT_H
#define IRIS_MYRESTRICT_H 1

#include <lawa/flensforlawa.h>
#include <lawa/lawa.h>

namespace lawa {

template <typename Operator, typename Precond, typename TI, typename CRS>
    void
    myRestrict(const Operator       &A,
               const Precond        &P,
               const IndexSet<TI>   &Lambda,
               SparseGeMatrix<CRS>  &B);

template <typename VX, typename TI, typename VY>
    void
    myRestrict(const DenseVector<VX>  &x,
               const IndexSet<TI>     &Lambda,
               DenseVector<VY>        &y);

template <typename VX, typename TI, typename VY>
    void
    myRestrictSub(const DenseVector<VX>  &x,
                  const IndexSet<TI>     &Lambda,
                  DenseVector<VY>        &y);

} // namespace lawa

#endif // IRIS_MYRESTRICT_H