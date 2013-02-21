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

} // namespace lawa

#endif // IRIS_MYRESTRICT_H