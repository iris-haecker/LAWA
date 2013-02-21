#ifndef IRIS_RHS_H
#define IRIS_RHS_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <vector>

namespace lawa {

template <typename T, typename Precond>
struct MyRhs
{
    MyRhs(const Function<T>    &f,
          const MyOperator<T>  &myOperator,
          const Precond        &P);

    template <typename VX>
        int
        filter(const T &tol, DenseVector<VX> &x) const;

    int                     j0, j1;
    DenseVector<Array<T> >  rhsData;
    T                       norm;
};

} // namespace lawa

#endif // IRIS_RHS_H