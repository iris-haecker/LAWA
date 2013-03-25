#ifndef IRIS_MYRHS2_H
#define IRIS_MYRHS2_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <vector>

namespace lawa {

template <typename T, typename Precond>
struct MyRhs2
{
    MyRhs2(const Function<T>     &f,
           const MyOperator<T>   &myOperator,
           const Precond         &P,
           T                     eps);

    template <typename VX>
        int
        filter(const T          &tol,
               DenseVector<VX>  &x) const;

    template <typename TI, typename VX>
        int
        filter(const T             &tol,
               const IndexSet<TI>  &Lambda,
               DenseVector<VX>     &x) const;

    int                     j0, j1;
    DenseVector<Array<T> >  rhsData;
    T                       norm;
};

} // namespace lawa

#endif // IRIS_MYRHS2_H