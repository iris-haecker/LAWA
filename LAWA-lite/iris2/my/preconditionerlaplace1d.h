#ifndef IRIS2_MY_PRECONDITIONERLAPLACE1D_H
#define IRIS2_MY_PRECONDITIONERLAPLACE1D_H 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T>
struct PreconditionerLaplace1D
{

    PreconditionerLaplace1D(const Laplace1D<T> &A);

    T
    operator()(XType e, int j, int k) const;

    T
    operator()(const Index1D &index) const;

    const Laplace1D<T>    &A;
};

template <typename T>
struct PreconditionerLaplace1D_M2
{

    PreconditionerLaplace1D_M2(const Laplace1D<T> &A);

    T
    operator()(XType xtype1, int j1, int k1,
               XType xtype2, int j2, int k2) const;

    T
    operator()(const Index1D &index1, const Index1D &index2) const;

    T
    operator()(const Index2D &index) const;

    const Laplace1D<T>    &A;
};


} // namespace lawa

#endif // IRIS2_MY_PRECONDITIONERLAPLACE1D_H
