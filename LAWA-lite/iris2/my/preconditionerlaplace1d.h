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

} // namespace lawa

#endif // IRIS2_MY_PRECONDITIONERLAPLACE1D_H
