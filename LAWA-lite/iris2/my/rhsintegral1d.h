#ifndef IRIS2_MY_RHSINTEGRAL1D_H
#define IRIS2_MY_RHSINTEGRAL1D_H 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T>
struct RHSIntegral1D
{

    RHSIntegral1D(const CompoundBasis<T> &V, Function<T> f, int order);

    T
    operator()(XType e, int j, int k) const;

    T
    operator()(const Index1D &lambda) const;

    typedef typename CompoundBasis<T>::PrimalBasis  _Basis;

    const CompoundBasis<T>      &V;
    Function<T>                 f;
    IntegralF<Gauss, _Basis>    leftIntegralF;
    IntegralF<Gauss, _Basis>    rightIntegralF;

};

} // namespace lawa

#endif // IRIS2_MY_RHSINTEGRAL1D_H
