#ifndef IRIS2_MY_RHSINTEGRAL1D_TCC
#define IRIS2_MY_RHSINTEGRAL1D_TCC 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T>
RHSIntegral1D<T>::RHSIntegral1D(const CompoundBasis<T> &_V,
                                Function<T> _f,
                                int order)
    : V(_V), f(_f),
      leftIntegralF(f, V.basisLeft),
      rightIntegralF(f, V.basisRight)
{
    leftIntegralF.quadrature.setOrder(order);
    rightIntegralF.quadrature.setOrder(order);
}

template <typename T>
T
RHSIntegral1D<T>::operator()(XType e, int j, int k) const
{
    if (V.isLeft(j,k,e)) {
        return leftIntegralF(j,k,e,0);
    }
    return rightIntegralF(j,k,e,0);
}

template <typename T>
T
RHSIntegral1D<T>::operator()(const Index1D &lambda) const
{
    return operator()(lambda.xtype, lambda.j, lambda.k);
}

} // namespace lawa

#endif // IRIS2_MY_RHSINTEGRAL1D_TCC
