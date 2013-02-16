#ifndef IRIS_MYINTEGRAL_TCC
#define IRIS_MYINTEGRAL_TCC 1

#include <lawa/lawa.h>

namespace lawa {

template <typename T>
MyIntegral<T>::MyIntegral(const MyBasis<T> &_U, const MyBasis<T> &_V)
    : U(_U), V(_V)
{
}

template <typename T>
T
MyIntegral<T>::operator()(int j1, int k1, int deriv1,
                          int j2, int k2, int deriv2) const
{
    const PrimalBasis &_U = U.getBasis(j1, k1);
    const PrimalBasis &_V = V.getBasis(j2, k2);
    
    XType e1 = U.getXType(j1);
    XType e2 = V.getXType(j2);

    j1 = U.getActualLevel(j1);
    j2 = V.getActualLevel(j2);

    Integral<Gauss,PrimalBasis,PrimalBasis>  integral(_U, _V);

    return integral(j1, k1, e1, deriv1, j2, k2, e2, deriv2);
}

template <typename T>
T
MyIntegral<T>::operator()(int absoluteIndex1, int deriv1,
                          int absoluteIndex2, int deriv2) const
{
    const int j1 = U.getLevel(absoluteIndex1);
    const int k1 = U.getIndex(absoluteIndex1);
    
    const int j2 = V.getLevel(absoluteIndex2);
    const int k2 = V.getIndex(absoluteIndex2);

    return operator()(j1, k1, deriv1, j2, k2, deriv2);
}

} // namespace lawa

#endif // IRIS_MYINTEGRAL_TCC