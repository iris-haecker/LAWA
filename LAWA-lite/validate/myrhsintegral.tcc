#ifndef VALIDATE_MYRHSINTEGRAL_TCC
#define VALIDATE_MYRHSINTEGRAL_TCC 1

#include <lawa/lawa.h>

namespace lawa {

template <typename T>
MyRhsIntegral<T>::MyRhsIntegral(const Function<T> &_f, const MyBasis<T> &_V)
    : f(_f), V(_V)
{
}

template <typename T>
T
MyRhsIntegral<T>::operator()(int j, int k, int deriv) const
{
    const PrimalBasis &_V = V.getBasis(j, k);

    XType e = V.getXType(j);

    j = V.getActualLevel(j);

    IntegralF<Gauss, PrimalBasis>   rhsIntegral(f, _V);

    return rhsIntegral(j, k, e, deriv);
}

template <typename T>
T
MyRhsIntegral<T>::operator()(int absoluteIndex, int deriv) const
{
    const int j = V.getLevel(absoluteIndex);
    const int k = V.getIndex(absoluteIndex);

    return operator()(j, k, deriv);
}

} // namespace lawa

#endif // VALIDATE_MYRHSINTEGRAL_TCC