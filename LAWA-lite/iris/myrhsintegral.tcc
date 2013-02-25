#ifndef IRIS_MYRHSINTEGRAL_TCC
#define IRIS_MYRHSINTEGRAL_TCC 1

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

    const T value = rhsIntegral(j, k, e, deriv);

/*
    std::cerr << "(j, k, e) = (" << j << ", " << k << ", " << e << ") : "
              << value << std::endl;
*/
    return value;
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

#endif // IRIS_MYRHSINTEGRAL_TCC