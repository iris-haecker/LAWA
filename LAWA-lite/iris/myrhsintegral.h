#ifndef IRIS_MYRHSINTEGRAL_H
#define IRIS_MYRHSINTEGRAL_H 1

#include <lawa/lawa.h>
#include <iris/mybasis.h>

namespace lawa {

template <typename T>
struct MyRhsIntegral
{
    typedef typename MyBasis<T>::PrimalBasis  PrimalBasis;
    
    MyRhsIntegral(const Function<T> &f, const MyBasis<T> &V);

    T
    operator()(int j, int k, int deriv) const;

    T
    operator()(int absoluteIndex, int deriv = 0) const;
    
    MyBasis<T>      V;
    Function<T>     f;
};

} // namespace lawa

#endif // IRIS_MYRHSINTEGRAL_H