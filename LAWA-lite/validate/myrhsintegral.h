#ifndef VALIDATE_MYRHSINTEGRAL_H
#define VALIDATE_MYRHSINTEGRAL_H 1

#include <lawa/lawa.h>
#include <validate/mybasis.h>

namespace lawa {

template <typename T>
struct MyRhsIntegral
{
    typedef typename MyBasis<T>::PrimalBasis  PrimalBasis;
    
    MyRhsIntegral(const Function<T> &f, const MyBasis<T> &V);

    T
    operator()(int j, int k, int deriv) const;

    T
    operator()(int absoluteIndex, int deriv) const;
    
    MyBasis<T>      V;
    Function<T>     f;
};

} // namespace lawa

#endif // VALIDATE_MYRHSINTEGRAL_H