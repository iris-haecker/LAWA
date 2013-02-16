#ifndef IRIS_MYINTEGRAL_H
#define IRIS_MYINTEGRAL_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename T>
struct MyIntegral
{
    typedef typename MyBasis<T>::PrimalBasis  PrimalBasis;
    
    MyIntegral(const MyBasis<T> &U, const MyBasis<T> &V);

    T
    operator()(int j1, int k1, int deriv1,
               int j2, int k2, int deriv2) const;

    T
    operator()(int absoluteIndex1, int deriv1,
               int absoluteIndex2, int deriv2) const;
    
    MyBasis<T>  U, V;
};

} // namespace lawa

#endif // IRIS_MYINTEGRAL_H