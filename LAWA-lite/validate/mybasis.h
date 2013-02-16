#ifndef IRIS_MYBASIS_H
#define IRIS_MYBASIS_H 1

#include <lawa/lawa.h>

namespace lawa {

template <typename T>
struct MyBasis
{
    typedef Basis<T, Primal, Interval, Dijkema>  PrimalBasis;

    MyBasis(int d, int d_, int j0=-1);

    int
    firstIndex(int j) const;

    int
    lastIndex(int j) const;

    int
    length(int j) const;
    
    int
    getAbsoluteIndex(int j, int k)const;
    
    int
    getFirstAbsoluteIndex(int j) const;

    int
    getLastAbsoluteIndex(int j) const;
    
    int
    getLevel(int absoluteIndex) const;

    int
    getIndex(int absoluteIndex) const;

    const PrimalBasis &
    getBasis(int j, int k) const;

    PrimalBasis &
    getBasis(int j, int k);

    XType
    getXType(int j) const;
    
    int
    getActualLevel(int j) const;

    T
    operator()(T x, int j, int k, int deriv) const;

    T
    operator()(T x, int absoluteIndex, int deriv) const;

    int
    j0() const;

    PrimalBasis  basisLeft, basisRight;
};

} // namespace lawa

#endif // IRIS_MYBASIS_H