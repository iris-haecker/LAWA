#ifndef VALIDATE_MYBASIS_H
#define VALIDATE_MYBASIS_H 1

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

    Support<T>
    support(int absoluteIndex) const;
    
    int
    minK(int j, const T &x) const;

    int
    maxK(int j, const T &x) const;

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

    int
    d() const;

    PrimalBasis  basisLeft, basisRight;
};

} // namespace lawa

#endif // VALIDATE_MYBASIS_H