#ifndef IRIS2_MY_COMPOUNDBASIS_H
#define IRIS2_MY_COMPOUNDBASIS_H 1

namespace lawa {

template <typename T>
struct CompoundBasis
{
    typedef Basis<T, Primal, Interval, Dijkema>     PrimalBasis;
    typedef PrimalBasis::BasisFunctionType          BasisFunctionType;

    CompoundBasis(int d, int d_, int j0=-1);

    virtual long
    minK(XType e, int j, const T &x) const;

    virtual long
    maxK(XType e, int j, const T &x) const;

    virtual Support<T>
    support(XType e, int j, long k) const;

    const PrimalBasis &
    getBasis(XType e, int j, int k) const;

    T
    operator()(XType e, int j, int k, T x, int deriv) const;

    int
    j0() const;

    int
    d() const;

    PrimalBasis  basisLeft, basisRight;
};

} // namespace lawa

#endif // IRIS2_MY_COMPOUNDBASIS_H