#ifndef IRIS2_MY_COMPOUNDBASIS_H
#define IRIS2_MY_COMPOUNDBASIS_H 1

#include <lawa/constructions/interval/dijkema/primal/basis.h>

namespace lawa {

template <typename T>
struct CompoundBasis
{
    typedef Basis<T, Primal, Interval, Dijkema>      PrimalBasis;
    typedef typename PrimalBasis::BasisFunctionType  BasisFunctionType;

    CompoundBasis(int d, int d_, int j0=-1);

    template <BoundaryCondition Left, BoundaryCondition Right>
        void
        enforceBoundaryCondition();

    int
    minK(int j, XType e, const T &x) const;

    int
    maxK(int j, XType e, const T &x) const;

    Support<T>
    support(int j, long k, XType e) const;

    DenseVector<Array<T> >
    singularSupport(int j, int k, XType e) const;

    const PrimalBasis &
    getBasis(int j, int k, XType e) const;

    bool
    isLeft(int j, int k, XType e) const;

    bool
    isRight(int j, int k, XType e) const;

//
//  Range of scaling functions
//
    Range<int>
    rangeI(int j) const;

    Range<int>
    rangeIL(int j=0) const;

    Range<int>
    rangeII(int j) const;

    Range<int>
    rangeIR(int j) const;

//
//  Range of wavelet functions
//
    Range<int>
    rangeJ(int j) const;

    Range<int>
    rangeJL(int j=0) const;

    Range<int>
    rangeJI(int j) const;

    Range<int>
    rangeJR(int j) const;

    T
    operator()(XType e, int j, int k, T x, int deriv) const;

    int
    j0() const;

    int          d, d_;
    PrimalBasis  basisLeft, basisRight;
};

} // namespace lawa

#endif // IRIS2_MY_COMPOUNDBASIS_H