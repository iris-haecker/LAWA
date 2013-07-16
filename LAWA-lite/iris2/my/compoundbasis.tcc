#ifndef IRIS2_MY_COMPOUNDBASIS_TCC
#define IRIS2_MY_COMPOUNDBASIS_TCC 1

namespace lawa {

template <typename T>
CompoundBasis<T>::CompoundBasis(int _d, int _d_, int j0)
    : d(_d), d_(_d_), basisLeft(d, d_, j0), basisRight(d, d_, j0)
{
}

template <typename T>
template <BoundaryCondition LeftBC, BoundaryCondition RightBC>
void
CompoundBasis<T>::enforceBoundaryCondition()
{
    if (LeftBC==DirichletBC) {
        basisLeft.enforceBoundaryCondition<LeftBC>();
    }

    if (RightBC==DirichletBC) {
        basisRight.enforceBoundaryCondition<RightBC>();
    }
}

template <typename T>
int
CompoundBasis<T>::minK(int j, XType e, const T &x) const
{
    if (x<T(0.5)) {
        return basisLeft.generator(e).minK(j, x);
    }
    return basisRight.generator(e).minK(j, x);
}

template <typename T>
int
CompoundBasis<T>::maxK(int j, XType e, const T &x) const
{
    if (x<T(0.5)) {
        return basisLeft.generator(e).maxK(j, x);
    }
    return basisRight.generator(e).maxK(j, x);
}

template <typename T>
Support<T>
CompoundBasis<T>::support(int j, long k, XType e) const
{
    return getBasis(j, k, e).generator(e).support(j,k);
}

template <typename T>
const typename CompoundBasis<T>::PrimalBasis &
CompoundBasis<T>::getBasis(int j, int k, XType e) const
{
    int firstIndex = (e==XBSpline) ? basisLeft.mra.rangeI(j).firstIndex()
                                   : basisLeft.rangeJ(j).firstIndex();

    int lastIndex  = (e==XBSpline) ? basisRight.mra.rangeI(j).lastIndex()
                                   : basisRight.rangeJ(j).lastIndex();

    int middle = (firstIndex + lastIndex)/2;

    if (k<=middle) {
        return basisLeft;
    }
    return basisRight;
}

template <typename T>
bool
CompoundBasis<T>::isLeft(int j, int k, XType e) const
{
    int firstIndex = (e==XBSpline) ? basisLeft.mra.rangeI(j).firstIndex()
                                   : basisLeft.rangeJ(j).firstIndex();

    int lastIndex  = (e==XBSpline) ? basisRight.mra.rangeI(j).lastIndex()
                                   : basisRight.rangeJ(j).lastIndex();

    int middle = (firstIndex + lastIndex)/2;

    return (k<=middle);
}

template <typename T>
bool
CompoundBasis<T>::isRight(int j, int k, XType e) const
{
    int firstIndex = (e==XBSpline) ? basisLeft.mra.rangeI(j).firstIndex()
                                   : basisLeft.rangeJ(j).firstIndex();

    int lastIndex  = (e==XBSpline) ? basisRight.mra.rangeI(j).lastIndex()
                                   : basisRight.rangeJ(j).lastIndex();

    int middle = (firstIndex + lastIndex)/2;

    return (k>middle);
}

//
//  Range of scaling functions
//
template <typename T>
Range<int>
CompoundBasis<T>::rangeI(int j) const
{
    return Range<int>(basisLeft.mra.rangeIL(j).firstIndex(),
                      basisRight.mra.rangeIR(j).lastIndex());
}

template <typename T>
Range<int>
CompoundBasis<T>::rangeIL(int j) const
{
    return basisLeft.mra.rangeIL(j);
}

template <typename T>
Range<int>
CompoundBasis<T>::rangeII(int j) const
{
    return basisRight.mra.rangeII(j);
}

template <typename T>
Range<int>
CompoundBasis<T>::rangeIR(int j) const
{
    return basisRight.mra.rangeIR(j);
}

//
//  Range of wavelet functions
//
template <typename T>
Range<int>
CompoundBasis<T>::rangeJ(int j) const
{
    return Range<int>(basisLeft.rangeJL(j).firstIndex(),
                      basisLeft.rangeJR(j).lastIndex());
}

template <typename T>
Range<int>
CompoundBasis<T>::rangeJL(int j) const
{
    return basisLeft.rangeJL(j);
}

template <typename T>
Range<int>
CompoundBasis<T>::rangeJI(int j) const
{
    return basisLeft.rangeJI(j);
}

template <typename T>
Range<int>
CompoundBasis<T>::rangeJR(int j) const
{
    return basisRight.rangeJR(j);
}

template <typename T>
T
CompoundBasis<T>::operator()(XType e, int j, int k, T x, int deriv) const
{
    if (e==XBSpline) {
        return getBasis(j, k, e).mra.phi(x, j, k, deriv);
    }
    return getBasis(j, k, e).psi(x, j, k, deriv);
}

template <typename T>
int
CompoundBasis<T>::j0() const
{
    return basisLeft.j0;
}

} // namespace lawa

#endif // IRIS2_MY_COMPOUNDBASIS_TCC