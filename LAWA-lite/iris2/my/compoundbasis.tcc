#ifndef IRIS2_MY_COMPOUNDBASIS_TCC
#define IRIS2_MY_COMPOUNDBASIS_TCC 1

namespace lawa {

template <typename T>
CompoundBasis<T>::CompoundBasis(int d, int d_, int j0)
    : basisLeft(d, d_, j0), basisRight(d, d_, j0)
{
}

template <typename T>
int
CompoundBasis<T>::minK(XType e, int j, const T &x) const
{
    if (x<T(0.5)) {
        return basisLeft.generator(e).minK(j, x);
    }
    return basisRight.generator(e).minK(j, x);
}

template <typename T>
int
CompoundBasis<T>::maxK(XType e, int j, const T &x) const
{
    if (x<T(0.5)) {
        return basisLeft.generator(e).maxK(j, x);
    }
    return basisRight.generator(e).maxK(j, x);
}

template <typename T>
Support<T>
CompoundBasis<T>::support(XType e, int j, long k) const
{
    return getBasis(e, j, k).generator(e).support(j,k);
}

template <typename T>
const typename CompoundBasis<T>::PrimalBasis &
CompoundBasis<T>::getBasis(XType e, int j, int k) const
{
    int firstIndex = (e==XBSpline) ? basisLeft.mra.rangeI(j).firstIndex()
                                   : basisLeft.rangeJ(j).firstIndex();

    int lastIndex  = (e==XBSpline) ? basisRight.mra.rangeI(j).lastIndex();
                                   : basisRight.rangeJ(j).lastIndex();

    int middle = (firstIndex + lastIndex)/2;

    if (k<=middle) {
        return basisLeft;
    }
    return basisRight;
}

template <typename T>
T
CompoundBasis<T>::operator()(XType e, int j, int k, T x, int deriv) const
{
    if (e==XBSpline) {
        return getBasis(e, j, k).mra.phi(x, j, k, deriv);
    }
    return getBasis(e, j, k).psi(x, j, k, deriv);
}

template <typename T>
int
CompoundBasis<T>::j0() const
{
    return basisLeft.j0;
}

template <typename T>
int
CompoundBasis<T>::d() const
{
    return basisLeft.d;
}

} // namespace lawa

#endif // IRIS2_MY_COMPOUNDBASIS_TCC