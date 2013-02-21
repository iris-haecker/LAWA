#ifndef IRIS_MYBASIS_TCC
#define IRIS_MYBASIS_TCC 1

namespace lawa {

template <typename T>
MyBasis<T>::MyBasis(int d, int d_, int j0)
    : basisLeft(d, d_, j0), basisRight(d, d_, j0)
{
    std::cerr << "j0 = " << basisLeft.j0 << std::endl;
}
    
template <typename T>
int
MyBasis<T>::firstIndex(int j) const
{
    if (j<basisLeft.j0) {
        return basisLeft.mra.rangeI(j+1).firstIndex();
    }
    return basisLeft.rangeJ(j).firstIndex();
}

template <typename T>
int
MyBasis<T>::lastIndex(int j) const
{
    if (j<basisLeft.j0) {
        return basisRight.mra.rangeI(j+1).lastIndex();
    }
    return basisRight.rangeJ(j).lastIndex();
}

template <typename T>
int
MyBasis<T>::length(int j) const
{
    return lastIndex(j) - firstIndex(j) + 1;
}
    
template <typename T>
int
MyBasis<T>::getAbsoluteIndex(int j, int k) const
{
    return getFirstAbsoluteIndex(j) + k - firstIndex(j);
}
    
template <typename T>
int
MyBasis<T>::getFirstAbsoluteIndex(int j) const
{
    int result = 1;
    for (int l=basisLeft.j0-1; l<j; ++l) {
        result += length(l);
    }
    return result;
}

template <typename T>
int
MyBasis<T>::getLastAbsoluteIndex(int j) const
{
    return getFirstAbsoluteIndex(j) + length(j) - 1;
}

template <typename T>
int
MyBasis<T>::getLevel(int absoluteIndex) const
{
    int j = j0()-1;

    while (absoluteIndex>0) {
        absoluteIndex -= length(j);
        if (absoluteIndex>0) {
            ++j;
        }
    }
    return j;
}

template <typename T>
int
MyBasis<T>::getIndex(int absoluteIndex) const
{
    int j = j0()-1;
    int k = absoluteIndex;

    while (absoluteIndex>0) {
        absoluteIndex -= length(j);
        if (absoluteIndex>0) {
            ++j;
            k = absoluteIndex;
        }
    }
    return firstIndex(j) + k - 1;
}

template <typename T>
Support<T>
MyBasis<T>::support(int absoluteIndex) const
{
    int  j = getLevel(absoluteIndex);
    int  k = getIndex(absoluteIndex);
    
    int    _j    = getActualLevel(j);
    XType  xType = getXType(j);
    
    return getBasis(j, k).generator(xType).support(_j,k);
}

template <typename T>
int
MyBasis<T>::minK(int j, const T &x) const
{
    int    _j    = getActualLevel(j);
    XType  xType = getXType(j);

    if (x<T(0.5)) {
        return basisLeft.generator(xType).minK(_j, x);
    }
    return basisRight.generator(xType).minK(_j, x);
}

template <typename T>
int
MyBasis<T>::maxK(int j, const T &x) const
{
    int    _j    = getActualLevel(j);
    XType  xType = getXType(j);

    if (x<T(0.5)) {
        return basisLeft.generator(xType).maxK(_j, x);
    }
    return basisRight.generator(xType).maxK(_j, x);
}

template <typename T>
typename MyBasis<T>::PrimalBasis &
MyBasis<T>::getBasis(int j, int k)
{
    int middle = (firstIndex(j) + lastIndex(j))/2;
    if (k<=middle) {
        return basisLeft;
    }
    return basisRight;
}

template <typename T>
const typename MyBasis<T>::PrimalBasis &
MyBasis<T>::getBasis(int j, int k) const
{
    int middle = (firstIndex(j) + lastIndex(j))/2;
    if (k<=middle) {
        return basisLeft;
    }
    return basisRight;
}

template <typename T>
XType
MyBasis<T>::getXType(int j) const
{
    if (j<basisLeft.j0) {
        return XBSpline;        
    }
    return XWavelet;
}

template <typename T>
int
MyBasis<T>::getActualLevel(int j) const
{
    return j + 1 - getXType(j);
}

template <typename T>
T
MyBasis<T>::operator()(T x, int j, int k, int deriv) const
{
    if (getXType(j)==XBSpline) {
        return getBasis(j+1,k).mra.phi(x, j+1, k, deriv);
    }
    return getBasis(j,k).psi(x, j, k, deriv);
}

template <typename T>
T
MyBasis<T>::operator()(T x, int absoluteIndex, int deriv) const
{
    int j = getLevel(absoluteIndex);
    int k = getIndex(absoluteIndex);

    return operator()(x, j, k, deriv);
}

template <typename T>
int
MyBasis<T>::j0() const
{
    return basisLeft.j0;
}

template <typename T>
int
MyBasis<T>::d() const
{
    return basisLeft.d;
}


} // namespace lawa

#endif // IRIS_MYBASIS_TCC