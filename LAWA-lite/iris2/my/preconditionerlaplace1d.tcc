#ifndef IRIS2_MY_PRECONDITIONERLAPLACE1D_TCC
#define IRIS2_MY_PRECONDITIONERLAPLACE1D_TCC 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

//-- PreconditionerLaplace1D ---------------------------------------------------

template <typename T>
PreconditionerLaplace1D<T>::PreconditionerLaplace1D(const Laplace1D<T> &_A)
    : A(_A)
{
}

template <typename T>
T
PreconditionerLaplace1D<T>::operator()(XType e, int j, int k) const
{
    return operator()(Index1D(j,k,e));
}

template <typename T>
T
PreconditionerLaplace1D<T>::operator()(const Index1D &index) const
{
    return T(1) / A.pH1(index,index);
}

//-- PreconditionerLaplace1D_M2 ------------------------------------------------

template <typename T>
PreconditionerLaplace1D_M2<T>::PreconditionerLaplace1D_M2(const Laplace1D<T> &_A)
    : A(_A)
{
}

template <typename T>
T
PreconditionerLaplace1D_M2<T>::operator()(XType e1, int j1, int k1,
                                          XType e2, int j2, int k2) const
{
    return operator()(Index1D(j1,k1,e1), Index1D(j2,k2,e2));
}

template <typename T>
T
PreconditionerLaplace1D_M2<T>::operator()(const Index1D &index1,
                                          const Index1D &index2) const
{
    return T(1) / (A.pH1(index1,index1)*A.pH1(index2,index2));
}

template <typename T>
T
PreconditionerLaplace1D_M2<T>::operator()(const Index2D &lambda) const
{
    return operator()(lambda.index1, lambda.index2);
}

} // namespace lawa

#endif // IRIS2_MY_PRECONDITIONERLAPLACE1D_TCC
