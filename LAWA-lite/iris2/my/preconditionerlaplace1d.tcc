#ifndef IRIS2_MY_PRECONDITIONERLAPLACE1D_TCC
#define IRIS2_MY_PRECONDITIONERLAPLACE1D_TCC 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

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
    return T(1) / A(index,index);
    //return T(1) / A.pH1(index,index);
}

} // namespace lawa

#endif // IRIS2_MY_PRECONDITIONERLAPLACE1D_TCC
