#ifndef IRIS2_MY_RHSINTEGRAL2D_H
#define IRIS2_MY_RHSINTEGRAL2D_H 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T>
struct RHSIntegral2D
{

    RHSIntegral2D(const CompoundBasis<T> &V1,
                  const CompoundBasis<T> &V2,
                  Function2D<T> f, int order);

    T
    operator()(XType e1, int j1, int k1,
               XType e2, int j2, int k2) const;

    T
    operator()(const Index1D &lambda1, const Index1D &lambda2) const;

    T
    operator()(const Index2D &lambda) const;

    typedef typename CompoundBasis<T>::PrimalBasis  _Basis;

    const CompoundBasis<T>      &V1;
    const CompoundBasis<T>      &V2;
    Function2D<T>               f;
    Integral2D<FullGridGL, _Basis, _Basis>  leftBottomIntegral2D;   // x=0, y=0
    Integral2D<FullGridGL, _Basis, _Basis>  leftTopIntegral2D;      // x=0, y=1
    Integral2D<FullGridGL, _Basis, _Basis>  rightBottomIntegral2D;  // x=1, y=0
    Integral2D<FullGridGL, _Basis, _Basis>  rightTopIntegral2D;     // x=1, y=1
};

} // namespace lawa

#endif // IRIS2_MY_RHSINTEGRAL2D_H
