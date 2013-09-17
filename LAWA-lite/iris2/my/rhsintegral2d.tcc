#ifndef IRIS2_MY_RHSINTEGRAL2D_TCC
#define IRIS2_MY_RHSINTEGRAL2D_TCC 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T>
RHSIntegral2D<T>::RHSIntegral2D(const CompoundBasis<T> &_V1,
                                const CompoundBasis<T> &_V2,
                                Function2D<T> _f,
                                int order)
    : V1(_V1), V2(_V2), f(_f),
      leftBottomIntegral2D(f, V1.basisLeft, V2.basisLeft),      // x=0, y=0
      leftTopIntegral2D(f, V1.basisLeft, V2.basisRight),        // x=0, y=1
      rightBottomIntegral2D(f, V1.basisRight, V2.basisLeft),    // x=1, y=0
      rightTopIntegral2D(f, V1.basisRight, V2.basisRight)       // x=1, y=1
{
    leftBottomIntegral2D.quadrature.setOrder(order);
    leftTopIntegral2D.quadrature.setOrder(order);
    rightBottomIntegral2D.quadrature.setOrder(order);
    rightTopIntegral2D.quadrature.setOrder(order);
}

template <typename T>
T
RHSIntegral2D<T>::operator()(XType e1, int j1, int k1,
                             XType e2, int j2, int k2) const
{
    if (V1.isLeft(j1,k1,e1)) {  // x = 0
        if (V2.isLeft(j2,k2,e2)) { // y = 0
            return leftBottomIntegral2D(j1,k1,e1,0,j2,k2,e2,0);
        } else {                   // y = 1
            return leftTopIntegral2D(j1,k1,e1,0,j2,k2,e2,0);
        }
    } else {                    // x = 1
        if (V2.isLeft(j2,k2,e2)) { // y = 0
            return rightBottomIntegral2D(j1,k1,e1,0,j2,k2,e2,0);
        } else {                   // y = 1
            return rightTopIntegral2D(j1,k1,e1,0,j2,k2,e2,0);
        }
    }
}

template <typename T>
T
RHSIntegral2D<T>::operator()(const Index1D &lambda1,
                             const Index1D &lambda2) const
{
    return operator()(lambda1.xtype, lambda1.j, lambda1.k,
                      lambda2.xtype, lambda2.j, lambda2.k);
}

template <typename T>
T
RHSIntegral2D<T>::operator()(const Index2D &lambda) const
{
    return operator()(lambda.index1, lambda.index2);
}

} // namespace lawa

#endif // IRIS2_MY_RHSINTEGRAL2D_TCC
