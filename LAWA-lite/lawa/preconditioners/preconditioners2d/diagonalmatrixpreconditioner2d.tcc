#include <cmath>

namespace lawa {

template <typename T, typename Basis2D, typename BilinearForm>
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::DiagonalMatrixPreconditioner2D(const BilinearForm &a)
    : _a(a)
{
}

template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator()(XType xtype1, int j1, int k1,
                                                                   XType xtype2, int j2, int k2) const
{
    return 1./std::sqrt(fabs(_a(xtype1,j1,k1,xtype2,j2,k2,
                                xtype1,j1,k1,xtype2,j2,k2)));
}

template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}

}   // namespace lawa

