#include <cmath>

namespace lawa {

template <typename T, typename Basis>
H1NormPreconditioner2D<T,Basis>::H1NormPreconditioner2D(const Basis &basis)
    : _integral_x(basis.first, basis.first), _integral_y(basis.second, basis.second)
{
}

template <typename T, typename Basis>
T
H1NormPreconditioner2D<T,Basis>::operator()(XType xtype1, int j1, int k1,
                                            XType xtype2, int j2, int k2) const
{
    T dd_x = _integral_x(j1,k1,xtype1,1,j1,k1,xtype1,1);
    T id_x = _integral_x(j1,k1,xtype1,0,j1,k1,xtype1,0);
    T dd_y = _integral_y(j2,k2,xtype2,1,j2,k2,xtype2,1);
    T id_y = _integral_y(j2,k2,xtype2,0,j2,k2,xtype2,0);
    return 1./std::sqrt(dd_x*id_y + id_x*dd_y + id_x*id_y);
}

template <typename T, typename Basis>
T
H1NormPreconditioner2D<T,Basis>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype,index.index1.j,index.index1.k,
                            index.index2.xtype,index.index2.j,index.index2.k);
}

}   // namespace lawa

