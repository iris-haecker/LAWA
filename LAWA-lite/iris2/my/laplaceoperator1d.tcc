#ifndef IRIS2_MY_LAPLACEOPERATOR_TCC
#define IRIS2_MY_LAPLACEOPERATOR_TCC 1

#include <cassert>
#include <iostream>

namespace lawa {


template <typename T>
LaplaceOperator1D<T>::LaplaceOperator1D(int d, int d_, int jMin)
    : U(d, d_, jMin), V(d, d_, jMin), j0(std::min(U.j0(), V.j0()))
{
    //
    // Functions from U and V must vanish on the left and right boundary
    //
    V.basisRight.enforceBoundaryCondition<lawa::DirichletBC>();
    V.basisLeft.enforceBoundaryCondition<lawa::DirichletBC>();

    U.basisRight.enforceBoundaryCondition<lawa::DirichletBC>();
    U.basisLeft.enforceBoundaryCondition<lawa::DirichletBC>();
}

template <typename T>
LaplaceOperator1D<T>::LaplaceOperator1D(const MyLaplaceOperator1D<T> &A)
    : U(A.U), V(A.V), j0(A.j0)
{
}


template <typename T>
T
LaplaceOperator1D<T>::operator()(XType e1, int j1, int k1,
                                 XType e2, int j2, int k2) const
{
    typedef typename MyBasis<T>::PrimalBasis  PrimalBasis;

    const PrimalBasis &_U = U.getBasis(e1, j1, k1);
    const PrimalBasis &_V = V.getBasis(e2, j2, k2);

    Integral<Gauss,PrimalBasis,PrimalBasis>  integral(_V, _U);

    // v_x * u_x
    return integral(j1, k1, e1, 1, j2, k2, e2, 1);
}

template <typename T>
T
LaplaceOperator1D<T>::operator()(const Index1D &row_index,
                                 const Index1D &col_index) const
{
    return operator()(row_index.xtype, row_index.j, row_index.k,
                      col_index.xtype, col_index.j, col_index.k);
}

} // namespace lawa

#endif // IRIS2_MY_LAPLACEOPERATOR_TCC
