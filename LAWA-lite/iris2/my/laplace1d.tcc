#ifndef IRIS2_MY_LAPLACE1D_TCC
#define IRIS2_MY_LAPLACE1D_TCC 1

#include <cassert>
#include <iostream>

#include <lawa/lawa.h>
#include <iris2/iris2.h>

namespace lawa {

template <typename T>
Laplace1D<T>::Laplace1D(int d, int d_)
    : U(d, d_), V(d, d_), j0(std::min(U.j0(), V.j0()))
{
    using namespace std;
    using namespace lawa;

    //
    // Functions from U and V must vanish on the left and right boundary
    //

    U.basisLeft.enforceBoundaryCondition<DirichletBC>();
    U.basisRight.enforceBoundaryCondition<DirichletBC>();

    V.basisLeft.enforceBoundaryCondition<DirichletBC>();
    V.basisRight.enforceBoundaryCondition<DirichletBC>();
}

template <typename T>
Laplace1D<T>::Laplace1D(int d, int d_, int jMin)
    : U(d, d_, jMin), V(d, d_, jMin), j0(std::min(U.j0(), V.j0()))
{
    //
    // Functions from U and V must vanish on the left and right boundary
    //
    U.basisLeft.enforceBoundaryCondition<DirichletBC>();
    U.basisRight.enforceBoundaryCondition<DirichletBC>();

    V.basisLeft.enforceBoundaryCondition<DirichletBC>();
    V.basisRight.enforceBoundaryCondition<DirichletBC>();
}

template <typename T>
Laplace1D<T>::Laplace1D(const Laplace1D<T> &A)
    : U(A.U), V(A.V), j0(A.j0)
{
}

template <typename T>
int
Laplace1D<T>::minK2(int j1, int k1, XType e1, int j2, XType e2) const
{
    const T x = V.support(j1,k1,e1).l1;

    return U.minK(j2, e2, x);
}

template <typename T>
int
Laplace1D<T>::maxK2(int j1, int k1, XType e1, int j2, XType e2) const
{
    const T x = V.support(j1,k1,e1).l2;

    return U.maxK(j2, e2, x);
}

template <typename T>
T
Laplace1D<T>::operator()(XType e1, int j1, int k1,
                         XType e2, int j2, int k2) const
{
    typedef typename CompoundBasis<T>::PrimalBasis  PrimalBasis;

    const PrimalBasis &_U = U.getBasis(j1, k1, e1);
    const PrimalBasis &_V = V.getBasis(j2, k2, e2);

    Integral<Gauss,PrimalBasis,PrimalBasis>  integral(_V, _U);

    // v_x * u_x
    return integral(j1, k1, e1, 1, j2, k2, e2, 1);
}

template <typename T>
T
Laplace1D<T>::operator()(const Index1D &row_index,
                         const Index1D &col_index) const
{
    return operator()(row_index.xtype, row_index.j, row_index.k,
                      col_index.xtype, col_index.j, col_index.k);
}

template <typename T>
int
Laplace1D<T>::d() const
{
    return U.d;
}


} // namespace lawa

#endif // IRIS2_MY_LAPLACE1D_TCC
