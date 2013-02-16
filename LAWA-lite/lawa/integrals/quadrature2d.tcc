/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <cmath>

namespace lawa {

template <typename Integral2D>
Quadrature2D<SparseGridGP,Integral2D>::Quadrature2D(const Integral2D &integral)
   : _integral(integral), _level(1), numGridPoints(0)
{
    _initSparseGrid();
}

template <typename Integral2D>
const typename Integral2D::T
Quadrature2D<SparseGridGP,Integral2D>::operator()(T ax, T bx, T ay, T by) const
{
    if ((ax == bx) || (ay == by)) {
        return 0.;
    }
    
    T result = 0.;
    for (int i=1; i<=_weights.numRows(); ++i) {
        T x1 = 0.5*( (bx-ax)*_knots(i,1) + (bx+ax));
        T x2 = 0.5*( (by-ay)*_knots(i,2) + (by+ay));
        result += _weights(i,1) * _integral.integrand(x1,x2);
    }
    result *= 0.25*(bx-ax)*(by-ay);
    return result;
}

template <typename Integral2D>
void
Quadrature2D<SparseGridGP,Integral2D>::setOrder(int order)
{
    _level = ceil(log2(order+1));
    _initSparseGrid();
}

template <typename Integral2D>
void
Quadrature2D<SparseGridGP,Integral2D>::setLevel(int level)
{
    _level = level;
    _initSparseGrid();
}

/*
 * The following routine contains code fragments from
 * file:     "sparse_grid_mixed_dataset.C"
 * routine:  "sparse_grid_mixed_dataset_handle"
 *
 * written by John Burkardt,
 * published under GNU LGPL and available online
 * location: "http://people.sc.fsu.edu/~jburkardt/cpp_src/cpp_src.html"
 */
template <typename Integral2D>
void
Quadrature2D<SparseGridGP,Integral2D>::_initSparseGrid()
{
    int dim_num = 2;
    int rule[2] = {3,3};     //3 = identifier for GP in libsparsegrid
    T alpha[2]  = {0.,0.}; //weight function = 1
    T beta[2]   = {0.,0.}; //weight function = 1
    int max_level = _level-1;
    T tol = 1e-16;
    int point_num;
    int point_total_num;
    int *sparse_index;
    int *sparse_order;
    double *sparse_point;
    int *sparse_unique_index;
    double *sparse_weight;

    point_total_num = webbur::sparse_grid_mixed_size_total ( dim_num, max_level, rule );

    point_num = webbur::sparse_grid_mixed_size ( dim_num, max_level, rule, alpha,
                                                 beta, tol );

    numGridPoints = point_num;

    sparse_unique_index = new int[point_total_num];

    webbur::sparse_grid_mixed_unique_index ( dim_num, max_level, rule, alpha, beta,
                                  tol, point_num, point_total_num, sparse_unique_index );

    sparse_order = new int[dim_num*point_num];
    sparse_index = new int[dim_num*point_num];

    webbur::sparse_grid_mixed_index ( dim_num, max_level, rule, point_num,
                    point_total_num, sparse_unique_index, sparse_order, sparse_index );

    //  Compute points and weights.
    sparse_point = new double [ dim_num * point_num ];

    webbur::sparse_grid_mixed_point ( dim_num, max_level, rule, alpha, beta,
                                  point_num, sparse_order, sparse_index, sparse_point );

    sparse_weight = new double[point_num];

    webbur::sparse_grid_mixed_weight ( dim_num, max_level, rule, alpha, beta,
                       point_num, point_total_num, sparse_unique_index, sparse_weight );

    _knots.engine().resize(point_num, dim_num);
    _weights.engine().resize(point_num, 1);

    for (int j = 0; j < point_num; j++ ) {
        for (int dim = 0; dim < dim_num; dim++ ) {
            _knots(j+1,dim+1) = sparse_point[dim+j*dim_num];
            _weights(j+1,1)   = sparse_weight[j];
        }
    }
    
    delete[] sparse_unique_index;
    delete[] sparse_order;
    delete[] sparse_index;
    delete[] sparse_point;
    delete[] sparse_weight;
}



template <typename Integral2D>
Quadrature2D<FullGridGL,Integral2D>::Quadrature2D(const Integral2D &integral)
   : _integral(integral), _order(1)
{
    _initFullGrid();
}

template <typename Integral2D>
const typename Integral2D::T
Quadrature2D<FullGridGL,Integral2D>::operator()(T ax, T bx, T ay, T by) const
{
    if ((ax == bx) || (ay == by))   return 0.;
    T result = 0.;
    for (int i=1; i<=_weights.numRows(); ++i) {
        T x1 = 0.5*( (bx-ax)*_knots(i,1) + (bx+ax));
        T x2 = 0.5*( (by-ay)*_knots(i,2) + (by+ay));
        result += _weights(i,1) * _integral.integrand(x1,x2);
    }
    result *= 0.25*(bx-ax)*(by-ay);
    return result;
}

template <typename Integral2D>
void
Quadrature2D<FullGridGL,Integral2D>::setOrder(int order)
{
    _order = order;
    _initFullGrid();
}

template <typename Integral2D>
void
Quadrature2D<FullGridGL,Integral2D>::_initFullGrid()
{
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots1d(_order,_order);
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights1d(_order,_order);
    T eps = Const<T>::EQUALITY_EPS;
    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=_order; ++k) {
       int     m = (k+1)/2;
       T xm = 0.5 * (x2+x1),
         xl = 0.5 * (x2-x1);
       for (int i=1; i<=m; ++i) {
           T z = cos(M_PI*(i-0.25)/(k+0.5)),
             z1, pp;
           do {
               T p1 = 1.0,
                 p2 = 2.0;
               for (int j=1; j<=k; ++j) {
                   T p3 = p2;
                     p2 = p1;
                     p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
               }
               pp = k * (z*p1-p2)/(z*z-1.0);
               z1 = z;
               z = z1-p1/pp;
           } while (fabs(z-z1) > eps);
           _knots1d(k,i)     = xm - xl*z;
           _knots1d(k,k+1-i) = xm + xl*z;
           _weights1d(k,i)     = 2.0*xl/((1.0-z*z)*pp*pp);
           _weights1d(k,k+1-i) = _weights1d(k,i);
        }
    }
    _knots.engine().resize(_order*_order, 2);
    _weights.engine().resize(_order*_order, 1);

    int count=1;
    for (int i = 1; i <= _order; ++i) {
        for (int j = 1; j <= _order; ++j) {
            _knots(count,1)   = _knots1d(_order,i);
            _knots(count,2)   = _knots1d(_order,j);
            _weights(count,1) = _weights1d(_order,i)*_weights1d(_order,j);
            ++count;
        }
    }
}

}    //namespace lawa

