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

#ifndef LAWA_INTEGRALS_QUADRATURE2D_H
#define LAWA_INTEGRALS_QUADRATURE2D_H 1

#include <extensions/sparsegrid/sparse_grid_mixed.h>
#include <lawa/math/math.h>

namespace lawa {

template <QuadratureType Quad, typename Integral2D>
class Quadrature2D
{
};

template <typename Integral2D>
class Quadrature2D<SparseGridGP,Integral2D>
{
    public:
        typedef typename Integral2D::T T;

        Quadrature2D(const Integral2D &integral);

        const T
        operator()(T ax, T bx, T ay, T by) const;

        void setOrder(int order);
        void setLevel(int level);

        int numGridPoints;

    private:
        const Integral2D &_integral;
        int _level;

        void
        _initSparseGrid();

        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;
};

template <typename Integral2D>
class Quadrature2D<FullGridGL,Integral2D>
{
    public:
        typedef typename Integral2D::T T;
        
        Quadrature2D(const Integral2D &integral);

        const T
        operator()(T ax, T bx, T ay, T by) const;

        void setOrder(int order);

    private:
        const Integral2D &_integral;
        int _order;

        void
        _initFullGrid();

        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;
};

}   //namespace lawa

#include <lawa/integrals/quadrature2d.tcc>

#endif   //LAWA_RIGHTHANDSIDES_QUADRATURE2D_H

