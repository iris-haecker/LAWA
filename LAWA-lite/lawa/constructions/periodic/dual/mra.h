/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_DUAL_MRA_H
#define LAWA_CONSTRUCTIONS_PERIODIC_DUAL_MRA_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/periodic/dual/bspline.h>
#include <lawa/constructions/periodic/refinementmatrix.h>

namespace lawa {

template <typename _T>
class MRA<_T,Dual,Periodic,CDF>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Dual;
        static const DomainType Domain = Periodic;
        static const Construction Cons = CDF;

        typedef BasisFunction<T,Dual,Periodic,CDF> BasisFunctionType;
        typedef BSpline<T,Dual,Periodic,CDF> BSplineType;

        MRA(int _d, int _d_, int j=0);

        int
        level() const;

        void
        setLevel(int j) const;

        int
        cardI_(int j) const;

        Range<int>
        rangeI_(int j) const;

        const int d, d_, j0;
        BSpline<T,Dual,Periodic,CDF> phi_;
        RefinementMatrix<T,Periodic,CDF> M0_;
        
    protected:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/periodic/dual/mra.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_DUAL_MRA_H

