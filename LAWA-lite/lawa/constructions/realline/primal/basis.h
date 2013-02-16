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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_H
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_H 1

#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/settings/enum.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Primal,R,CDF>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = R;
        static const Construction Cons = CDF;

        typedef BasisFunction<T,Primal,R,CDF> BasisFunctionType;
        typedef BSpline<T,Primal,R,CDF> BSplineType;
        typedef Wavelet<T,Primal,R,CDF> WaveletType;

        Basis(int _d, int _d_, int j=0);

        int 
        level() const;

        void
        setLevel(int j) const;

        const BasisFunctionType &
        generator(XType xtype) const;

        const int d, d_, j0;
        MRA<T,Primal,R,CDF> mra;
        MRA<T,Dual,R,CDF> mra_;
        Wavelet<T,Primal,R,CDF> psi;
        RefinementMatrix<T,R,CDF> M1;
        
    private:
        mutable int _j;
};

} // namespace lawa

#include <lawa/constructions/realline/primal/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_H

