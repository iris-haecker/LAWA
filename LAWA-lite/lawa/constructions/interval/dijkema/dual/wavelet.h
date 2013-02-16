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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_WAVELET_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_WAVELET_H 1

#include <lawa/constructions/support.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
struct Wavelet<_T,Dual,Interval,Dijkema>
    : public BasisFunction<_T,Dual,Interval,Dijkema>
{
    typedef _T T;
    static const FunctionSide Side = Dual;
    static const DomainType Domain = Interval;
    static const Construction Cons = Dijkema;

    Wavelet(const Basis<T,Dual,Interval,Dijkema> &_basis_);

    T
    operator()(T x, int j, long k, unsigned short deriv=0) const;

    Support<T>
    support(int j, long k) const;

    const Basis<T,Dual,Interval,Dijkema> &basis_;
};

} // namespace lawa

#include <lawa/constructions/interval/dijkema/dual/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_WAVELET_H

