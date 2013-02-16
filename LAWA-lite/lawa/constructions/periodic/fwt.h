/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_FWT_H
#define LAWA_CONSTRUCTIONS_PERIODIC_FWT_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/periodic/primal/basis.h>
#include <lawa/constructions/periodic/dual/basis.h>

namespace lawa {

template <typename X, typename Y>
    void
    decompose(const DenseVector<X> &x, 
              const Basis<typename X::ElementType,Dual,Periodic,CDF> &basis_, int j,
              DenseVector<Y> &y);

template <typename X, typename Y>
    void
    reconstruct(const DenseVector<X> &x, 
                const Basis<typename X::ElementType,Primal,Periodic,CDF> &basis, int j,
                DenseVector<Y> &y);

template <typename X, typename Y>
    void
    fwt(const DenseVector<X> &x, 
        const Basis<typename X::ElementType,Dual,Periodic,CDF> &basis_, int j,
        DenseVector<Y> &y);

template <typename X, typename Y>
    void
    ifwt(const DenseVector<X> &x, 
         const Basis<typename X::ElementType,Primal,Periodic,CDF> &basis, int j,
         DenseVector<Y> &y);

} // namespace lawa

#include <lawa/constructions/periodic/fwt.tcc>

#endif // LAWA_CONSTRUCTIONS_PERIODIC_FWT_H

