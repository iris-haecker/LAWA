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

#ifndef LAWA_MATH_POW2_H
#define LAWA_MATH_POW2_H 1

#include <lawa/settings/typetraits.h>

namespace lawa {

template <typename T>
    typename RestrictTo<flens::IsSame<T,int>::value
                     || flens::IsSame<T,long>::value, T>::Type
    pow2i(int expo);

template <typename T>
    typename RestrictTo<flens::IsSame<T,double>::value, T>::Type
    pow2i(int expo);

template <typename T>
T
pow2ih(int expo);

} // namespace lawa

#include <lawa/math/pow2.tcc>

#endif //  LAWA_MATH_POW2_H

