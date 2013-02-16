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

#ifndef LAWA_AUX_ARROW_H
#define LAWA_AUX_ARROW_H 1

#include <lawa/flensforlawa.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <cxxblas::StorageOrder Order, typename T>
    void
    arrow(const flens::GeMatrix<flens::FullStorage<T, Order> > &In,
                flens::GeMatrix<flens::FullStorage<T, Order> > &Out);

} // namespace lawa

#include <lawa/aux/arrow.tcc>

#endif // LAWA_AUX_ARROW_H

