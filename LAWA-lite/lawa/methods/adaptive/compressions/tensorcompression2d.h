/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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


#ifndef  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_TENSORCOMPRESSION2D_H
#define  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_TENSORCOMPRESSION2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/aux/timer.h>

namespace lawa {

template <typename T, typename Compression_x, typename Compression_y>
struct TensorCompression2D
{
    Compression_x &c_x;
    Compression_y &c_y;

    IndexSet<Index1D> LambdaRow_x;
    IndexSet<Index1D> LambdaRow_y;
    IndexSet<Index2D> LambdaRow;

    TensorCompression2D(Compression_x &_c_x, Compression_y &_c_y);

    void
    setParameters(const IndexSet<Index2D> &_LambdaRow);

    IndexSet<Index2D>
    SparsityPattern(const Index2D &lambda_col, const IndexSet<Index2D> &_Lambda);
};

} // namespace lawa

#include <lawa/methods/adaptive/compressions/tensorcompression2d.tcc>

#endif //  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_TENSORCOMPRESSION2D_H

