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

#ifndef ADAPTIVE_POSTPROCESSING_PLOTTING_H
#define ADAPTIVE_POSTPROCESSING_PLOTTING_H 1

#include <iostream>
#include <fstream>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template <typename T, typename Basis>
void
getSingularPoints(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
                  DenseVector<Array<T> > &sing_pts);

//Plot solution only on singular  points of the solution
template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
     const Preconditioner &P, T (*u)(T), const char* filename);

//Plot solution on a fixed grid
template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
     const Preconditioner &P, T (*u)(T), T (*du)(T), T a, T b, T h, const char* filename);

template <typename T, typename Basis2D, typename Preconditioner>
void
plot2D(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
       const Preconditioner &P, T (*u)(T,T), T a1, T b1, T a2, T b2, T h, const char* filename);
       
template <typename T, typename Basis2D, typename Preconditioner>
void
plot2D(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
       const Preconditioner &P, T (*u)(T,T), T (*dy_u)(T,T), T a1, T b1, T a2, T b2, 
       T h1, T h2, const char* filename);

template <typename T, DomainType Domain, Construction Cons>
void
plotCoeff(const Coefficients<AbsoluteValue,T,Index1D > &coeff, const Basis<T,Primal,Domain,Cons> &basis, const char* filename);

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotCoeff2D(const Coefficients<AbsoluteValue,T,Index> &coeff, const Basis_x &basis_x, const Basis_y &basis_y, const char* filename);

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotScatterCoeff2D(const Coefficients<AbsoluteValue,T,Index> &coeff, const Basis_x &basis_x,
                   const Basis_y &basis_y, const char* filename);

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotScatterCoeff2D(const Coefficients<Lexicographical,T,Index> &coeff, const Basis_x &basis_x,
                   const Basis_y &basis_y, const char* filename);


}  // namespace lawa

#include <lawa/methods/adaptive/postprocessing/plotting.tcc>

#endif // ADAPTIVE_POSTPROCESSING_PLOTTING_H

