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

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV_TENSOR_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV_TENSOR_H 1

#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/algorithms/apply1d.h>

namespace lawa {

template <typename T, typename Basis, typename APPLY_TENSOR, typename RHS>
class GHS_ADWAV_TENSOR {

        typedef typename IndexSet<Index2D>::iterator                             set_it;
        typedef typename IndexSet<Index2D>::const_iterator                       const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
        typedef typename Coefficients<AbsoluteValue,T,Index2D >::const_iterator  const_coeff_abs_it;
        typedef typename Coefficients<Lexicographical,T,Index2D>::value_type     val_type;

        typedef typename APPLY_TENSOR::MAType MA;


    public:

        GHS_ADWAV_TENSOR(const Basis &_basis, APPLY_TENSOR &Apply, RHS &_F);

        Coefficients<Lexicographical,T,Index2D>
        SOLVE(T nuM1, T _eps, int NumOfIterations=100, T H1norm=T(0));

        std::vector<Coefficients<Lexicographical,T,Index2D> > solutions;
        std::vector<T>                                        residuals;
        std::vector<T>                                        times;
        std::vector<int>                                      linsolve_iterations;


    private:
        const Basis         &basis;
        APPLY_TENSOR        &Apply;
        MA                  &A;
        RHS                 &F;
        T                   cA, CA, kappa;
        T                   alpha, omega, gamma, theta;
        T                   eps;

        IndexSet<Index2D>
        GROW(const Coefficients<Lexicographical,T,Index2D>      &w,
             T                                                  nu_bar,
             T                                                  &nu);

        Coefficients<Lexicographical,T,Index2D>
        GALSOLVE(const IndexSet<Index2D>                        &Lambda,
                 const Coefficients<Lexicographical,T,Index2D>  &g,
                 const Coefficients<Lexicographical,T,Index2D>  &w,
                 T                                              delta,
                 T                                              tol);
};

}    //namespace lawa

#include <lawa/methods/adaptive/solvers/ghsadwav_tensor.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV_TENSOR_H

