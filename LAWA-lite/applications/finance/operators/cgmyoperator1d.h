/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Sebastian Kestler, Kristina Steih, Mario Rometsch, Alexander Stippler.

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

#ifndef APPLICATIONS_FINANCE_OPERATORS_CGMYOPERATOR1D_H
#define APPLICATIONS_FINANCE_OPERATORS_CFGMYOPERATOR1D_H 1

#include <applications/finance/kernels/cgmykernel.h>
#include <applications/finance/operators/financeoperator.h>
#include <applications/finance/processes/processes.h>

namespace lawa {

template <typename T, typename Basis>
struct FinanceOperator1D<T, CGMY, Basis>
{
    typedef typename Basis::BSplineType PrimalSpline;
    typedef typename Basis::WaveletType PrimalWavelet;

    FinanceOperator1D(const Parameters<T,CGMY> &_params, const Basis& _basis, T _R1=0., T _R2=1.);

    T
    operator()(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const;

    T
    operator()(const Index1D &row_index, const Index1D &col_index) const;

    const Parameters<T,CGMY>   &params;
    const Basis               &basis;
    Kernel<T,CGMY>            kernel;
    PrimalSpline phi, d_phi, delta_phi;
    PrimalWavelet psi, d_psi, delta_psi;
    T R1, R2;
    T OneDivSqrtR2pR1, OneDivR2pR1, R1DivR1pR2;

    Integral<Gauss, Basis, Basis> integral;

    mutable std::map<T,T> values_tailintegral;

};

}   // namespace lawa

#include <applications/finance/operators/cgmyoperator1d.tcc>

#endif  // APPLICATIONS_FINANCE_OPERATORS_FINANCEOPERATOR1D_H
