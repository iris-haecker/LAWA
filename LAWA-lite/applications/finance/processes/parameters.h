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


#ifndef APPLICATIONS_FINANCE_PROCESSES_PARAMETERS_H
#define APPLICATIONS_FINANCE_PROCESSES_PARAMETERS_H 1


#include <applications/finance/processes/parameters.h>

namespace lawa {

template < typename T, ProcessType Type>
struct Parameters
{

};

template <typename T>
struct Parameters<T,CGMY>
{
  Parameters(T _r, T _k_C, T _k_G, T _k_M, T _k_Y);

  T r;
  T k_C;
  T k_G;
  T k_M;
  T k_Y;
};

}   // namespace lawa

#include <applications/finance/processes/parameters.tcc>

#endif  // APPLICATIONS_FINANCE_PROCESSES_PARAMETERS_H
