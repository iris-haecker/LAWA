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

#ifndef LAWA_SETTINGS_ENUM_H
#define LAWA_SETTINGS_ENUM_H 1

namespace lawa {

enum FunctionSide {
    Primal,
    Dual,
    Orthogonal
};

enum DomainType {
    Periodic,
    R,
    Interval
};

enum Dimension {
    OneD = 1,
    TwoD = 2
};

enum XType {
    XBSpline = 0,
    XWavelet = 1
};

enum Construction {
    CDF,
    AnyInterval,
    DKU,
    Primbs,
    Dijkema,
    Multi
};

enum BoundaryCondition {
    NoBC = 0,
    DirichletBC = 1
};

enum QuadratureType { 
    Gauss, 
    Trapezoidal,
    SparseGridGP,    //Gauss-Patterson quadrature on a sparse grid
    FullGridGL        //Gauss-Legendre  quadrature on a full grid
};

enum SortingCriterion {
    AbsoluteValue,
    Lexicographical,
    Lexicographical2
    //Uniform
};

enum MethodType {
    Adaptive,
    SparseGrid,
    Uniform
};

} // namespace lawa

#endif // LAWA_SETTINGS_ENUM_H

