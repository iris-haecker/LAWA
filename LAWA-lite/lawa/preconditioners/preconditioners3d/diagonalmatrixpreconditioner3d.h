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


#ifndef LAWA_PRECONDITIONERS_PRECONDITIONERS3D_DIAGONALMATRIXPRECONDITIONER3D_H
#define LAWA_PRECONDITIONERS_PRECONDITIONERS3D_DIAGONALMATRIXPRECONDITIONER3D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integrals.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename Basis3D, typename BilinearForm>
class DiagonalMatrixPreconditioner3D
{

    public:
        DiagonalMatrixPreconditioner3D(const BilinearForm &a);

        T
        operator()(XType xtype1, int j1, int k1, 
                   XType xtype2, int j2, int k2,
                   XType xtype3, int j3, int k3) const;

        T
        operator()(const Index3D &index) const;

    private:
        const BilinearForm &_a;
};

}   // namespace lawa

#include <lawa/preconditioners/preconditioners3d/diagonalmatrixpreconditioner3d.tcc>

#endif // LAWA_PRECONDITIONERS_PRECONDITIONERS3D_DIAGONALMATRIXPRECONDITIONER3D_H

