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

#include <cassert>
#include <cmath>
#include <list>

#include <lawa/aux/arrow.h>
#include <lawa/math/math.h>
#include <lawa/constructions/realline/primal/bspline.h>
#include <extensions/flens/lapack_flens.h>
#include <lawa/constructions/interval/primbs/primal/spline_helper.h>

namespace lawa {

template <typename T>
MRA<T,Primal,Interval,Primbs>::MRA(int _d, int j)
    : d(_d), mu(d&1),
      min_j0(iceil(log2(2*d))),
      j0((j==-1) ? min_j0 : j), phiR(d),
      l1((mu-d)/2), l2((mu+d)/2),
      _bc(2,0), _j(j0), phi(*this)
{
    assert(d>1);
    assert(_j>=min_j0);
    
    _calcM0();
    setLevel(_j);
}

template <typename T>
MRA<T,Primal,Interval,Primbs>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
int
MRA<T,Primal,Interval,Primbs>::cardI(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) + d - 1 - (_bc(0)+_bc(1));
}

template <typename T>
int
MRA<T,Primal,Interval,Primbs>::cardIL(int /*j*/) const
{
    return d - 1 -_bc(0);
}

template <typename T>
int
MRA<T,Primal,Interval,Primbs>::cardII(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) - d - 1;
}

template <typename T>
int
MRA<T,Primal,Interval,Primbs>::cardIR(int /*j*/) const
{
    return d - 1 - _bc(1);
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
Range<int>
MRA<T,Primal,Interval,Primbs>::rangeI(int j) const
{
    assert(j>=min_j0);
    return Range<int>(1 + _bc(0), pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
Range<int>
MRA<T,Primal,Interval,Primbs>::rangeIL(int /*j*/) const
{
    return Range<int>(1 + _bc(0), d - 1);
}

template <typename T>
Range<int>
MRA<T,Primal,Interval,Primbs>::rangeII(int j) const
{
    assert(j>=min_j0);
    return Range<int>(d, pow2i<T>(j));
}

template <typename T>
Range<int>
MRA<T,Primal,Interval,Primbs>::rangeIR(int j) const
{
    assert(j>=min_j0);
    return Range<int>(pow2i<T>(j) + 1, pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
int
MRA<T,Primal,Interval,Primbs>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,Interval,Primbs>::setLevel(int j) const
{
    assert(j>=min_j0);
    _j = j;
    M0.setLevel(_j);
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Primal,Interval,Primbs>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;
    _calcM0();
}

template <typename T>
void
MRA<T,Primal,Interval,Primbs>::_calcM0()
{
    typedef GeMatrix<FullStorage<T,cxxblas::ColMajor> > FullColMatrix;
    FullColMatrix R(_(1-l2, pow2i<T>(min_j0)-l1-1),
                    _(1-l2, pow2i<T>(min_j0)-l1-1));
    R.diag(0) = 1.;

    // replace this part as soon as well understood ;-)
    DenseVector<Array<T> > knots = linspace(-d+1., d-1., 2*d-1);
    knots.engine().changeIndexBase(1);
    FullColMatrix Transformation(knots.length()-d, knots.length()-d);
    Transformation.diag(0) = 1.;
    for (int i=1; i<d; ++i) {
        FullColMatrix Tmp = insertKnot(d-1,knots,0.), Tmp2;
        blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,
                 1.,Tmp,Transformation,0.,Tmp2);
        Transformation = Tmp2;
    }

    FullColMatrix InvTrans = Transformation(_(Transformation.lastRow()-(d-1)+1,
                                              Transformation.lastRow()), _ );
    //--- inverse(InvTrans)
    FullColMatrix TransTmp = InvTrans, Trans, TransTmp2;
    Trans = inv(InvTrans);

    Trans.engine().changeIndexBase(R.firstRow(),R.firstCol());
    R(_(Trans.firstRow(),Trans.lastRow()),
      _(Trans.firstCol(),Trans.lastCol())) = Trans;

    arrow(Trans,TransTmp);
    TransTmp.engine().changeIndexBase(R.lastRow()-TransTmp.numRows()+1,
                                      R.lastCol()-TransTmp.numCols()+1);
    R(_(TransTmp.firstRow(),TransTmp.lastRow()),
      _(TransTmp.firstCol(),TransTmp.lastCol())) = TransTmp;

    FullColMatrix ExtM0(_(-l2+1,pow2i<T>(min_j0+1)-l1-1),
                        _(-l2+1,pow2i<T>(min_j0)-l1-1));
    for (int q=ExtM0.firstCol(); q<=ExtM0.lastCol(); ++q) {
        for (int p=std::max(l1+2*q,ExtM0.firstRow());
             p<=std::min(l2+2*q,ExtM0.lastRow()); ++p) {
            ExtM0(p,q) = phiR.a(p-2*q);
        }
    }

    FullColMatrix RjPlus1(_(-l2+1,pow2i<T>(min_j0+1)-l1-1),
                          _(-l2+1,pow2i<T>(min_j0+1)-l1-1));
    RjPlus1.diag(0) = 1.;
    Trans.engine().changeIndexBase(RjPlus1.firstRow(),RjPlus1.firstCol());
    RjPlus1(Trans) = Trans;
    TransTmp.engine().changeIndexBase(RjPlus1.lastRow()-TransTmp.numRows()+1,
                                      RjPlus1.lastCol()-TransTmp.numCols()+1);
    RjPlus1(TransTmp) = TransTmp;

    FullColMatrix Mj0;
    FullColMatrix M0Tmp, Tmp = RjPlus1;
    FullColMatrix TmpTmp = Tmp;
    
    Tmp = inv(TmpTmp);

    blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,Tmp,ExtM0,0.,M0Tmp);
    blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,1.,M0Tmp,R,0.,Mj0);
    blas::scal(Const<T>::R_SQRT2, Mj0);

    // TODO: integrate boundary conditions in calculation process.
    if ((_bc(0)==DirichletBC) && (_bc(1)==DirichletBC)) {
        FullColMatrix Mj0withBC = Mj0( _(Mj0.firstRow()+1,Mj0.lastRow()-1) , 
                                       _(Mj0.firstCol()+1,Mj0.lastCol()-1));
        Mj0withBC.engine().changeIndexBase(Mj0.firstRow()+1, Mj0.firstCol()+1);
        M0 = RefinementMatrix<T,Interval,Primbs>(d-1-_bc(0), d-1-_bc(1), 
                                                 Mj0withBC, min_j0, min_j0);
    } else {
        M0 = RefinementMatrix<T,Interval,Primbs>(d-1, d-1, Mj0, min_j0, min_j0);
    }
    M0.setLevel(_j);
}

} // namespace lawa

