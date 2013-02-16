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
#include <lawa/constructions/multi/interval/_wavelet_evaluator.h>

namespace lawa {

template <typename T>
Basis<T,Orthogonal,Interval,Multi>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _bc(2,0), _j(j0), psi(*this)
{
    assert(d>=2);
    
    switch (d) {
        case 2:
            // left wavelets
            _numLeftParts = 3;
            _leftEvaluator = new Evaluator[3];
            _leftEvaluator[0] = _linear_wavelet_left_evaluator0;
            _leftEvaluator[1] = _linear_wavelet_left_evaluator1;
            _leftEvaluator[2] = _linear_wavelet_inner_evaluator2;
            
            _leftSupport = new Support<T>[3];
            _leftSupport[0] = Support<T>(0,1);
            _leftSupport[1] = Support<T>(0,0.5);
            _leftSupport[2] = Support<T>(0,1);
            
            _leftSingularSupport = new DenseVector<Array<T> >[3];
            _leftSingularSupport[0].engine().resize(7,0);
            _leftSingularSupport[0] = 0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _leftSingularSupport[1] = linspace(0.0,1.0,5);
            _leftSingularSupport[2] = linspace(0.0,1.0,9);
            
            _leftScalingFactors.engine().resize(3,0);
            _leftScalingFactors = 1.0,1.0,1.0;
            
            // inner wavelets
            _numInnerParts = 3; 
            _innerEvaluator = new Evaluator[3];
            _innerEvaluator[0] = _linear_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _linear_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _linear_wavelet_inner_evaluator2;
            
            _innerSupport = new Support<T>[3];
            _innerSupport[0] = Support<T>(-1,1);
            _innerSupport[1] = Support<T>(-1,1);
            _innerSupport[2] = Support<T>( 0,1);
            
            _innerSingularSupport = new DenseVector<Array<T> >[3];
            _innerSingularSupport[0].engine().resize(13,0);
            _innerSingularSupport[0] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _innerSingularSupport[1].engine().resize(13,0);
            _innerSingularSupport[1] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _innerSingularSupport[2] = linspace(0.0,1.0,9);
            
            // right wavelets
            _numRightParts = 2;
            _rightEvaluator = new Evaluator[2];
            _rightEvaluator[0] = _linear_wavelet_right_evaluator0;
            _rightEvaluator[1] = _linear_wavelet_right_evaluator1;
            
            _rightSupport = new Support<T>[2];
            _rightSupport[0] = Support<T>(-1.0,0.0);
            _rightSupport[1] = Support<T>(-0.5,0.0);
            
            _rightSingularSupport = new DenseVector<Array<T> >[2];
            _rightSingularSupport[0].engine().resize(7,0);
            _rightSingularSupport[0] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0;
            _rightSingularSupport[1] = linspace(-1.0,0.0,5);
            
            _rightScalingFactors.engine().resize(2,0);
            _rightScalingFactors = 1.0,1.0;
            
        case 4:
            break;
            
        default: std::cerr << "Wavelet<T,Orthogonal,Interval,Multi> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }

    setLevel(_j);
}
    
template <typename T>
Basis<T,Orthogonal,Interval,Multi>::~Basis()
{
    delete[] _leftEvaluator;
    delete[] _innerEvaluator;
    delete[] _rightEvaluator;
    delete[] _leftSupport;
    delete[] _innerSupport;
    delete[] _rightSupport;
    delete[] _leftSingularSupport;
    delete[] _innerSingularSupport;
    delete[] _rightSingularSupport;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Orthogonal,Interval,Multi>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Orthogonal,Interval,Multi>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);
    _bc = 1,1;
    
    switch (d) {
        case 2:
            // left wavelets
            _numLeftParts = 2;
            delete[] _leftEvaluator;
            _leftEvaluator = new Evaluator[2];
            _leftEvaluator[0] = _linear_wavelet_inner_evaluator1;
            _leftEvaluator[1] = _linear_wavelet_inner_evaluator2;
            
            delete[] _leftSupport;
            _leftSupport = new Support<T>[2];
            _leftSupport[0] = Support<T>(0,1);
            _leftSupport[1] = Support<T>(0,1);
            
            delete[] _leftSingularSupport;
            _leftSingularSupport = new DenseVector<Array<T> >[2];
            _leftSingularSupport[0].engine().resize(7,0);
            _leftSingularSupport[0] = 0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _leftSingularSupport[1] = linspace(0.0,1.0,9);
            
            _leftScalingFactors.engine().resize(2,0);
            _leftScalingFactors = 1.231176897368409556535,1.0;
                      
            // right wavelets
            _numRightParts = 1;
            delete[] _rightEvaluator;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _linear_wavelet_inner_evaluator1;
            
            delete[] _rightSupport;
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(-1.0,0.0);

            delete[] _rightSingularSupport;
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(7,0);
            _rightSingularSupport[0] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0;
            
            _rightScalingFactors.engine().resize(1,0);
            _rightScalingFactors = 1.714277953522051005391;
            
            break;
        
        case 4:            
            // left wavelets
            _numLeftParts = 4;
            _leftEvaluator = new Evaluator[4];
            _leftEvaluator[0] = _cubic_wavelet_left_evaluator1;
            _leftEvaluator[1] = _cubic_wavelet_left_evaluator2;
            _leftEvaluator[2] = _cubic_wavelet_inner_evaluator4;
            _leftEvaluator[3] = _cubic_wavelet_inner_evaluator5;

            _leftSupport = new Support<T>[4];
            _leftSupport[0] = Support<T>(0,1);
            _leftSupport[1] = Support<T>(0,1);
            _leftSupport[2] = Support<T>(0,1);
            _leftSupport[3] = Support<T>(0,1);

            _leftSingularSupport = new DenseVector<Array<T> >[4];
            _leftSingularSupport[0].engine().resize(9,0);
            _leftSingularSupport[0] = 0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0;
            _leftSingularSupport[1].engine().resize(9,0);
            _leftSingularSupport[1] = 0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0;
            _leftSingularSupport[2].engine().resize(9,0);
            _leftSingularSupport[2] = 0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0;
            _leftSingularSupport[3] = linspace(0.0,1.0,9);

            _leftScalingFactors.engine().resize(4,0);
            _leftScalingFactors = 1.0,1.0,1.0,1.0;

            // inner wavelets
            _numInnerParts = 6; 
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _cubic_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _cubic_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _cubic_wavelet_inner_evaluator2;
            _innerEvaluator[3] = _cubic_wavelet_inner_evaluator3;
            _innerEvaluator[4] = _cubic_wavelet_inner_evaluator4;
            _innerEvaluator[5] = _cubic_wavelet_inner_evaluator5;

            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(-1,1);
            _innerSupport[1] = Support<T>(-1,1);
            _innerSupport[2] = Support<T>(-1,1);
            _innerSupport[3] = Support<T>(-1,1);
            _innerSupport[4] = Support<T>(-1,1);
            _innerSupport[5] = Support<T>( 0,1);

            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(13,0);
            _innerSingularSupport[0] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _innerSingularSupport[1].engine().resize(13,0);
            _innerSingularSupport[1] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _innerSingularSupport[2].engine().resize(13,0);
            _innerSingularSupport[2] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _innerSingularSupport[3].engine().resize(13,0);
            _innerSingularSupport[3] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,0.125,0.25,0.375,0.5,0.75,1.0;
            _innerSingularSupport[4] = linspace(0.0,1.0,9);
            _innerSingularSupport[5] = linspace(0.0,1.0,9);

            // right wavelets
            _numRightParts = 2;
            _rightEvaluator = new Evaluator[2];
            _rightEvaluator[0] = _cubic_wavelet_right_evaluator1;
            _rightEvaluator[1] = _cubic_wavelet_right_evaluator2;

            _rightSupport = new Support<T>[2];
            _rightSupport[0] = Support<T>(-1.0,0.0);
            _rightSupport[1] = Support<T>(-1.0,0.0);

            _rightSingularSupport = new DenseVector<Array<T> >[2];
            _rightSingularSupport[0] = linspace(-1.0,0.0,9);
            _rightSingularSupport[1] = linspace(-1.0,0.0,9);

            _rightScalingFactors.engine().resize(2,0);
            _rightScalingFactors = 1.0,1.0;
            
            break;
            
        default: std::cerr << "Boundary conditions not realized yet "
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
    
    mra.enforceBoundaryCondition<BC>();
}

template <typename T>
const BasisFunction<T,Orthogonal,Interval,Multi> &
Basis<T,Orthogonal,Interval,Multi>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi; 
    } else {
        return psi;
    }
}

// cardinalities of whole, left, inner, right index sets (primal).
template <typename T>
long
Basis<T,Orthogonal,Interval,Multi>::cardJ(int j) const
{
    assert(j>=j0);
    return _numLeftParts + (pow2i<long>(j)-1)*_numInnerParts + _numRightParts;
}

template <typename T>
long
Basis<T,Orthogonal,Interval,Multi>::cardJL(int j) const
{
    assert(j>=j0 or j==-1);
    return _numLeftParts;
}

template <typename T>
long
Basis<T,Orthogonal,Interval,Multi>::cardJI(int j) const
{
    assert(j>=j0);
    return (pow2i<long>(j)-1)*_numInnerParts;
}

template <typename T>
long
Basis<T,Orthogonal,Interval,Multi>::cardJR(int j) const
{
    assert(j>=j0);
    return _numRightParts;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJ(int j) const
{
    assert(j>=j0);
    return Range<int>(1,cardJ(j));
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJL(int j) const
{
    assert(j>=j0 or j==-1);
    return Range<int>(1,cardJL());
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJI(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL()+1,cardJL()+cardJI(j));
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJR(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL()+cardJI(j)+1,cardJ(j));
}

} // namespace lawa
