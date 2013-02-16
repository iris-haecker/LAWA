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
#include <lawa/constructions/multi/interval/_bspline_evaluator.h>

namespace lawa {

template <typename T>
MRA<T,Orthogonal,Interval,Multi>::MRA(int _d, int j)
    : d(_d), j0((j==-1) ? 0 : j), _bc(2,0), _j(j0), phi(*this)
{
    assert(d>=2);
    assert(j0>=0);
    
    switch (d) {
        case 2: 
            // left B-splines 
            _numLeftParts = 3;
            _leftEvaluator = new Evaluator[3];
            _leftEvaluator[0] = _linear_bspline_left_evaluator0;
            _leftEvaluator[1] = _linear_bspline_inner_evaluator1;
            _leftEvaluator[2] = _linear_bspline_inner_evaluator2;
            
            _leftSupport = new Support<T>[3];
            _leftSupport[0] = Support<T>(0,1);
            _leftSupport[1] = Support<T>(0,1);
            _leftSupport[2] = Support<T>(0,1);
            
            _leftSingularSupport = new DenseVector<Array<T> >[3];
            _leftSingularSupport[0] = linspace(0.0,1.0,5);
            _leftSingularSupport[1] = linspace(0.0,1.0,3);
            _leftSingularSupport[2] = linspace(0.0,1.0,5);
            
            _leftScalingFactors.engine().resize(3,0);
            _leftScalingFactors = 1.714277953522051005391,1.0,1.0;
            
            // inner B-splines 
            _numInnerParts = 3;
            _innerEvaluator = new Evaluator[3];
            _innerEvaluator[0] = _linear_bspline_inner_evaluator0;
            _innerEvaluator[1] = _linear_bspline_inner_evaluator1;
            _innerEvaluator[2] = _linear_bspline_inner_evaluator2;            
            
            _innerSupport = new Support<T>[3];
            _innerSupport[0] = Support<T>(-1,1);
            _innerSupport[1] = Support<T>( 0,1);
            _innerSupport[2] = Support<T>( 0,1);
            
            _innerSingularSupport = new DenseVector<Array<T> >[3];
            _innerSingularSupport[0] = linspace(-1.0,1.0,9);
            _innerSingularSupport[1] = linspace(0.0,1.0,3);
            _innerSingularSupport[2] = linspace(0.0,1.0,5);
            
            // right B-splines
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _linear_bspline_right_evaluator0;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(-1,0);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0] = linspace(-1.0,0.0,5);
            
            _rightScalingFactors.engine().resize(1,0);
            _rightScalingFactors(0) = 1.231176897368409556535;
            
            break;
            
        case 3:/* // Works only for Dirichlet BC so far !!!
                // just inner parts set here!!!          
            // inner B-splines 
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _quadratic_bspline_inner_evaluator4;
            _innerEvaluator[1] = _quadratic_bspline_inner_evaluator5;
            _innerEvaluator[2] = _quadratic_bspline_inner_evaluator0;            
            _innerEvaluator[3] = _quadratic_bspline_inner_evaluator1;
            _innerEvaluator[4] = _quadratic_bspline_inner_evaluator2;
            _innerEvaluator[5] = _quadratic_bspline_inner_evaluator3;            
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(-1,1);
            _innerSupport[1] = Support<T>(-1,1);
            _innerSupport[2] = Support<T>( 0,1);
            _innerSupport[3] = Support<T>( 0,1);
            _innerSupport[4] = Support<T>( 0,1);
            _innerSupport[5] = Support<T>( 0,1);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0] = linspace(-1.0,1.0,9);
            _innerSingularSupport[1] = linspace(0.0,1.0,3);
            _innerSingularSupport[2] = linspace(0.0,1.0,5);
            _innerSingularSupport[3] = linspace(0.0,1.0,5);
            _innerSingularSupport[4] = linspace(0.0,1.0,5);
            _innerSingularSupport[5] = linspace(0.0,1.0,5);
            */
        case 4:
            break;
            
        default: std::cerr << "BSpline<T,Orthogonal,Interval,Multi> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
    
    setLevel(_j);
}

template <typename T>
MRA<T,Orthogonal,Interval,Multi>::~MRA()
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

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
long
MRA<T,Orthogonal,Interval,Multi>::cardI(int j) const
{
    assert(j>=j0);
    return _numLeftParts + (pow2i<long>(j)-1)*_numInnerParts + _numRightParts;
}

template <typename T>
long
MRA<T,Orthogonal,Interval,Multi>::cardIL(int /*j*/) const
{
    return _numLeftParts;
}

template <typename T>
long
MRA<T,Orthogonal,Interval,Multi>::cardII(int j) const
{
    assert(j>=j0);
    return (pow2i<long>(j)-1)*_numInnerParts;
}

template <typename T>
long
MRA<T,Orthogonal,Interval,Multi>::cardIR(int /*j*/) const
{
    return _numRightParts;
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeI(int j) const
{
    assert(j>=j0);
    return Range<int>(0,cardI(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeIL(int /*j*/) const
{
    return Range<int>(0,cardIL() - 1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeII(int j) const
{
    assert(j>=j0);
    return Range<int>(cardIL(), cardIL()+cardII(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeIR(int j) const
{
    assert(j>=j0);
    return Range<int>(cardIL()+cardII(j),cardI(j)-1);
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Orthogonal,Interval,Multi>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Orthogonal,Interval,Multi>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;
    
    switch (d) {
        case 2: 
            // left B-splines 
            _numLeftParts = 2;
            delete[] _leftEvaluator;
            _leftEvaluator = new Evaluator[2];
            _leftEvaluator[0] = _linear_bspline_inner_evaluator1;
            _leftEvaluator[1] = _linear_bspline_inner_evaluator2;
            
            delete[] _leftSupport;
            _leftSupport = new Support<T>[2];
            _leftSupport[0] = Support<T>(0,1);
            _leftSupport[1] = Support<T>(0,1);

            delete[] _leftSingularSupport;
            _leftSingularSupport = new DenseVector<Array<T> >[2];
            _leftSingularSupport[0] = linspace(0.0,1.0,3);
            _leftSingularSupport[1] = linspace(0.0,1.0,5);
            
            _leftScalingFactors.engine().resize(2,0);
            _leftScalingFactors = 1.0,1.0;
            
            // right B-splines
            _numRightParts = 0;
            
            break;
            
        case 3:
            // left B-splines 
            _numLeftParts = 2;
            _leftEvaluator = new Evaluator[2];
            _leftEvaluator[0] = _linear_bspline_inner_evaluator1;
            _leftEvaluator[1] = _linear_bspline_inner_evaluator2;
            
            _leftSupport = new Support<T>[2];
            _leftSupport[0] = Support<T>(0,1);
            _leftSupport[1] = Support<T>(0,1);
            
            _leftSingularSupport = new DenseVector<Array<T> >[2];
            _leftSingularSupport[0] = linspace( 0.0,1.0,3);
            _leftSingularSupport[1] = linspace( 0.0,1.0,5);
            
            _leftScalingFactors.engine().resize(2,0);
            _leftScalingFactors = 1.0,1.0;
            
            // right B-splines
            _numRightParts = 0;
            
            break;
            
        case 4:
            // left B-splines 
            _numLeftParts = 5;
            _leftEvaluator = new Evaluator[5];
            _leftEvaluator[0] = _cubic_bspline_left_evaluator1;
            _leftEvaluator[1] = _cubic_bspline_inner_evaluator2;
            _leftEvaluator[2] = _cubic_bspline_inner_evaluator3;
            _leftEvaluator[3] = _cubic_bspline_inner_evaluator4;
            _leftEvaluator[4] = _cubic_bspline_inner_evaluator5;
            
            _leftSupport = new Support<T>[5];
            _leftSupport[0] = Support<T>(0,1);
            _leftSupport[1] = Support<T>(0,1);
            _leftSupport[2] = Support<T>(0,1);
            _leftSupport[3] = Support<T>(0,1);
            _leftSupport[4] = Support<T>(0,1);
            
            _leftSingularSupport = new DenseVector<Array<T> >[5];
            _leftSingularSupport[0] = linspace( 0.0,1.0,5);
            _leftSingularSupport[1] = linspace( 0.0,1.0,3);
            _leftSingularSupport[2] = linspace( 0.0,1.0,3);
            _leftSingularSupport[3] = linspace( 0.0,1.0,5);
            _leftSingularSupport[4] = linspace( 0.0,1.0,5);
            
            _leftScalingFactors.engine().resize(5,0);
            _leftScalingFactors = 1.0,1.0,1.0,1.0,1.0;
            
            // inner B-splines 
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _cubic_bspline_inner_evaluator0;
            _innerEvaluator[1] = _cubic_bspline_inner_evaluator1;
            _innerEvaluator[2] = _cubic_bspline_inner_evaluator2;            
            _innerEvaluator[3] = _cubic_bspline_inner_evaluator3;
            _innerEvaluator[4] = _cubic_bspline_inner_evaluator4;
            _innerEvaluator[5] = _cubic_bspline_inner_evaluator5;            
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(-1,1);
            _innerSupport[1] = Support<T>(-1,1);
            _innerSupport[2] = Support<T>( 0,1);
            _innerSupport[3] = Support<T>( 0,1);
            _innerSupport[4] = Support<T>( 0,1);
            _innerSupport[5] = Support<T>( 0,1);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0] = linspace(-1.0,1.0,9);
            _innerSingularSupport[1] = linspace( 0.0,1.0,9);
            _innerSingularSupport[2] = linspace( 0.0,1.0,3);
            _innerSingularSupport[3] = linspace(-1.0,1.0,3);
            _innerSingularSupport[4] = linspace( 0.0,1.0,5);
            _innerSingularSupport[5] = linspace( 0.0,1.0,5);
            
            // right B-splines
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _cubic_bspline_right_evaluator1;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(-1,0);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0] = linspace(-1.0,0.0,5);
            
            _rightScalingFactors.engine().resize(1,0);
            _rightScalingFactors(0) = 1.0;
            
            break;

        default: std::cerr << "Boundary conditions not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }    
}

} // namespace lawa
