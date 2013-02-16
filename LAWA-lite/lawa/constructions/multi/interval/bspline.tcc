#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
BSpline<T,Orthogonal,Interval,Multi>::BSpline(const MRA<T,Orthogonal,Interval,Multi> &_mra)
    : mra(_mra), d(_mra.d)
{
}
    
template <typename T>
BSpline<T,Orthogonal,Interval,Multi>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2ih<T>(2*j*deriv+j) * mra._leftScalingFactors((int)k) *
               mra._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }
    
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil((k+1.-mra._numLeftParts)/mra._numInnerParts);
        return pow2ih<T>(2*j*deriv+j) * 
               mra._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
    
    // right part
    int type  = (int)(k+1 - (mra.cardI(j) - mra._numRightParts + 1));
    long shift = iceil((k+1. - mra._numLeftParts)/mra._numInnerParts);
    return pow2ih<T>(2*j*deriv+j) * mra._rightScalingFactors(type) *
           mra._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}
    
template <typename T>
Support<T>
BSpline<T,Orthogonal,Interval,Multi>::support(int j, long k) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSupport[k];
    }
    
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil((k+1.-mra._numLeftParts)/mra._numInnerParts);
        return pow2i<T>(-j) * (mra._innerSupport[type]+shift);
    }
    
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    long shift = iceil((k+1.-mra._numLeftParts)/mra._numInnerParts);
    return pow2i<T>(-j) * (mra._rightSupport[type]+shift);
}

template <typename T>
DenseVector<Array<T> >
BSpline<T,Orthogonal,Interval,Multi>::singularSupport(int j, long k) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSingularSupport[k];
    }
    
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil((k+1.-mra._numLeftParts)/mra._numInnerParts);
        DenseVector<Array<T> > result = mra._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    long shift = iceil((k+1. - mra._numLeftParts)/mra._numInnerParts);
    DenseVector<Array<T> > result = mra._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}
        
} // namespace lawa
