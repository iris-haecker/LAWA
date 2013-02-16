#ifndef LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BASIS_H
#define LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BASIS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Orthogonal,Interval,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Multi;
        
        typedef BasisFunction<T,Orthogonal,Interval,Multi> BasisFunctionType;
        typedef BSpline<T,Orthogonal,Interval,Multi> BSplineType;
        typedef Wavelet<T,Orthogonal,Interval,Multi> WaveletType;

        Basis(const int d, const int j=-1);
    
        virtual
        ~Basis();
    
        int
        level() const;
    
        void
        setLevel(const int j) const;
    
        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();
    
        const BasisFunctionType &
        generator(XType xtype) const;

        //--- cardinalities of whole, left, inner, right index set.
        long
        cardJ(const int j) const;
        
        long
        cardJL(const int j=-1) const;

        long
        cardJI(const int j) const;

        long
        cardJR(const int j=-1) const;
    
        //--- ranges of whole, left, inner, right index set.
        const flens::Range<int>
        rangeJ(const int j) const;
    
        const flens::Range<int>
        rangeJL(const int j=-1) const;

        const flens::Range<int>
        rangeJI(const int j) const;

        const flens::Range<int>
        rangeJR(const int j=-1) const;
    
        MRA<T,Orthogonal,Interval,Multi> mra;
    
        const int d;
        const int j0;          // minimal used(!) level.
    
    private:
        DenseVector<Array<int> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
                                       // bc(1) = 1 -> Dirichlet BC right.
        
        mutable int _j;                // the current level.
    
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        friend class Wavelet<T,Orthogonal,Interval,Multi>;

        unsigned int _numLeftParts, 
                     _numInnerParts, 
                     _numRightParts;
        Evaluator *_leftEvaluator, 
                  *_innerEvaluator, 
                  *_rightEvaluator;
        Support<T> *_leftSupport, 
                   *_innerSupport, 
                   *_rightSupport;
        DenseVector<Array<T> > *_leftSingularSupport, 
                               *_innerSingularSupport, 
                               *_rightSingularSupport;
        DenseVector<Array<T> > _leftScalingFactors, _rightScalingFactors;
        
    public:
        Wavelet<T,Orthogonal,Interval,Multi> psi;
};

} // namespace lawa

#include <lawa/constructions/multi/interval/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BASIS_H
