#ifndef LAWA_CONSTRUCTIONS_MULTI_REALLINE_BSPLINE_H
#define LAWA_CONSTRUCTIONS_MULTI_REALLINE_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Orthogonal,R,Multi>
    : public BasisFunction<_T,Orthogonal,R,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = Multi;
        
        BSpline(const int _d);

        //TODO    BSpline(MRA<T,Orthogonal,R,Multi> &mra); 
        
        virtual
        ~BSpline();
        
        T
        operator()(T x, int j, long k, unsigned short deriv) const;
        
        Support<T>
        support(int j, long k) const;
        
        DenseVector<Array<T> >
        singularSupport(int j, long k) const;

        const unsigned int d;

    private:
        typedef T (*Evaluator)(T x, unsigned short deriv);

        long
        _shift(long k) const;

        int
        _type(long k) const;

        unsigned int _numSplines;
        Evaluator *_evaluator;
        Support<T> *_support;
        DenseVector<Array<T> > *_singularSupport;
    //        T
    //TODO    tic(int j) const;
    //    int polynomialOrder;
};

} // namespace lawa

#include <lawa/constructions/multi/realline/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_MULTI_REALLINE_BSPLINE_H
