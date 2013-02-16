#ifndef LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BSPLINE_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
class BSpline<_T,Orthogonal,Interval,Multi>
    : public BasisFunction<_T,Orthogonal,Interval,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Multi;

        BSpline(const MRA<T,Orthogonal,Interval,Multi> &mra);

        virtual
        ~BSpline();

        T
        operator()(T x, int j, long k, unsigned short deriv) const;

        Support<T>
        support(int j, long k) const;

        DenseVector<Array<T> >
        singularSupport(int j, long k) const;

        T
        tic(int j) const;

        const MRA<T,Orthogonal,Interval,Multi> &mra;
        const unsigned int d;
};

} // namespace lawa

#include <lawa/constructions/multi/interval/bspline.tcc>

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BSPLINE_H
