#ifndef IRIS2_MY_SYMMETRICAPPLY1D_H
#define IRIS2_MY_SYMMETRICAPPLY1D_H 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T, typename MA>
class SymmetricApply1D
{

    typedef Coefficients<AbsoluteValue,T,Index1D>       CoefficientsAbs;
    typedef Coefficients<Lexicographical,T,Index1D>     CoefficientsLex;

    typedef typename CoefficientsAbs::const_iterator    CoeffAbsIt;
    typedef typename IndexSet<Index1D>::const_iterator  IndexSetIt;

    public:
        typedef MA MAType;

        const ParametersLaplace1D<T>    &parameters;
        MA                              &A;

        SymmetricApply1D(const ParametersLaplace1D<T>  &parameters,
                         MA                            &A);

        CoefficientsLex
        operator()(const CoefficientsLex &v, int k);

        CoefficientsLex
        operator()(const CoefficientsLex &v, T eps);

    private:

        int
        findK(const CoefficientsAbs &v, T eps);

};

}   //namespace lawa

#endif // IRIS2_MY_SYMMETRICAPPLY1D_H
