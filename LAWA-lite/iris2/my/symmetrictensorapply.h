#ifndef IRIS2_MY_SYMMETRICTENSORAPPLY_H
#define IRIS2_MY_SYMMETRICTENSORAPPLY_H 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T, typename Apply1, typename Apply2>
class SymmetricTensorApply
{
    typedef Coefficients<AbsoluteValue,T,Index2D>       CoefficientsAbs;
    typedef Coefficients<Lexicographical,T,Index2D>     CoefficientsLex;

    typedef typename CoefficientsAbs::const_iterator    CoeffAbsIt;
    typedef typename IndexSet<Index2D>::const_iterator  IndexSetIt;

    public:

        Apply1     &A1;
        Apply2     &A2;

        SymmetricTensorApply(Apply1 &A1, Apply2 &A2);

        CoefficientsLex
        operator()(const CoefficientsLex &v, T eps);
};

}   //namespace lawa

#endif // IRIS2_MY_SYMMETRICTENSORAPPLY_H
