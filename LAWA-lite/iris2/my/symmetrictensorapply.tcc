#ifndef IRIS2_MY_SYMMETRICTENSORAPPLY_TCC
#define IRIS2_MY_SYMMETRICTENSORAPPLY_TCC 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T, typename Apply1, typename Apply2>
SymmetricTensorApply<T,Apply1,Apply2>::SymmetricTensorApply(Apply1 &_A1,
                                                            Apply2 &_A2)
    : A1(_A1), A2(_A2)
{
}

template <typename T, typename Apply1, typename Apply2>
typename SymmetricTensorApply<T,Apply1,Apply2>::CoefficientsLex
SymmetricTensorApply<T,Apply1,Apply2>::operator()(const CoefficientsLex  &u,
                                                  T                      eps)
{
    CoefficientsLex v;
    CoefficientsLex w;

    //
    // Compute V = A2*U with v=vec(V), u=vec(U)
    //
    auto U = splitCols(u);
    for (auto itU=U.begin(); itU!=U.end(); ++itU) {
        Index1D col = itU->first;
        auto    _u  = itU->second;

        auto _v = A2(_u, eps);
        join(v, _v, col);
    }

    //
    // Compute W' = A1*V' with w=vec(W), v=vec(V)
    //
    auto V = splitRows(v);
    for (auto itV=V.begin(); itV!=V.end(); ++itV) {
        Index1D row = itV->first;
        auto    _v  = itV->second;
        
        auto _w = A1(_v, eps);
        join(w, row, _w);
    }

    return w;
}

} // namespace lawa

#endif // IRIS2_MY_SYMMETRICTENSORAPPLY_TCC
