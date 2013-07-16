#ifndef IRIS2_MY_PARAMETERSLAPLACE1D_H
#define IRIS2_MY_PARAMETERSLAPLACE1D_H 1

#include <iris2/my/laplace1d.h>

namespace lawa {

template <typename T>
struct ParametersLaplace1D
{
    const Laplace1D<T>      &A;

    int j0;
    T cA, CA, kappa;
    T alpha, omega, gamma, theta;    //GHSADWAV parameters (from [GHS:2007])

    ParametersLaplace1D(const Laplace1D<T> &A);

    void
    getGHSADWAVParameters(T &_alpha, T &_omega, T &_gamma, T &_theta) const;

};

} // namespace lawa

#endif // IRIS2_MY_PARAMETERSLAPLACE1D_H