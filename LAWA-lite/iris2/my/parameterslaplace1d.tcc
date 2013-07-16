#ifndef IRIS2_MY_PARAMETERSLAPLACE1D_TCC
#define IRIS2_MY_PARAMETERSLAPLACE1D_TCC 1

#include <iris2/iris2.h>

namespace lawa {

template <typename T>
ParametersLaplace1D<T>::ParametersLaplace1D(const Laplace1D<T> &_A)
    : A(_A)
{
    j0 = A.j0;

    cA = 0.58;
    CA = 1.86;

    kappa = CA/cA;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    omega = 0.01;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;
}

template <typename T>
void
ParametersLaplace1D<T>::getGHSADWAVParameters(T &_alpha,
                                              T &_omega,
                                              T &_gamma,
                                              T &_theta) const
{
    _alpha = alpha;
    _omega = omega;
    _gamma = gamma;
    _theta = theta;
}

} // namespace lawa

#endif // IRIS2_MY_PARAMETERSLAPLACE1D_TCC