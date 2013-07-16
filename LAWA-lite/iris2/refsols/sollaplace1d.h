#ifndef IRIS2_REFSOLS_LAPLACE1D_H
#define IRIS2_REFSOLS_LAPLACE1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

/*
 * Reference solutions u and corresponding righthand sides for second order PDEs
 * with constant coefficients on R:
 *       - alpha * u'' = f
 */

template<typename T>
struct SolLaplace1D
{
    static int nr;

    static T alpha;

    static DenseVector<Array<T> > sing_pts;

    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas;

    static void
    setExample(int nr, T alpha = T(1));

    static T
    exact(T x, int deriv);

    static T
    u(T x);

    static T
    d_u(T x);

    static T
    rhs(T x);

    static T
    H1norm();

    static
    int
    getMinimalLevel(int d, int d_);

    static
    void
    getRHS_W_XBSplineParameters(int d, int d_, T &_left_bound, T &_right_bound,
                                int &_J_plus_smooth, int &_J_plus_singular,
                                bool &_singular_integral, T eps=1e-5);

    static
    void
    getRHS_WO_XBSplineParameters(int d, int d_, T &_left_bound, T &_right_bound,
                                 int &_J_plus_smooth, int &_J_minus_smooth,
                                 int &_J_plus_singular, int &_J_minus_singular,
                                 bool &_singular_integral, T eps=1e-5);
};


} // namespace lawa

#endif // IRIS2_REFSOLS_LAPLACE1D_H
