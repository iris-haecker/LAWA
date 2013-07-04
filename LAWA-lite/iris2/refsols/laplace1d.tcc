#ifndef IRIS2_REFSOLS_LAPLACE1D_TCC
#define IRIS2_REFSOLS_LAPLACE1D_TCC 1

#include <iris2/refsols/laplace.h>

namespace lawa {

template <typename T>
int
Laplace1D<T>::nr;

template <typename T>
T
Laplace1D<T>::alpha;


template <typename T>
flens::DenseVector<Array<T> >
Laplace1D<T>::sing_pts;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
Laplace1D<T>::deltas;


template <typename T>
void
Laplace1D<T>::setExample(int _nr, T _alpha)
{
    nr    = _nr;
    alpha = _alpha;

    assert(nr>=1);
    assert(nr<=1);
    if (nr==1) {

    }
}

template <typename T>
T
Laplace1D<T>::exact(T x, int deriv)
{

    if (nr==1) {
        if (deriv==0) {
            return sin(M_PI*x);
        } else if (deriv==1) {
            return M_PI*cos(M_PI*x);
        } else if (deriv==1) {
            return -M_PI*M_PI*sin(M_PI*x);
        }
    }
    assert(0);
    return 0;
}

template <typename T>
T
Laplace1D<T>::u(T x)
{
    return Laplace1D<T>::exact(x, 0);
}

template <typename T>
T
Laplace1D<T>::d_u(T x)
{
    return Laplace1D<T>::exact(x, 1);
}

template <typename T>
T
Laplace1D<T>::rhs(T x)
{
    return -alpha*exact(x, 2);
}

template <typename T>
T
Laplace1D<T>::H1norm()
{
    if (nr==1) {
        return sqrt((M_PI*M_PI+1)/2);
    }
    assert(0);
    return 0;
}

template <typename T>
int
Laplace1D<T>::getMinimalLevel(int d, int d_)
{
    return d - 1;
}

template <typename T>
void
Laplace1D<T>::getRHS_W_XBSplineParameters(int  d,
                                          int  d_,
                                          T    &_left_bound,
                                          T    &_right_bound,
                                          int  &_J_plus_smooth,
                                          int  &_J_plus_singular,
                                          bool &_singular_integral,
                                          T    /*eps*/)
{
}

template <typename T>
void
Laplace1D<T>::getRHS_WO_XBSplineParameters(int  d,
                                           int  d_,
                                           T    &_left_bound,
                                           T    &_right_bound,
                                           int  &_J_plus_smooth,
                                           int  &_J_minus_smooth,
                                           int  &_J_plus_singular,
                                           int  &_J_minus_singular,
                                           bool &_singular_integral,
                                           T    /*eps*/)
{
}

} //namespace lawa

#endif // IRIS2_REFSOLS_LAPLACE1D_TCC
