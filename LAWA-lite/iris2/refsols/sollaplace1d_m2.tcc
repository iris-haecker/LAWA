#ifndef IRIS2_REFSOLS_SOLLAPLACE1D_M2_TCC
#define IRIS2_REFSOLS_SOLLAPLACE1D_M2_TCC 1

#include <iris2/refsols/sollaplace1d.h>

namespace lawa {

template <typename T>
int
SolLaplace1D_M2<T>::nr;

template <typename T>
T
SolLaplace1D_M2<T>::alpha;


template <typename T>
flens::DenseVector<Array<T> >
SolLaplace1D_M2<T>::sing_pts;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
SolLaplace1D_M2<T>::deltas;


template <typename T>
void
SolLaplace1D_M2<T>::setExample(int _nr, T _alpha)
{
    nr    = _nr;
    alpha = _alpha;

    assert(nr==1);
}

template <typename T>
T
SolLaplace1D_M2<T>::exact(T x, T y)
{

    if (nr==1) {
        return sin(M_PI*x)*sin(M_PI*y);
    }
    
    assert(0);
    return 0;
}

template <typename T>
T
SolLaplace1D_M2<T>::u(T x, T y)
{
    return SolLaplace1D<T>::exact(x, y);
}

template <typename T>
T
SolLaplace1D_M2<T>::rhs(T x, T y)
{
    if (nr==1) {
        return M_PI*M_PI*sin(M_PI*x)*M_PI*M_PI*sin(M_PI*y);
    }

    assert(0);
    return 0;
}

template <typename T>
T
SolLaplace1D_M2<T>::H1norm()
{
    assert(0);
    return 0;
}

template <typename T>
int
SolLaplace1D_M2<T>::getMinimalLevel(int d, int d_)
{
    return d - 1;
}

template <typename T>
void
SolLaplace1D_M2<T>::getRHS_W_XBSplineParameters(int  d,
                                                int  d_,
                                                T    &_left_bound,
                                                T    &_right_bound,
                                                int  &_J_plus_smooth,
                                                int  &_J_plus_singular,
                                                bool &_singular_integral,
                                                T    /*eps*/)
{
    // TODO:  So far these parameters are hard coded in the examples ...
}

template <typename T>
void
SolLaplace1D_M2<T>::getRHS_WO_XBSplineParameters(int  d,
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
    // TODO:  So far these parameters are hard coded in the examples ...
}

} //namespace lawa

#endif // IRIS2_REFSOLS_SOLLAPLACE1D_M2_TCC
