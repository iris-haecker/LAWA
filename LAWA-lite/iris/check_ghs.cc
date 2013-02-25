#include <iostream>
#define SOLVER_DEBUG


#include <lawa/flensforlawa.h>
#include <extensions/flens/sparsematrix.h>
#include <extensions/flens/crs.h>

namespace flens {

// sparse_gemv
template <typename T, typename VX, typename VY>
void
mv(cxxblas::Transpose trans, T alpha, const SparseGeMatrix<CRS<T> > &A,
   const DenseVector<VX> &x,
   typename DenseVector<VY>::ElementType beta, DenseVector<VY> &y);

}


#include <iris/iris.cxx>


using namespace flens;
using namespace lawa;
using namespace std;

double
g(double t)
{
    const int l  = 100;
    const int u0 = 0;

    return sin(2*M_PI*l*t) + 2*M_PI*l*cos(2*M_PI*l*t) + u0;
}

int
main()
{
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::DenseVector<Array<int> >                  IntDenseVector;
    typedef MyOperator<double>                               Operator;
    typedef MyPrecond<Operator>                              Precond;
    typedef MyRhs<double, Precond>                           Rhs;
    typedef MyGHS<Operator, Rhs, Precond>                    GHS;

    const int d  = 2;
    const int d_ = 4;

    const int    jMax = 5;
    const double eps = 0.00001;
    const int    numOfIterations = 20;

    Operator              operatorA(d, d_, jMax);
    Precond               P(operatorA);
    Rhs                   rhs(Function<double>(g), operatorA, P);
    RealDenseVector       w;

    GHS                   ghs(operatorA, rhs, P, 0.4, 0.012618, 0.009581, 2./7);

    ghs.solve(rhs.norm, eps, numOfIterations, w);

}