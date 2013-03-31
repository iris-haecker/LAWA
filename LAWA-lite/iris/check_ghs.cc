#include <iostream>
//#define SOLVER_DEBUG


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
    const int l  = 1;
    const int u0 = 0;

    //return sin(2*M_PI*l*t) + 2*M_PI*l*cos(2*M_PI*l*t) + u0;
    //return 4*M_PI*M_PI*sin(M_PI*t);
    return 1;
}

int
main()
{
    typedef flens::DenseVector<Array<double> >              RealDenseVector;
    typedef flens::DenseVector<Array<int> >                 IntDenseVector;
    typedef MyOperator<double>                              Operator;
    typedef MyPrecond<Operator>                             Precond;
    typedef MyRhs<double, Precond>                          Rhs;
    typedef MyGHS<Operator, Rhs, Precond>                   GHS;
    typedef MyEval<double>                                  Eval;

    const int d  = 3;
    const int d_ = 5;

    const int    jMax = 14;
    const double eps = 0.0000001;
    const int    numOfIterations = 40;

    Operator              operatorA(d, d_, jMax, jMax);
    Precond               P(operatorA);
    Rhs                   rhs(Function<double>(g), operatorA, P);
    RealDenseVector       w;


    double  alpha = 0.4;
    double  omega = 0.012618;
    double  gamma = 0.009581;
    double  theta = 2./7;

    alpha = 0.9;

    GHS     ghs(operatorA, rhs, P, alpha, omega, gamma, theta);

    ghs.solve(rhs.norm, eps, numOfIterations, w);

    Eval sol(operatorA.U, w);

    sol.dump(1000, "ghs.dat");


    std::fstream   out("ghs_sol_raw.dat", std::fstream::out);
    for (int k=w.firstIndex(); k<=w.lastIndex(); ++k) {
        out << w(k) << " ";
    }
}
