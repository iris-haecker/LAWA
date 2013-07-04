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
    using std::pow;

    const int l  = 8;
    const int u0 = 0;

    /*
    if (t<1.0/3) {
        return 1;
    }
    return 20;
    */

    //return pow(2*M_PI*l, 2) * sin(2*M_PI*l*t);

    return 1;
    //return sin(2*M_PI*l*t);

    // return sin(2*M_PI*l*t) + 2*M_PI*l*cos(2*M_PI*l*t) + u0;
}

double
sol_(double x)
{
    const int l  = 8;

    // return sin(2*M_PI*l*x);
    return 0.5*x*(1-x);
}


int
main()
{
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::DenseVector<Array<int> >                  IntDenseVector;
    typedef MyOperator<double>                               Operator;
    typedef MyPrecond<Operator>                              Precond;
    typedef MyPrecondId<double>                              PrecondId;
    typedef MyRhs2<double, Precond>                          Rhs;
    typedef MyGHS2<Operator, Rhs, Precond>                   GHS;
    typedef MyEval<double>                                   Eval;

    const int d  = 3;
    const int d_ = 5;

    const int    jMaxV = 9;
    const int    jMaxU = 9;
    const double eps = 0.000000001;
    const int    numOfIterations = 25;

    Operator              operatorA(d, d_, jMaxV, jMaxU);
    Precond               P(operatorA);
    PrecondId             Id;
    //Rhs                   rhs(Function<double>(g), operatorA, Id, eps);
    Rhs                   rhs(Function<double>(g), operatorA, P, eps);
    Function<double>      sol(sol_);

    RealDenseVector       w;

    w.engine().resize(operatorA.numCols());

    double  alpha = 0.5;
    double  omega = 0.001;
    double  gamma = 0.01;
    double  theta = 0.28;

    std::cerr.precision(20);
    std::cerr << "rhs = " << rhs.rhsData << std::endl;
    std::cerr.precision(20);
    //std::cerr << "operatorA = " << operatorA << std::endl;

    GHS  ghs(operatorA, rhs, P, alpha, omega, gamma, theta);

    ghs.solve(rhs.norm, eps, numOfIterations, w, sol);

}
