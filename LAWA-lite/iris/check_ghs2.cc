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

    return 1;
    //return sin(2*M_PI*l*t);

    // return sin(2*M_PI*l*t) + 2*M_PI*l*cos(2*M_PI*l*t) + u0;
}

double
sol_(double x)
{
    return 0.5*x*(1-x);
}


int
main()
{
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::DenseVector<Array<int> >                  IntDenseVector;
    typedef MyOperator<double>                               Operator;
    typedef MyPrecond3<double>                               Precond;
    typedef MyPrecondId<double>                              PrecondId;
    typedef MyRhs2<double, PrecondId>                        Rhs;
    typedef MyGHS2<Operator, Rhs, Precond>                   GHS;
    typedef MyEval<double>                                   Eval;

    const int d  = 2;
    const int d_ = 4;

    const int    jMaxV = 9;
    const int    jMaxU = 9;
    const double eps = 0.000000001;
    const int    numOfIterations = 400;

    Operator              operatorA(d, d_, jMaxV, jMaxU);
    Precond               P;
    PrecondId             Id;
    Rhs                   rhs(Function<double>(g), operatorA, Id, eps);
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

    GHS     ghs(operatorA, rhs, alpha, omega, gamma, theta);

    ghs.solve(rhs.norm, eps, numOfIterations, w, sol);



/*
    double           nu_kM1 = rhs.norm;
    double           nu_k;
    IndexSet<int>    Lambda_kP1;
    RealGeMatrix     B;
    RealDenseVector  g_kP1;

    rhs.filter(0, g_kP1);

//
//  GROW & GALSOLVE
//
    std::cerr << std::endl
              << "-------------------------------" << std::endl
              << "GROW & GALSOLVE" << std::endl
              << "-------------------------------" << std::endl;

    ghs.grow(w, theta*nu_kM1, eps, nu_k, Lambda_kP1);

    std::cerr << "Lambda_kP1 = " << Lambda_kP1 << std::endl;

    ghs.galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);

    std::cerr << "w = " << w << std::endl;


//
//  GROW & GALSOLVE
//
    std::cerr << std::endl
              << "-------------------------------" << std::endl
              << "GROW & GALSOLVE" << std::endl
              << "-------------------------------" << std::endl;

    ghs.grow(w, theta*nu_kM1, eps, nu_k, Lambda_kP1);

    std::cerr << "Lambda_kP1 = " << Lambda_kP1 << std::endl;

    ghs.galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);

    std::cerr << "w = " << w << std::endl;

//
//  GROW & GALSOLVE
//
    std::cerr << std::endl
              << "-------------------------------" << std::endl
              << "GROW & GALSOLVE" << std::endl
              << "-------------------------------" << std::endl;

    ghs.grow(w, theta*nu_kM1, eps, nu_k, Lambda_kP1);

    std::cerr << "Lambda_kP1 = " << Lambda_kP1 << std::endl;

    ghs.galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);

    std::cerr << "w = " << w << std::endl;

//
//  GROW & GALSOLVE
//
    std::cerr << std::endl
              << "-------------------------------" << std::endl
              << "GROW & GALSOLVE" << std::endl
              << "-------------------------------" << std::endl;

    ghs.grow(w, theta*nu_kM1, eps, nu_k, Lambda_kP1);

    std::cerr << "Lambda_kP1 = " << Lambda_kP1 << std::endl;

    ghs.galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);

    std::cerr << "w = " << w << std::endl;

//
//  GROW & GALSOLVE
//
    std::cerr << std::endl
              << "-------------------------------" << std::endl
              << "GROW & GALSOLVE" << std::endl
              << "-------------------------------" << std::endl;

    ghs.grow(w, theta*nu_kM1, eps, nu_k, Lambda_kP1);

    std::cerr << "Lambda_kP1 = " << Lambda_kP1 << std::endl;

    ghs.galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);

    std::cerr << "w = " << w << std::endl;

*/


/*

    ghs.galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);

    std::cerr << "2) after galsolve:" << std::endl;
    std::cerr << "w = " << w << std::endl;

    ghs.grow(w, theta*nu_kM1, eps, nu_k, Lambda_kP1);

    std::cerr << "Lambda_kP1 = " << Lambda_kP1 << std::endl;
*/


}
