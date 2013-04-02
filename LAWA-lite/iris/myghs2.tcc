#ifndef IRIS_MYGHS2_TCC
#define IRIS_MYGHS2_TCC 1

#include <limits>
#include <map>
#include <iris/iris.h>
#include <cmath>
#include <sstream>

namespace lawa {

template <typename Operator, typename Rhs, typename Precond>
MyGHS2<Operator, Rhs, Precond>::MyGHS2(const Operator  &_opA,
                                       const Rhs       &_rhs,
                                       double          _alpha,
                                       double          _omega,
                                       double          _gamma,
                                       double          _theta)
    : opA(_opA), rhs(_rhs), alpha(_alpha), omega(_omega), gamma(_gamma),
      theta(_theta)
{
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS2<Operator, Rhs, Precond>::grow(const RealDenseVector  &w,
                                     double                 nu_,
                                     double                 epsilon,
                                     double                 &nu,
                                     IndexSet<int>          &Lambda) const
{
    std::cerr << "BEGIN: grow" << std::endl;

    using std::sqrt;

    static RealDenseVector      r, Aw, AtAw;
    static IntegerDenseVector   rAbsSorted;
    double                      rNorm = 0;
    double                      rLambdaNormSquare = 0;

    double zeta = 2*omega*nu_/(1-omega);
    /*
    std::cerr << "zeta =  " << zeta << std::endl;
    std::cerr << "nu_ =   " << nu_ << std::endl;
    std::cerr << "omega = " << omega << std::endl;
    */

    do {
        zeta = zeta/2;
//
//      r = RHS(zeta/2) - APPLY(w, zeta/2)
//
        rhs.filter(zeta/2, r);

#       ifdef  USE_APPLY
        MyApply<Operator, PrecondId>  A(opA, Id, zeta/2);

        Aw   = A*w;
        AtAw = transpose(A)*Aw;
        r -= AtAw;
#       else

        Aw   = opA*w;
        AtAw = Aw; // !!! transpose(opA)*Aw;
        r    -= AtAw;
#       endif

        rNorm = sqrt(r*r);
        nu    = rNorm + zeta;

        std::cerr << "rNorm =       " << rNorm << std::endl;
        std::cerr << "zeta =        " << zeta << std::endl;
        std::cerr << "nu =          " << nu << std::endl;
        std::cerr << "omega =       " << omega << std::endl;
        std::cerr << "omega*rNorm = " << omega*rNorm << std::endl;

        if ((nu <= epsilon) || (zeta <= omega*rNorm)) {
            break;
        }

    } while (true);

    //std::cerr << "r = " << r << std::endl;
    std::cerr << "rNorm = " << rNorm << std::endl;
    std::cerr << "nu = " << nu << std::endl;

    if (nu>epsilon) {
        mySupport(w, Lambda);

        rLambdaNormSquare = 0;

//
//      Compute || P_Lambda r||^2
//
        IndexSet<int>::const_iterator it;
        for (it=Lambda.begin(); it!=Lambda.end(); ++it) {
            rLambdaNormSquare += pow(r(*it), 2);
        }

        std::cerr << "residual sqr-norm of r on old set: rLambdaNormSquare ="
                  << rLambdaNormSquare
                  << std::endl;

//
//      If necessary: Increase index set Lambda until
//                    || P_Lambda r||^2  >=  alpha ||r||^2
//
        if (rLambdaNormSquare<alpha*pow(rNorm, 2)) {

//
//          We should use a bucket sort here
//
            //std::cerr << "START: myAbsSort(r, rAbsSorted);" << std::endl;
            myAbsSort(r, rAbsSorted);
            //std::cerr << "END: myAbsSort(r, rAbsSorted);" << std::endl;
            // std::cerr << "rAbsSorted = " << rAbsSorted << std::endl;

            for (int k=1; k<=rAbsSorted.length(); ++k) {

//
//              Check if index is already in Lambda
//
                if (Lambda.count(rAbsSorted(k))>0) {
                    continue;
                }
//
//              Update Lambda and || P_Lambda r||^2
//
                Lambda.insert(rAbsSorted(k));
                // std::cerr << "insert: " << rAbsSorted(k) << std::endl;

                rLambdaNormSquare += pow(r(rAbsSorted(k)), 2);
                if (rLambdaNormSquare>=alpha*pow(rNorm,2)) {
                    break;
                }
            }
            // std::cerr << "Lambda.size() = " << Lambda.size() << std::endl;
            std::cerr << "residual sqr-norm of r on new set: rLambdaNormSquare ="
                      << rLambdaNormSquare
                      << std::endl;
        }

    }
    std::cerr << "END: grow" << std::endl;
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS2<Operator, Rhs, Precond>::galsolve(const IndexSet<int>    &Lambda,
                                         const RealDenseVector  &g,
                                         RealDenseVector        &w,
                                         double                 delta,
                                         double                 epsilon) const
{
    std::cerr << "BEGIN: galsolve" << std::endl;

    using namespace std;


    const int N = Lambda.size();

    if (N==0) {
        return;
    }

    RealDenseVector  Aw, AtAw, r, x;
    Precond          P;

    SparseGeMatrix<CRS<double, CRS_General> >  B;
    //RealGeMatrix     B;
    static int       k = 0;


    myRestrict(opA, P, Lambda, B);
    // myRestrict(opA, P, Lambda, k, B, epsilon);
    std::cerr << "galsolve: k = " << k << std::endl;

#   ifdef  USE_APPLY

    MyApply<Operator, PrecondId>  A(opA, Id, zeta/3);


    Aw = A*w;
    AtAw = transpose(A)*Aw;

#   else

    Aw   = opA*w;
    AtAw = Aw; // !!! transpose(opA)*Aw;

#   endif

    myRestrict(g, Lambda, r);
    myRestrictSub(AtAw, Lambda, r);

    //std::cerr << "galsolve: r = " << r << std::endl;

    x.engine().resize(B.numCols());
    int numIt = lawa::cg(B, x, r);

    myExpandAdd(x, Lambda, w);

    std::cerr << "numIt = " << numIt << std::endl;
    //std::cerr << "x = " << x << std::endl;

    std::cerr << "END: galsolve" << std::endl;
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS2<Operator, Rhs, Precond>::solve(double           nuM1,
                                      double           epsilon,
                                      int              numOfIterations,
                                      RealDenseVector  &w) const
{
    double           nu_kM1 = nuM1;
    double           nu_k;
    RealDenseVector  g_kP1;
    IndexSet<int>    Lambda_kP1;

    if (w.length()!=opA.numCols()) {
        w.engine().resize(opA.numCols());
    }

    for (int k=0; k<numOfIterations; ++k) {
        grow(w, theta*nu_kM1, epsilon, nu_k, Lambda_kP1);

        std::cerr << "k = " << k
                  << ": nu_k = " << nu_k
                  << std::endl;
        std::cerr << "Lambda_kP1 = "
                  << Lambda_kP1
                  << std::endl;

        if (nu_k<=epsilon) {
            break;
        }

        rhs.filter(0, g_kP1);

        galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);
        nu_kM1 = nu_k;

        // std::cerr << "w = " << w << std::endl;

    }
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS2<Operator, Rhs, Precond>::solve(double           nuM1,
                                      double           epsilon,
                                      int              numOfIterations,
                                      RealDenseVector  &w,
                                      Function<double> &sol) const
{
    using std::sqrt;

    double             nu_kM1 = nuM1;
    double             nu_k;
    RealDenseVector    g_kP1;
    IndexSet<int>      Lambda_kP1;

    if (w.length()!=opA.numCols()) {
        w.engine().resize(opA.numCols());
    }

    for (int k=0; k<numOfIterations; ++k) {
        grow(w, theta*nu_kM1, epsilon, nu_k, Lambda_kP1);

        std::cerr << "Lambda_kP1 = "
                  << Lambda_kP1
                  << std::endl;

        MyEval<double> eval(opA.U, w);

        double errorH1 = eval.error_HNorm(opA, rhs, w, sqrt(330.0)/60.0);

        std::cerr << "GHS: " << std::endl
                  << "    k =          " << k  << std::endl
                  << "    N =          " << Lambda_kP1.size() << std::endl
                  << "    nu_k =       " << nu_k << std::endl
                  << "    Error-L1 =   " << eval.diff_L1(1000, sol) << std::endl
                  << "    Error-L2 =   " << eval.diff_L2(1000, sol) << std::endl
                  << "    Error-LInf = " << eval.diff_LInf(1000, sol) << std::endl
                  << "    Error-H1   = " << errorH1 << std::endl
                  << std::endl;

        std::stringstream  filename;

        filename << "snap/snapshot_" << k << ".dat";
        eval.dump(1000, sol, filename.str().c_str());

        if (nu_k<=epsilon) {
            break;
        }

        rhs.filter(0, g_kP1);

        galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);
        nu_kM1 = nu_k;

        // std::cerr << "w = " << w << std::endl;

    }

    MyEval<double> eval(opA.U, w);

    

}



} // namespace lawa

#endif // IRIS_MYGHS2_H
