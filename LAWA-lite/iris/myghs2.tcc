#ifndef IRIS_MYGHS2_TCC
#define IRIS_MYGHS2_TCC 1

#include <limits>
#include <map>
#include <iris/iris.h>

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
    static RealDenseVector      r, Aw, AtAw, Aw_, AtAw_;
    static IntegerDenseVector   rAbsSorted;
    double                      rNorm;
    double                      rLambdaNormSquare;

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
        //MyApply<Operator, PrecondId>  A(opA, Id, zeta/2000000000000);
        double tol = std::numeric_limits<double>::epsilon();
        MyApply<Operator, PrecondId>  A(opA, Id, tol);

        rhs.filter(zeta/2, r);
        Aw   = A*w;
        AtAw = transpose(A)*Aw;
        r -= AtAw;

        //std::cerr << "Aw = " << Aw << std::endl;
        /*

        Aw_   = opA*w;
        AtAw_ = transpose(opA)*Aw_;
        r    -= AtAw_;

        //std::cerr << "Aw_ = " << Aw_ << std::endl;

        Aw_   -= Aw;
        AtAw_ -= AtAw;
        //std::cerr << "diff(Aw) = " << Aw_ << std::endl;
        std::cerr << "diff(AtAw) = " << AtAw_ << std::endl;
        */

        // std::cerr << "r = " << r << std::endl;

        rNorm = sqrt(r*r);
        nu    = rNorm + zeta;

        /*
        std::cerr << "nu =    " << nu << std::endl;
        std::cerr << "omega = " << omega << std::endl;
        std::cerr << "rNorm = " << rNorm << std::endl;
        std::cerr << "omega*rNorm = " << omega*rNorm << std::endl;
        std::cerr << "zeta =  " << zeta << std::endl;
        */

        if ((nu <= epsilon) || (zeta <= omega*rNorm)) {
            break;
        }

    } while (true);

    std::cerr << "r = " << r << std::endl;
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
        }

    }
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS2<Operator, Rhs, Precond>::galsolve(const IndexSet<int>    &Lambda,
                                         const RealDenseVector  &g,
                                         RealDenseVector        &w,
                                         double                 delta,
                                         double                 epsilon) const
{
    using namespace std;

    const int N = Lambda.size();

    if (N==0) {
        return;
    }

    Precond        P;
    RealGeMatrix   B;
    static int     k = 0;

    myRestrict(opA, P, Lambda, k, B, epsilon);

    std::cerr << "galsolve: k = " << k << std::endl;
    std::cerr.precision(20);
    //std::cerr << "galsolve: B = " << B << std::endl;
    //std::cerr << "galsolve: P = " << P << std::endl;


    //MyApply<Operator, PrecondId>  A(opA, Id, zeta/2000000000000);
    double tol = std::numeric_limits<double>::epsilon();
    MyApply<Operator, PrecondId>  A(opA, Id, tol);

    RealDenseVector Aw, AtAw;

    //std::cerr << "w = " << w << std::endl;

    Aw = A*w;
    //Aw = opA*w;

    //std::cerr << "Aw = " << Aw << std::endl;

    AtAw = transpose(A)*Aw;
    //AtAw = transpose(opA)*Aw;

    //std::cerr << "AtAw = " << AtAw << std::endl;
    //std::cerr << "g = " << g << std::endl;

    RealDenseVector r;
    myRestrict(g, Lambda, r);
    myRestrictSub(AtAw, Lambda, r);

    std::cerr << "galsolve: r = " << r << std::endl;

    RealDenseVector x(B.numCols());
    //std::cerr << "START: lawa::cg(B, x, r, epsilon, 10*N*N);" << std::endl;
    int numIt = lawa::cg(B, x, r);
    //int numIt = lawa::cg(B, x, r, epsilon/3000000000, 10*N*N);
    std::cerr << "numIt = " << numIt << std::endl;
    //std::cerr << "END: lawa::cg(B, x, r, epsilon, 10*N*N);" << std::endl;

    std::cerr << "x = " << x << std::endl;

    myExpandAdd(x, Lambda, w);
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


} // namespace lawa

#endif // IRIS_MYGHS2_H
