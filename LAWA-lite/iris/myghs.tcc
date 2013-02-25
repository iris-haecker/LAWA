#ifndef IRIS_MYGHS_TCC
#define IRIS_MYGHS_TCC 1

#include <map>
#include <iris/iris.h>

namespace lawa {

template <typename Operator, typename Rhs, typename Precond>
MyGHS<Operator, Rhs, Precond>::MyGHS(const Operator  &_opA,
                                     const Rhs       &_rhs,
                                     const Precond   &_P,
                                     double          _alpha,
                                     double          _omega,
                                     double          _gamma,
                                     double          _theta)
    : opA(_opA), rhs(_rhs), P(_P), alpha(_alpha), omega(_omega), gamma(_gamma),
      theta(_theta)
{
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS<Operator, Rhs, Precond>::grow(const RealDenseVector  &w,
                                    double                 nu_,
                                    double                 epsilon,
                                    double                 &nu,
                                    IndexSet<int>          &Lambda) const
{
    static RealDenseVector      r;
    static IntegerDenseVector   rAbsSorted;
    double rNorm, rLambdaNormSquare;

    std::cerr << "nu_ =   " << nu_ << std::endl;


    double zeta = 2*omega*nu_/(1-omega);
    std::cerr << "zeta =  " << zeta << std::endl;
    std::cerr << "nu_ =   " << nu_ << std::endl;
    std::cerr << "omega = " << omega << std::endl;


    do {
        zeta = zeta/2;
//
//      r = RHS(zeta/2) - APPLY(w, zeta/2)
//
        MyApply<Operator, Precond>  A(opA, P, zeta/2);
        rhs.filter(zeta/2, r);
        r -= A*w;

        rNorm = sqrt(r*r);
        nu    = rNorm + zeta;

        if ((nu <= epsilon) || (zeta <= omega*rNorm)) {
            break;
        }
        std::cerr << "nu =    " << nu << std::endl;
        std::cerr << "zeta =  " << zeta << std::endl;
        std::cerr << "omega = " << omega << std::endl;
        std::cerr << "rNorm = " << rNorm << std::endl;
        std::cerr << "omega*rNorm = " << omega*rNorm << std::endl;

    } while (true);

    std::cerr << "rNorm = " << rNorm << std::endl;
    std::cerr << "r = " << r << std::endl;

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
        if (rLambdaNormSquare<pow(alpha*rNorm, 2)) {

//
//          We should use a bucket sort here
//
            myAbsSort(r, rAbsSorted);
            std::cerr << "rAbsSorted = " << rAbsSorted << std::endl;

            for (int k=1; k<=rAbsSorted.length(); ++k) {

                std::cerr << "rLambdaNormSquare = "
                          << rLambdaNormSquare << std::endl;
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
                rLambdaNormSquare += pow(r(rAbsSorted(k)), 2);
                if (rLambdaNormSquare>=pow(alpha*rNorm,2)) {
                    break;
                }
            }

        }

    }
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS<Operator, Rhs, Precond>::galsolve(const IndexSet<int>    &Lambda,
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

    cerr << "1" << endl;

    SparseGeMatrix<CRS<double ,CRS_General> >  B;
    myRestrict(opA, P, Lambda, B);

    MyApply<Operator, Precond>  A(opA, P, epsilon/3);
    RealDenseVector Aw;
    cerr << "2" << endl;

    Aw = A*w;
    cerr << "3" << endl;

    RealDenseVector r;
    myRestrict(g, Lambda, r);
    cerr << "4" << endl;
    myRestrictSub(Aw, Lambda, r);
    cerr << "5" << endl;

    RealDenseVector x(B.numCols());
    lawa::cg(B, x, r, epsilon, 4);
    std::cerr << "x =" << x << std::endl;

    myExpandAdd(x, Lambda, w);
}

template <typename Operator, typename Rhs, typename Precond>
void
MyGHS<Operator, Rhs, Precond>::solve(double           nuM1,
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
        if (nu_k<=epsilon) {
            break;
        }

        rhs.filter(0, g_kP1);

        galsolve(Lambda_kP1, g_kP1, w, (1+gamma)*nu_k, gamma*nu_k);
        nu_kM1 = nu_k;
    }
}


} // namespace lawa

#endif // IRIS_MYGHS_H