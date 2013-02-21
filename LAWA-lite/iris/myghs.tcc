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
                                     double          _theta)
    : opA(_opA), rhs(_rhs), P(_P), alpha(_alpha), omega(_omega), theta(_theta)
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

    double zeta = 2*omega*nu_/(1-omega);

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

} // namespace lawa

#endif // IRIS_MYGHS_H