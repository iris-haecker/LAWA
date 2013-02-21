#ifndef IRIS_MYGHS_H
#define IRIS_MYGHS_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename Operator, typename Rhs, typename Precond>
struct MyGHS
{
    typedef flens::DenseVector<Array<double> >      RealDenseVector;
    typedef flens::DenseVector<Array<int> >         IntegerDenseVector;

    MyGHS(const Operator  &opA,
          const Rhs       &rhs,
          const Precond   &P,
          double          alpha,
          double          omega,
          double          theta);

    void
    grow(const RealDenseVector  &w,
         double                 nu_,
         double                 epsilon,
         double                 &nu,
         IndexSet<int>          &Lambda) const;

     void
     galsolve(const IndexSet<int>    &Lambda,
              const RealDenseVector  &g,
              RealDenseVector        &w,
              double                 delta,
              double                 epsilon) const;



    void
    solve() const;

    const Operator  &opA;
    const Rhs       &rhs;
    const Precond   &P;
    const double    alpha, omega, theta;
};

} // namespace lawa

#endif // IRIS_MYGHS_H