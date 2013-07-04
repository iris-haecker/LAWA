#ifndef IRIS_MYGHS2_H
#define IRIS_MYGHS2_H 1

#include <iris/mybasis.h>
#include <iris/myprecondid.h>

namespace lawa {

template <typename Operator, typename Rhs, typename Precond>
struct MyGHS2
{
    typedef flens::GeMatrix<FullStorage<double, ColMajor> > RealGeMatrix;
    typedef flens::DenseVector<Array<double> >              RealDenseVector;
    typedef flens::DenseVector<Array<int> >                 IntegerDenseVector;
    typedef MyPrecondId<double>                             PrecondId;

    MyGHS2(const Operator  &opA,
           const Rhs       &rhs,
           const Precond   &P,
           double          alpha,
           double          omega,
           double          gamma,
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
    solve(double           nuM1,
          double           epsilon,
          int              numOfIterations,
          RealDenseVector  &w) const;

    void
    solve(double           nuM1,
          double           epsilon,
          int              numOfIterations,
          RealDenseVector  &w,
          Function<double> &sol) const;

    const Operator      &opA;
    const Rhs           &rhs;
    const Precond       &P;
    const double        alpha, omega, gamma, theta;
    const PrecondId     Id;
};

} // namespace lawa

#endif // IRIS_MYGHS2_H
