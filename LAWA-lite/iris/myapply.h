#ifndef IRIS_MYAPPLY_H
#define IRIS_MYAPPLY_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <vector>

namespace lawa {

template <typename Operator, typename Precond>
struct MyApply
    : public GeneralMatrix<MyApply<Operator, Precond> >
{
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
    typedef flens::DenseVector<Array<int> >                  IntegerDenseVector;

    typedef RealDenseVector::ElementType                     ElementType;
    typedef RealDenseVector::IndexType                       IndexType;


    MyApply(const Operator &A, const Precond &P, double eps);

    void
    operator()(Transpose               trans,
               const double            alpha,
               const RealDenseVector   &v,
               const double            beta,
               RealDenseVector         &w)  const;

    int
    computeK(double                    eps,
             const RealDenseVector     &v,
             const IntegerDenseVector  &vSorted) const;

    int
    numRows() const;

    int
    numCols() const;

    const Operator                          &A;
    const Precond                           &P;
    double                                  eps;
};

} // namespace lawa


namespace flens {

template <typename Operator, typename Precond, typename VX, typename VY>
    void
    mv(Transpose transA, double alpha,
       const MyApply<Operator, Precond> &A, const DenseVector<VX> &x,
       double beta,
       DenseVector<VY> &y);

} // namespace flens


#endif // IRIS_MYAPPLY_H