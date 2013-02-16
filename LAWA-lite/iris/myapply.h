#ifndef IRIS_MYAPPLY_H
#define IRIS_MYAPPLY_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <vector>

namespace lawa {

using namespace lawa;
using namespace std;


template <typename Operator>
struct MyApply
    : public GeneralMatrix<MyApply<Operator> >
{
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
    typedef flens::DenseVector<Array<int> >                  IntegerDenseVector;

    typedef RealDenseVector::ElementType                     ElementType;
    typedef RealDenseVector::IndexType                       IndexType;


    MyApply(const Operator &A, double eps, double tol);

    void
    operator()(Transpose              trans,
               const double           alpha,
               const RealDenseVector  &v,
               int                    &N,
               const double           beta,
               RealDenseVector        &w)  const;

    void
    densify(int jMax, RealGeMatrix &MA) const;

    void
    precond(int jMax, RealDenseVector &MP) const;

    int
    cachePartialSort(const RealDenseVector  &v, double eps) const;

    int
    computeK(double eps, const RealDenseVector  &v, int N) const;

    int
    numRows() const;

    int
    numCols() const;

    const Operator                          &A;
    double                                  eps, tol;
    mutable IntegerDenseVector              vSorted;
};

} // namespace lawa


namespace flens {

template <typename Operator, typename VX, typename VY>
    void
    mv(Transpose transA, double alpha,
       const MyApply<Operator> &A, const DenseVector<VX> &x,
       double beta,
       DenseVector<VY> &y);

} // namespace flens


#endif // IRIS_MYAPPLY_H