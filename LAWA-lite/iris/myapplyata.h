#ifndef IRIS_MYAPPLYATA_H
#define IRIS_MYAPPLYATA_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <iris/myapply.h>
#include <vector>

namespace lawa {

using namespace lawa;
using namespace std;

template <typename MyApply>
struct MyApplyAtA
    : public SymmetricMatrix<MyApplyAtA<MyApply> >
{

    typedef typename MyApply::ElementType   ElementType;
    typedef typename MyApply::IndexType     IndexType;


    MyApplyAtA(const MyApply &ApplyA);

    int
    numRows() const;

    int
    numCols() const;

    int
    dim() const;

    const MyApply   &ApplyA;
};

} // namespace lawa

namespace flens {

template <typename MyApply, typename VX, typename VY>
   void
   mv(double alpha,
      const MyApplyAtA<MyApply> &AtA, const DenseVector<VX> &x,
      double beta,
      DenseVector<VY> &y);

} // namespace flens

#endif // IRIS_MYAPPLYATA_H