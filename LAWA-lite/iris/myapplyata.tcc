#ifndef IRIS_MYAPPLYATA_TCC
#define IRIS_MYAPPLYATA_TCC 1

#include <cassert>
#include <iris/myapplyata.h>

namespace lawa {

using namespace lawa;
using namespace std;

template <typename MyApply>
MyApplyAtA<MyApply>::MyApplyAtA(const MyApply &_ApplyA)
    : ApplyA(_ApplyA)
{
}

template <typename MyApply>
int
MyApplyAtA<MyApply>::numRows() const
{
    return dim();
}

template <typename MyApply>
int
MyApplyAtA<MyApply>::numCols() const
{
    return dim();
}

template <typename MyApply>
int
MyApplyAtA<MyApply>::dim() const
{
    return ApplyA.numCols();
}

} // namespace lawa


namespace flens {

template <typename MyApply, typename VX, typename VY>
void
mv(double alpha,
   const MyApplyAtA<MyApply> &AtA, const DenseVector<VX> &x,
   double beta,
   DenseVector<VY> &y)
{
    DenseVector<VX>  v;
    static int              N;

    AtA.ApplyA(NoTrans, double(1), x, N, double(0), v);
    AtA.ApplyA(Trans, alpha, v, N, beta, y);
}

} // namespace flens

#endif // IRIS_MYAPPLYATA_TCC