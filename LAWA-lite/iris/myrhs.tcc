#ifndef IRIS_MYRHS_TCC
#define IRIS_MYRHS_TCC 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <vector>

namespace lawa {

template <typename T, typename Precond>
MyRhs<T, Precond>::MyRhs(const Function<T> &f,
                         const MyOperator<T> &myOperator,
                         const Precond &P)
{
    using std::pow;
    using std::sqrt;

    MyRhsIntegral<T>   rhsIntegral(f, myOperator.V);

    rhsData.engine().resize(myOperator.numRows(), myOperator.firstRow());

    norm = 0;
    for (int k=rhsData.firstIndex(); k<=rhsData.lastIndex(); ++k) {
        rhsData(k) = P(k)*rhsIntegral(k);
        norm += pow(rhsData(k), 2);
    }
    norm = sqrt(norm);
}

template <typename T, typename Precond>
template <typename VX>
int
MyRhs<T, Precond>::filter(const T &tol, DenseVector<VX> &x) const
{
    using std::abs;

    x.engine().resize(rhsData.length()) || x.engine().fill(0);

    int count = 0;
    for (int k=rhsData.firstIndex(); k<=rhsData.lastIndex(); ++k) {
        if (abs(rhsData(k))>tol) {
            x(k) = rhsData(k);
            ++count;
        }
    }
    return count;
}

template <typename T, typename Precond>
template <typename TI, typename VX>
int
MyRhs<T, Precond>::filter(const T             &tol,
                          const IndexSet<TI>  &Lambda,
                          DenseVector<VX>     &x) const
{
    typedef IndexSet<int>::const_iterator  iterator;

    const int N = Lambda.size();
    x.engine().resize(N);

    int K=1;
    for (iterator k=Lambda.begin(); k!=Lambda.end(); ++k, ++K) {
        if (abs(rhsData(*k))>tol) {
            x(K) = rhsData(*k);
        }
    }
}


} // namespace lawa


#endif // IRIS_MYRHS_TCC