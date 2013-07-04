#ifndef IRIS_MYRHS2_TCC
#define IRIS_MYRHS2_TCC 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <vector>

namespace lawa {

template <typename T, typename Precond>
MyRhs2<T, Precond>::MyRhs2(const Function<T>    &f,
                           const MyOperator<T>  &myOperator,
                           const Precond        &P,
                           T                    eps)
{
    using std::abs;
    using std::pow;
    using std::sqrt;

    MyRhsIntegral<T>                rhsIntegral(f, myOperator.V);
    DenseVector<Array<T> >          b;

//
//  Eval right-hand side b
//
    b.engine().resize(myOperator.numRows(), myOperator.firstRow());

    int count = 0;

    norm = 0;
    for (int k=b.firstIndex(); k<=b.lastIndex(); ++k) {
        b(k) = rhsIntegral(k);
        //std::cerr << "b(" << k << ") = " << b(k) << std::endl;
        if (abs(b(k))<eps) {
            b(k) = 0;
            ++count;
        }
        norm += pow(b(k), 2);
    }
    norm = sqrt(norm);
    std::cerr << "computed b. Norm = " << norm
              << ", eps = " << eps
              << ", norm*eps = " << norm*eps
              << std::endl;
    std::cerr << "count = " << count << std::endl;

//
//  Compute an approximation of right-hand side A^T*b
//
#   ifdef  USE_APPLY
    MyApply<MyOperator<T>, Precond> A(myOperator, P, norm*eps);
    rhsData = transpose(A)*b;

#   else
    rhsData = transpose(myOperator)*b;

#   endif


//
//  Apply the preconditioner and compute the norm of P*A^T*b
//
    norm = 0;
    for (int k=rhsData.firstIndex(); k<=rhsData.lastIndex(); ++k) {
        rhsData(k) *= P(k);
        norm += pow(rhsData(k), 2);
    }
    norm = sqrt(norm);
    std::cerr << "computed norm" << std::endl;
}

template <typename T, typename Precond>
template <typename VX>
int
MyRhs2<T, Precond>::filter(const T &tol, DenseVector<VX> &x) const
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
MyRhs2<T, Precond>::filter(const T             &tol,
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


#endif // IRIS_MYRHS2_TCC
