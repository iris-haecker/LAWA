#ifndef IRIS_MYRESTRICT_TCC
#define IRIS_MYRESTRICT_TCC 1

#include <lawa/flensforlawa.h>
#include <lawa/lawa.h>

namespace lawa {


template <typename Operator, typename Precond, typename TI, typename CRS>
void
myRestrict(const Operator       &A,
           const Precond        &P,
           const IndexSet<TI>   &Lambda,
           SparseGeMatrix<CRS>  &B)
{
    typedef IndexSet<int>::const_iterator  iterator;

    const int N = Lambda.size();
    B.resize(N, N, 2*(A.U.d()+A.V.d()));

    int r=1, c;
    for (iterator row=Lambda.begin(); row!=Lambda.end(); ++row, ++r) {
        c = 1;
        for (iterator col=Lambda.begin(); col!=Lambda.end(); ++col, ++c) {

            std::cerr << "B(" << r << ", " << c
                      << ") = P(" << *row
                      << ") * A(" << *row << ", " << *col << ")"
                      << " = " << P(*row) * A(*row, *col)
                      << std::endl;
            B(r, c) = P(*row) * A(*row, *col);
        }
    }
    B.finalize();
}

template <typename VX, typename TI, typename VY>
void
myRestrict(const DenseVector<VX>  &x,
           const IndexSet<TI>     &Lambda,
           DenseVector<VY>        &y)
{
    typedef IndexSet<int>::const_iterator  iterator;

    const int N = Lambda.size();
    y.engine().resize(N);

    int K=1;
    for (iterator k=Lambda.begin(); k!=Lambda.end(); ++k, ++K) {
        y(K) = x(*k);
    }
}

template <typename VX, typename TI, typename VY>
void
myRestrictSub(const DenseVector<VX>  &x,
              const IndexSet<TI>     &Lambda,
              DenseVector<VY>        &y)
{
    typedef IndexSet<int>::const_iterator  iterator;

    const int N = Lambda.size();
    y.engine().resize(N);

    int K=1;
    for (iterator k=Lambda.begin(); k!=Lambda.end(); ++k, ++K) {
        y(K) -= x(*k);
    }
}

} // namespace lawa

#endif // IRIS_MYRESTRICT_TCC