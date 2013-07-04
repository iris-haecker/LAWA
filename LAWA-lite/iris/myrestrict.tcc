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
            B(r, c) = P(*row) * A(*row, *col);
        }
    }
    B.finalize();
}

template <typename CRS, typename MB>
void
myAAt(const SparseGeMatrix<CRS> &A, GeMatrix<MB> &B)
{
    typedef flens::DenseVector<Array<double> >  RealDenseVector;

    const int M = A.numRows();
    const int N = A.numCols();

    RealDenseVector  x(N), y(M);

    B.engine().resize(M, M);

    for (int r=1; r<=M; ++r) {
        x = 0;

        for (int k=A.engine().rows(r); k<A.engine().rows(r+1); ++k) {
            int c = A.engine().columns(k);
            x(c) = A.engine().values(k);
        }
        y = A*x;
        B(_,r) = y;
    }

}

template <typename Operator, typename TI, typename MB>
void
myRestrict(const Operator       &A,
           const IndexSet<TI>   &Lambda,
           int                  k,
           GeMatrix<MB>         &B)
{
    using std::max;
    using std::min;

    typedef IndexSet<int>::const_iterator               iterator;
    typedef SparseGeMatrix<CRS<double, CRS_General> >   SparseMatrix;

    const int N = Lambda.size();
    const int M = A.numRows();

//
//  In A_NM we store the transposed of A restricted to columns of Lambda.
//
    SparseMatrix A_NM(N, M, 2*k*(A.U.d()+A.V.d()));

    int c=1;
    for (iterator col=Lambda.begin(); col!=Lambda.end(); ++col, ++c) {

        int L  = A.getLevelOfCol(*col);
        int j0 = max(A.j0-1, L-k);
        int j1 = min(A.j1V, L+k);

        for (int j=j0; j<=j1; ++j) {

            int r0 = A.inCol_firstNonZeroWithLevel(*col, j);
            int r1 = A.inCol_lastNonZeroWithLevel(*col, j);
            r1 = std::min(r1, M);

            for (int r=r0; r<=r1; ++r) {
                /*
                std::cerr << "A_MN(" << r << ", " << c
                          << ") = A(" << r << ", " << *col << ")  "
                          << A(r, *col)
                          << std::endl;
                */
                A_NM(c, r) = A(r, *col);
            }
        }
    }
    A_NM.finalize();

    //std::cerr << "A_NM = " << A_NM << std::endl;

    //B = transpose(A_MN)*A_MN = A_NM*transpose(A_NM);
    myAAt(A_NM, B);
}

template <typename Operator, typename Precond, typename TI, typename MB>
void
myRestrict(const Operator       &A,
           Precond              &P,
           const IndexSet<TI>   &Lambda,
           int                  &k,
           GeMatrix<MB>         &B,
           double               eps)
{
    using std::abs;
    using std::max;

    typename GeMatrix<MB>::NoView  B_;

    double maxB, maxDiff, diff;

    do {
        myRestrict(A, Lambda, k, B_);
        myRestrict(A, Lambda, ++k, B);

        maxB    = 0;
        maxDiff = abs(B_(1,1) - B(1,1));
        for (int j=1; j<=B.numCols(); ++j) {
            for (int i=1; i<=B.numRows(); ++i) {
                diff = abs(B_(i,j) - B(i,j));
                if (diff>maxDiff) {
                    maxDiff = diff;
                }

                if (abs(B(i,j))>maxB) {
                    maxB = abs(B(i,j));
                }
            }
        }
        std::cerr << "->  maxDiff = " << maxDiff
                  << ",   maxB = " << maxB
                  << ",   maxDiff/maxB = " << (maxDiff/maxB)
                  << std::endl;

    } while (maxDiff>eps);

    --k;

    const int N = B.numRows();

    P.resize(N);
    for (int p=1; p<=N; ++p) {
        P(p) = B(p,p);
    }

    std::cerr << "k = " << k
              << ", maxDiff = " << maxDiff
              << std::endl;
}

template <typename VX, typename TI, typename VY>
void
myRestrict(const DenseVector<VX>  &x,
           const IndexSet<TI>     &Lambda,
           DenseVector<VY>        &y)
{
    typedef IndexSet<int>::const_iterator  iterator;

    const int N = Lambda.size();
    if (y.length()==0) {
        y.engine().resize(N);
    }
    assert(y.length()==N);

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
    if (y.length()==0) {
        y.engine().resize(N);
    }
    assert(y.length()==N);

    int K=1;
    for (iterator k=Lambda.begin(); k!=Lambda.end(); ++k, ++K) {
        y(K) -= x(*k);
    }
}

} // namespace lawa

#endif // IRIS_MYRESTRICT_TCC
