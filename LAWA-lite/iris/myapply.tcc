#ifndef IRIS_MYAPPLY_TCC
#define IRIS_MYAPPLY_TCC 1

#include <cassert>
#include <iris/myapply.h>
#include <iris/mysupport.h>
#include <iris/myabssort.h>

namespace lawa {

template <typename Operator, typename Precond>
MyApply<Operator, Precond>::MyApply(const Operator  &OperatorA,
                                    const Precond   &_P, 
                                    double          _eps)
    : A(OperatorA), P(_P), eps(_eps)
{
}

template <typename Operator, typename Precond>
void
MyApply<Operator, Precond>::operator()(Transpose              trans,
                                       const double           alpha,
                                       const RealDenseVector  &v,
                                       const double           beta,
                                       RealDenseVector        &w)  const
{
    using std::abs;
    using std::max;
    using std::min;

    if (beta==double(0)) {
        if ((trans==NoTrans) && (w.length()!=A.numRows())) {
            w.engine().resize(A.numRows(), 1);
        }
        if ((trans==Trans) && (w.length()!=A.numCols())) {
            w.engine().resize(A.numCols(), 1);
        }
    }

#   ifndef NDEBUG
    assert(v.firstIndex()==1);
    assert(w.firstIndex()==1);
    if (trans==NoTrans) {
        assert(A.numCols()==v.length());
        assert(A.numRows()==w.length());
    }
    if (trans==Trans) {
        assert(A.numRows()==v.length());
        assert(A.numCols()==w.length());
    }
#   endif

    if (beta!=double(1)) {
        for (int r=w.firstIndex(); r<=w.lastIndex(); ++r) {
            w(r) *= beta;
        }
    }

//
//  Sort v
//
    std::cerr << "APPLY: Start" << std::endl;

    static IntegerDenseVector  vSorted;
    static IndexSet<int>       vSupport;

    mySupport(v, vSupport);
    vSorted.engine().resize(int(vSupport.size()));
    myAbsSort(v, vSorted);

    int N = myNNZ(v, vSorted);

    //cerr << "v = " << v << endl;
    //cerr << "vSorted = " << vSorted << endl;

    int k      = computeK(eps, v, vSorted);
    int maxBin = int(log(double(N))/log(double(2))) + 1;

    maxBin = max(min(maxBin, k), 0);

    std::cerr << "APPLY: Sorted (N = " << N << ")" << std::endl;


    if (trans==NoTrans) {

        for (int bin=0; bin<=maxBin; ++bin) {

            int p0 = 1<<bin;
            int p1 = min(1<<(bin+1), N+1);

            for (int p=p0; p<p1; ++p) {

                int L  = A.getLevelOfCol(vSorted(p));

                int j0 = max(A.j0, L-(k-bin));
                int j1 = min(A.j1, L+(k-bin));

                for (int j=j0; j<j1; ++j) {

                    int r0 = A.inCol_firstNonZeroWithLevel(vSorted(p), j);
                    int r1 = A.inCol_lastNonZeroWithLevel(vSorted(p), j);
                    r1 = std::min(r1, A.numRows());

                    for (int r=r0; r<=r1; ++r) {
                        w(r) += alpha * P(r) * A(r, vSorted(p)) * v(vSorted(p));
                    }

                }
            }
        }

    } else if (trans==Trans) {

        for (int bin=0; bin<=maxBin; ++bin) {

            int p0 = 1<<bin;
            int p1 = min(1<<(bin+1), N+1);

            std::cerr << "APPLY: bin = " << bin
                      << ", k - bin = " << k - bin
                      << " (maxBin " << maxBin
                      << ", bin size " << p1-p0
                      << ", p0 = " << p0
                      << ", p1 = " << p1
                      << ")" << std::endl;


            for (int p=p0; p<p1; ++p) {

                if (v(vSorted(p))==0) {
                    continue;
                }

                int L  = A.getLevelOfRow(vSorted(p));

                int j0 = max(A.j0, L-(k-bin));
                int j1 = min(A.j1, L+(k-bin));

                for (int j=j0; j<j1; ++j) {

                    int c0 = A.inRow_firstNonZeroWithLevel(vSorted(p), j);
                    int c1 = A.inRow_lastNonZeroWithLevel(vSorted(p), j);
                    c1 = std::min(c1, A.numCols());

                    for (int c=c0; c<=c1; ++c) {
                        w(c) += alpha * P(c) * A(vSorted(p), c) * v(vSorted(p));
                    }

                }
            }
        }

    }
}

template <typename Operator, typename Precond>
int
MyApply<Operator, Precond>::computeK(double                    eps,
                                     const RealDenseVector     &v,
                                     const IntegerDenseVector  &vSorted) const
{
    using std::abs;
    using std::log;
    using std::min;
    using std::pow;
    using std::sqrt;

    const int    N   = vSorted.length();
    const double s   = A.smoothness();
    const double tau = double(1)/(s+double(0.5));

//
//  Compute k(eps) according to (7.28)
//
    double norm = 0;
    for (int n=1; n<=N; ++n) {
        norm += pow(abs(v(vSorted(n))), 2);
    }
    norm = sqrt(norm);

    double semiNorm = 0;
    for (int n=1; n<=N; ++n) {
        double tmp = pow(double(n), double(1)/tau) * abs(v(vSorted(n)));
        if (tmp>semiNorm) {
            semiNorm = tmp;
        }
    }
    
    norm += semiNorm;
    
    int k_eps = double(1)/s * log(norm/eps) / log(double(2));

//
//  Compute R_k according to (7.29)
//
    int maxBin = (N>=1) ? int(log(double(N))/log(double(2))) + 1
                        : 1;

    //std::cerr << "    N = " << N << std::endl;
    //std::cerr << "    maxBin = " << maxBin << std::endl;

    RealDenseVector normSec(_(0,maxBin));

    int bin = 0;
    for (int i=1; i<=N; ++i) {
        normSec(bin) += pow(v(vSorted(i)), 2);
        if (i==(1<<bin)) {
            ++bin;
        }
    }

    //std::cerr << "normSec = " << normSec << std::endl;

    int maxlevel=25;

    for (int k=1; k<=k_eps; ++k) {
        double R_k = 0;

        for (int i=k+1; i<=maxBin; ++i) {
            R_k += normSec(i);
        }
        R_k *= A.CA;
        R_k += pow(2.,-k*s) * normSec(0);

        for (int l=0; l<=k-1; ++l) {
            if (k-l<maxBin) {
                R_k += pow(double(2),-l*s) * normSec(k-l);
            }
        }

        if (R_k<=eps) {
            //std::cerr << "   computeK ==> k = " << k << ", k_eps = " << k_eps << std::endl;
            return std::min(k,maxlevel);
        }
    }

    // std::cerr << "    k_eps = " << k_eps << std::endl;
    // std::cerr << "    min(k_eps, maxlevel) = " << min(k_eps, maxlevel)
    //           << std::endl;
    // std::cerr << "END: computeK" << std::endl;
    return std::min(k_eps, maxlevel);
}

template <typename Operator, typename Precond>
int
MyApply<Operator, Precond>::numRows() const
{
    return A.numRows();
}

template <typename Operator, typename Precond>
int
MyApply<Operator, Precond>::numCols() const
{
    return A.numCols();
}

} // namespace lawa


namespace flens {

template <typename Operator, typename Precond, typename VX, typename VY>
void
mv(Transpose transA, double alpha,
   const MyApply<Operator, Precond> &A, const DenseVector<VX> &x,
   double beta,
   DenseVector<VY> &y)
{
    if (y.length()==0) {
        if (transA==NoTrans) {
            y.engine().resize(A.numRows(), 1);
        } else {
            y.engine().resize(A.numCols(), 1);
        }
    }

#   ifndef NDEBUG
    if (transA==NoTrans) {
        assert(A.numCols()==x.length());
        assert(A.numRows()==y.length());
        assert(x.firstIndex()==1);
        assert(y.firstIndex()==1);
    } else {
        assert(A.numCols()==y.length());
        assert(A.numRows()==x.length());
        assert(x.firstIndex()==1);
        assert(y.firstIndex()==1);
    }
#   endif

    int N;
    A(transA, alpha, x, beta, y);
}

} // namespace flens

#endif // IRIS_MYAPPLY_TCC