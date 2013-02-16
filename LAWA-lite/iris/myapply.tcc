#ifndef IRIS_MYAPPLY_TCC
#define IRIS_MYAPPLY_TCC 1

#include <cassert>
#include <iris/myapply.h>

namespace lawa {

using namespace lawa;
using namespace std;


template <typename Operator>
MyApply<Operator>::MyApply(const Operator &OperatorA, double _eps, double _tol)
    : A(OperatorA), eps(_eps), tol(_tol)
{
}

template <typename Operator>
void
MyApply<Operator>::operator()(Transpose              trans,
                              const double           alpha,
                              const RealDenseVector  &v,
                              int                    &N,
                              const double           beta,
                              RealDenseVector        &w)  const
{
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

    for (int r=w.firstIndex(); r<=w.lastIndex(); ++r) {
        w(r) *= beta;
    }

    using std::min;

    N     = cachePartialSort(v, tol);
    //cerr << "v = " << v << endl;
    //cerr << "vSorted = " << vSorted << endl;
    cerr << "N = " << N << endl;

    int k = computeK(eps, v, N);
    k = min(0, k);
    cerr << "k = " << k << endl;

    int maxBin = int(log(double(N))/log(double(2))) + 1;
    maxBin = min(maxBin, k);

    for (int bin=0; bin<=maxBin; ++bin) {
        int j0 = 1<<bin;
        int j1 = min(1<<(bin+1), N+1);

        /*
        std::cerr << "bin = " << bin
                  << ", j0 = " << j0
                  << ", j1 = " << j1
                  << std::endl;
        */

        for (int p=j0; p<j1; ++p) {

            /*
            std::cerr << "p = " << p
                      << ", vSorted(p) = " << vSorted(p)
                      << std::endl;
            */

            if (trans==NoTrans) {

                int L  = A.getLevelOfCol(vSorted(p));
//                std::cerr << "col = " << vSorted(p) << std::endl;

                for (int j=min(A.j0, L-(k-bin)); j<max(A.j1, L+(k-bin)); ++j) {

                    int r0 = A.inCol_firstNonZeroWithLevel(vSorted(p), j);
                    int r1 = A.inCol_lastNonZeroWithLevel(vSorted(p), j);
                    r1 = std::min(r1, A.numRows());

/*
                    std::cerr << "(" << j << ")  "
                              << r0 << " : " << r1 << std::endl;
*/

                    for (int r=r0; r<=r1; ++r) {
                        w(r) += alpha * A(r, vSorted(p)) * v(vSorted(p));
                    }
                }

            } else if (trans==Trans) {

                int j  = A.getLevelOfRow(vSorted(p));

                int c0 = A.firstColWithLevel(j-(k-bin));
                int c1 = A.lastColWithLevel(j+(k-bin));
                c1 = std::min(c1, A.numCols());

                for (int c=c0; c<=c1; ++c) {
                    w(c) += alpha * v(vSorted(p)) * A(vSorted(p), c);
                }

            } else {
                assert(0);
            }
        }
    }
}

template <typename Operator>
void
MyApply<Operator>::densify(int jMax, RealGeMatrix &MA) const
{
    int r1 = A.lastRowWithLevel(jMax);
    int c1 = A.lastColWithLevel(jMax);

    MA.engine().resize(r1, c1);

    for (int c=1; c<=c1; ++c) {

        if (c%100 == 0) {
            std::cerr << "c = " << c << std::endl;
        }

        for (int j=A.j0-1; j<=jMax; ++j) {
            int r1 = A.inCol_firstNonZeroWithLevel(c, j);
            int r2 = A.inCol_lastNonZeroWithLevel(c, j);
            
            // std::cerr << r1 << " " << r2 << std::endl;
            
            for (int r=r1; r<=r2; ++r) {
                MA(r, c) = A(r, c);
            }
        }
    }
}

template <typename Operator>
void
MyApply<Operator>::precond(int jMax, RealDenseVector &VP) const
{
    using std::min;

    int r1 = A.lastRowWithLevel(jMax);
    int c1 = A.lastColWithLevel(jMax);

    VP.engine().resize(min(r1, c1), 1);

    for (int p=1; p<=min(r1,c1); ++p) {
        int j = A.getLevelOfRow(p);
        VP(p) = 1/double(1<<(j+1));
    }
}



template <typename Operator>
int
MyApply<Operator>::cachePartialSort(const RealDenseVector  &v, double eps) const
{
    std::cerr << "start cachePartialSort ... ";
    using std::abs;

    int N = 0;

    vSorted.engine().resize(v.length());
    vSorted = 0;

    for (int i=1, I=1; i<=v.length(); ++i) {
        if (abs(v(i))>eps) {
            vSorted(I) = i;
            ++I;
            ++N;
        }
    }

    bool swapped;

    do {
        swapped = false;
        for (int i=1; i<N; ++i) {
            if (abs(v(vSorted(i)))<abs(v(vSorted(i+1)))) {
                swap(vSorted(i), vSorted(i+1));
                swapped = true;
            }
        }
    } while (swapped);
    std::cerr << "done." << endl;

    return N;
}

template <typename Operator>
int
MyApply<Operator>::computeK(double eps, const RealDenseVector  &v, int N) const
{
    using std::abs;
    using std::log;
    using std::min;
    using std::pow;
    using std::sqrt;

    double s   = A.smoothness();
    double tau = double(1)/(s+double(0.5));

//
//  Compute k(eps) according to (7.28)
//
    double norm = 0;
    for (int n=1; n<=N; ++n) {
        norm += pow(abs(v(vSorted(n))), 2);
    }
    norm = sqrt(norm);
    //std::cerr << "    norm = " << norm << std::endl;

    double semiNorm = 0;
    for (int n=1; n<=N; ++n) {
        double tmp = pow(double(n), double(1)/tau) * abs(v(vSorted(n)));
        if (tmp>semiNorm) {
            semiNorm = tmp;
        }
    }
    //std::cerr << "    semiNorm = " << semiNorm << std::endl;
    
    norm += semiNorm;
    
    //int k_eps = 10 * double(1)/s * log(norm/eps) / log(double(2));
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

        // std::cerr << "R_k = " << R_k << std::endl;

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

template <typename Operator>
int
MyApply<Operator>::numRows() const
{
    return A.numRows();
}

template <typename Operator>
int
MyApply<Operator>::numCols() const
{
    return A.numCols();
}

} // namespace lawa


namespace flens {

template <typename Operator, typename VX, typename VY>
void
mv(Transpose transA, double alpha,
   const MyApply<Operator> &A, const DenseVector<VX> &x,
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
    A(transA, alpha, x, N, beta, y);
}

} // namespace flens

#endif // IRIS_MYAPPLY_TCC