#ifndef IRIS_MYAPPLY_TCC
#define IRIS_MYAPPLY_TCC 1

#include <cassert>
#include <iris/myapply.h>

namespace lawa {

using namespace lawa;
using namespace std;


template <typename Operator>
MyApply<Operator>::MyApply(const Operator &OperatorA, double _eps, double _tol)
    : A(OperatorA), eps(_eps), tol(_tol),
      _col(A.numCols(), A.firstCol()),
      _row(A.numRows(), A.firstRow()),
      _colCache(A.numCols()), _rowCache(A.numRows())
{
    std::cerr << "_A: [" << A.firstRow() << "," << A.lastRow()
              << "] x [" << A.firstCol() << "," << A.lastCol()
              << "]" << std::endl;
    std::cerr << "_col: [" << _col.firstIndex() << "," << _col.lastIndex()
              << "]" << std::endl;
    std::cerr << "_row: [" << _row.firstIndex() << "," << _row.lastIndex()
              << "]" << std::endl;
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
    if (trans==NoTrans) {
        assert(A.numCols()==v.length());
        assert(A.numRows()==w.length());
    }
    if (trans==Trans) {
        assert(A.numRows()==v.length());
        assert(A.numCols()==w.length());
    }
#   endif
    
    
    std::cerr << "BEGIN: Apply operator" << std::endl;
    std::cerr << "trans = " << ((trans==Trans) ? 'T' : 'N')
              << ", alpha = " << alpha
              << ", beta =" << beta
              << std::endl;
    
    using std::min;

    N     = cachePartialSort(v, tol);
    cerr << "N = " << N << endl;

    int k = computeK(eps, v, N);
    cerr << "k = " << k << endl;


    int maxBin = int(log(double(N))/log(double(2))) + 1;
    maxBin = min(maxBin, k);

    std::cerr << "maxBin = " << maxBin << std::endl;

    for (int r=w.firstIndex(); r<=w.lastIndex(); ++r) {
        w(r) *= beta;
    }

    for (int bin=0; bin<=maxBin; ++bin) {
        std::cerr << "bin = " << bin << ": Caching ... ";
        int j0 = 1<<bin;
        int j1 = min(1<<(bin+1), N);

        for (int j=j0; j<j1; ++j) {
            const double _v = v(vSorted(j));

            if (trans==NoTrans) {
                cacheAkCol(k-bin, vSorted(j));
                std::cerr << "col " << vSorted(j) << " ... ";

                int r0 = _colCache[vSorted(j)-1].firstIndex();
                int r1 = _colCache[vSorted(j)-1].lastIndex();
                
                for (int r=r0; r<r1; ++r) {
                    const double _a = _colCache[vSorted(j)-1](r);
                    w(r) += alpha * _a * _v;
                }
            } else if (trans==Trans) {
                cacheAkRow(k-bin, vSorted(j));
                std::cerr << "row " << vSorted(j) << " ... ";

                int c0 = _rowCache[vSorted(j)-1].firstIndex();
                int c1 = _rowCache[vSorted(j)-1].lastIndex();
                
                for (int c=c0; c<c1; ++c) {
                    const double _a = _rowCache[vSorted(j)-1](c);
                    w(c) += alpha * _a * _v;
                }
            } else {
                assert(0);
            }
        }
        std::cerr << std::endl;
    }
    std::cerr << "END: Apply operator" << std::endl;
}

template <typename Operator>
bool
MyApply<Operator>::cacheAkCol(int k, int col) const
{
    if (_col(col)>=k+1) {
        std::cerr << " [CACHE HIT] ";
        return false;
    }

    int j  = A.getLevelOfCol(col);
    int j1 = std::min(A.j1, j+k);

    int r1 = A.firstRowWithLevel(j-k);
    int r2 = A.lastRowWithLevel(j1);

    _colCache[col-1].engine().resize(r2-r1+1, r1);
    for (int i=r1; i<=r2; ++i) {
        _colCache[col-1](i) = A(i,col);
    }
    _col(col) = k+1;
    return true;
}

template <typename Operator>
bool
MyApply<Operator>::cacheAkRow(int k, int row) const
{
    if (_row(row)>=k+1) {
        std::cerr << " [CACHE HIT] ";
        return false;
    }

    int j  = A.getLevelOfRow(row);
    int j1 = std::min(A.j1, j+k);

    int c1 = A.firstColWithLevel(j-k);
    int c2 = A.lastColWithLevel(j1);

    _rowCache[row-1].engine().resize(c2-c1+1, c1);
    for (int j=c1; j<=c2; ++j) {
        _rowCache[row-1](j) = A(row, j);
    }
    _row(row) = k+1;
    return true;
}

template <typename Operator>
int
MyApply<Operator>::cachePartialSort(const RealDenseVector  &v, double eps) const
{
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

    std::cerr << "BEGIN: computeK" << std::endl;

//
//  Compute k(eps) according to (7.28)
//
    double norm = 0;
    for (int n=1; n<=N; ++n) {
        norm += pow(abs(v(vSorted(n))), 2);
    }
    norm = sqrt(norm);
    std::cerr << "    norm = " << norm << std::endl;

    double semiNorm = 0;
    for (int n=1; n<=N; ++n) {
        double tmp = pow(double(n), double(1)/tau) * abs(v(vSorted(n)));
        if (tmp>semiNorm) {
            semiNorm = tmp;
        }
    }
    std::cerr << "    semiNorm = " << semiNorm << std::endl;
    
    norm += semiNorm;
    
    //int k_eps = 10 * double(1)/s * log(norm/eps) / log(double(2));
    int k_eps = double(1)/s * log(norm/eps) / log(double(2));

//
//  Compute R_k according to (7.29)
//
    int maxBin = (N>=1) ? int(log(double(N))/log(double(2))) + 1
                        : 1;

    std::cerr << "N = " << N << std::endl;

    std::cerr << "maxBin = " << maxBin << std::endl;

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
            std::cerr << "   computeK ==> k = " << k << ", k_eps = " << k_eps << std::endl;
            std::cerr << "END: computeK" << std::endl;
            return std::min(k,maxlevel);
        }
    }

    std::cerr << "    k_eps = " << k_eps << std::endl;
    std::cerr << "    min(k_eps, maxlevel) = " << min(k_eps, maxlevel)
              << std::endl;
    std::cerr << "END: computeK" << std::endl;

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


template <typename MyApply>
MyApply_AtA<MyApply>::MyApply_AtA(const MyApply &_ApplyA)
    : ApplyA(_ApplyA)
{
}

template <typename MyApply>
int
MyApply_AtA<MyApply>::numRows() const
{
    return dim();
}

template <typename MyApply>
int
MyApply_AtA<MyApply>::numCols() const
{
    return dim();
}

template <typename MyApply>
int
MyApply_AtA<MyApply>::dim() const
{
    return ApplyA.numCols();
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
    int N;
    A(transA, alpha, x, N, beta, y);
}


template <typename MyApply, typename VX, typename VY>
void
mv(double alpha,
   const MyApply_AtA<MyApply> &AtA, const DenseVector<VX> &x,
   double beta,
   DenseVector<VY> &y)
{
    DenseVector<VX>  v;
    static int              N;

    AtA.ApplyA(NoTrans, double(1), x, N, double(0), v);
    AtA.ApplyA(Trans, alpha, v, N, beta, y);
}


} // namespace flens

#endif // IRIS_MYAPPLY_TCC