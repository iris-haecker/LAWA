#ifndef IRIS_MYOPERATOR_H
#define IRIS_MYOPERATOR_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename T>
struct MyOperator
    : public GeneralMatrix<MyOperator<T> >
{
    typedef flens::GeMatrix<FullStorage<T, ColMajor> >  RealGeMatrix;
    typedef T                                           ElementType;
    typedef int                                         IndexType;

    MyOperator(int d, int d_, int _jMax);

    T
    operator()(int row, int col) const;

    int
    firstRow() const;

    int
    lastRow() const;

    int
    numRows() const;

    int
    firstCol() const;

    int
    lastCol() const;

    int
    numCols() const;

    int
    getLevelOfRow(int row) const;

    int
    getLevelOfCol(int col) const;

    int
    firstRowWithLevel(int j) const;

    int
    lastRowWithLevel(int j) const;

    int
    firstColWithLevel(int j) const;

    int
    lastColWithLevel(int j) const;

    double
    smoothness() const;

    int
    inCol_firstNonZeroWithLevel(int col, int j) const;

    int
    inCol_lastNonZeroWithLevel(int col, int j) const;

    void
    densify(RealGeMatrix &MA, int jMax = -1, bool brute = false) const;

    template <typename CRS>
        void
        restriction(const IndexSet<int> &indices, SparseGeMatrix<CRS> &A) const;

    MyBasis<double>  U, V;
    int              j0, j1;
    Range<int>       KU, KV;
    double           CA;
};

} // namespace lawa


namespace flens {

template <typename T, typename VX, typename VY>
    void
    mv(Transpose transA, double alpha,
       const MyOperator<T> &A, const DenseVector<VX> &x,
       double beta,
       DenseVector<VY> &y);

} // namespace flens

#endif // IRIS_MYOPERATOR_H