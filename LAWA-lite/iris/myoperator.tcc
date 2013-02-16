#ifndef IRIS_MYOPERATOR_TCC
#define IRIS_MYOPERATOR_TCC 1

#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
MyOperator<T>::MyOperator(int d, int d_, int jMax)
    : U(d, d_), V(d, d_),
      j0(std::min(U.j0(), V.j0())), j1(jMax),
      KU(U.getFirstAbsoluteIndex(j0-1), U.getLastAbsoluteIndex(j1)),
      KV(V.getFirstAbsoluteIndex(j0-1), V.getLastAbsoluteIndex(j1)),
      CA(1)
{
    assert(jMax>=j0);

    V.basisRight.enforceBoundaryCondition<lawa::DirichletBC>();
    V.basisLeft.enforceBoundaryCondition<lawa::DirichletBC>();
    U.basisRight.enforceBoundaryCondition<lawa::DirichletBC>();
    U.basisLeft.enforceBoundaryCondition<lawa::DirichletBC>();
}

template <typename T>
T
MyOperator<T>::operator()(int row, int col) const
{
    static MyIntegral<double> myIntegral(V, U);

    //return -myIntegral(row, 0, col, 1) + myIntegral(row, 0, col, 0);
    //return -myIntegral(row, 1, col, 1) + myIntegral(row, 0, col, 0);
    return -myIntegral(row, 1, col, 1);
    //return myIntegral(row, 0, col, 0);
}

template <typename T>
int
MyOperator<T>::firstRow() const
{
    return V.getFirstAbsoluteIndex(j0-1);
}

template <typename T>
int
MyOperator<T>::lastRow() const
{
    return V.getLastAbsoluteIndex(j1);
}

template <typename T>
int
MyOperator<T>::numRows() const
{
    return lastRow() - firstRow() + 1;
}

template <typename T>
int
MyOperator<T>::firstCol() const
{
    return U.getFirstAbsoluteIndex(j0-1);
}

template <typename T>
int
MyOperator<T>::lastCol() const
{
    return U.getLastAbsoluteIndex(j1);
}

template <typename T>
int
MyOperator<T>::numCols() const
{
    return lastCol() - firstCol() + 1;
}

template <typename T>
int
MyOperator<T>::getLevelOfRow(int row) const
{
    return V.getLevel(row);
}

template <typename T>
int
MyOperator<T>::getLevelOfCol(int col) const
{
    return U.getLevel(col);
}

template <typename T>
int
MyOperator<T>::firstRowWithLevel(int j) const
{
    return V.getFirstAbsoluteIndex(j);
}

template <typename T>
int
MyOperator<T>::lastRowWithLevel(int j) const
{
    return V.getLastAbsoluteIndex(j);
}

template <typename T>
int
MyOperator<T>::firstColWithLevel(int j) const
{
    return U.getFirstAbsoluteIndex(j);
}

template <typename T>
int
MyOperator<T>::lastColWithLevel(int j) const
{
    return U.getLastAbsoluteIndex(j);
}

template <typename T>
double
MyOperator<T>::smoothness() const
{
    return U.d()-1.5;
}

template <typename T>
int
MyOperator<T>::inCol_firstNonZeroWithLevel(int col, int j) const
{
    return V.getFirstAbsoluteIndex(j)
         + V.minK(j, U.support(col).l1)
         - V.firstIndex(j);
}

template <typename T>
int
MyOperator<T>::inCol_lastNonZeroWithLevel(int col, int j) const
{
    return V.getFirstAbsoluteIndex(j)
         + V.maxK(j, U.support(col).l2)
         - V.firstIndex(j);
}

template <typename T>
void
MyOperator<T>::densify(RealGeMatrix &MA, int jMax, bool brute) const
{
    if (jMax<j0-1) {
        jMax = j1;
    }

    int r1 = lastRowWithLevel(jMax);
    int c1 = lastColWithLevel(jMax);

    MA.engine().resize(r1, c1);

    for (int c=1; c<=c1; ++c) {

        if (brute) {

            for (int r=1; r<=r1; ++r) {
                MA(r, c) = operator()(r, c);
            }

        } else {

            for (int j=j0-1; j<=jMax; ++j) {
                int _r0 = inCol_firstNonZeroWithLevel(c, j);
                int _r1 = inCol_lastNonZeroWithLevel(c, j);

                for (int r=_r0; r<=_r1; ++r) {
                    MA(r, c) = operator()(r, c);
                }
            }

        }
    }
    
}


} // namespace lawa


namespace flens {

template <typename T, typename VX, typename VY>
void
mv(Transpose transA, double alpha,
   const MyOperator<T> &A, const DenseVector<VX> &x,
   double beta,
   DenseVector<VY> &y)
{
#   ifndef NDEBUG
    assert(x.firstIndex()==1);
    assert(y.firstIndex()==1);
    if (transA==NoTrans) {
        assert(A.numRows()==y.length());
        assert(A.numCols()==x.length());
    } else if (transA==Trans) {
        assert(A.numRows()==x.length());
        assert(A.numCols()==y.length());
    } else {
        assert(0);
    }
#   endif

    std::cerr << "DenseVector = MyOperator * DenseVector" << std::endl;

    if (transA==NoTrans) {

        for (int r=1; r<=y.length(); ++r) {
            y(r) *= beta;
        }

        long skip, nnz, skipTotal, nnzTotal;
        for (int c=1; c<=x.length(); ++c) {
            skip = 0;
            nnz  = 0;

            //std::cerr << c << "  ";
            for (int j=A.j0; j<=A.j1; ++j) {

                int r0 = A.inCol_firstNonZeroWithLevel(c, j);
                int r1 = A.inCol_lastNonZeroWithLevel(c, j);
                
                //std::cerr << r0 << ":" << (r1-r0+1) << " ";

                for (int r=r0; r<=r1; ++r) {
                    const T value = A(r,c);
                    y(r) += alpha*value*x(c);
                    if (std::abs(value)<100*std::numeric_limits<T>::epsilon()) {
                        ++skip;
                        ++skipTotal;
                    } else {
                        ++nnz;
                        ++nnzTotal;
                    }
                }

            }
            /*
            std::cerr << " (skip = " << skip
                      << ",  nnz = " << nnz
                      << ")" << std::endl;
            */
        }
        std::cerr << " (skipTotal = " << skipTotal
                  << ",  nnzTotal = " << nnzTotal
                  << ")" << std::endl;
    } else {
        assert(0);
    }
}

} // namespace flens

#endif // IRIS_MYOPERATOR_TCC