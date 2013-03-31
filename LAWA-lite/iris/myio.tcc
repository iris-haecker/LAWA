#ifndef IRIS_MYIO_TCC
#define IRIS_MYIO_TCC 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <iostream>

namespace flens {

template <typename T>
std::ostream &
operator<<(std::ostream &out,
           const SparseGeMatrix<CRS<T, CRS_General> > &A)
{
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;

    RealGeMatrix B(A.numRows(), A.numCols());
    for (int r=1; r<A.engine().rows.length(); ++r) {
        for (int k=A.engine().rows(r); k<A.engine().rows(r+1); ++k) {
            int c = A.engine().columns(k);
            B(r, c) = A.engine().values(k);
        }
    }
    out << B;

    return out;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out, const MyOperator<T> &opA)
{
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;

    RealGeMatrix A;

    opA.densify(A, opA.j1V, opA.j1U);

    out << A;

    return out;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out, const MyPrecond3<T> &P)
{
    out << P.values;
    return out;
}

} // namespace flens

#endif // IRIS_MYIO_tCC
