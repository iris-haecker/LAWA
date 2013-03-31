#ifndef IRIS_MYIO_H
#define IRIS_MYIO_H 1

#include <lawa/flensforlawa.h>
#include <iris/myoperator.h>
#include <iris/myprecond3.h>
#include <iostream>

namespace flens {

template <typename T>
    std::ostream &
    operator<<(std::ostream &out,
               const SparseGeMatrix<CRS<T, CRS_General> > &A);

template <typename T>
    std::ostream &
    operator<<(std::ostream &out, const MyOperator<T> &A);

template <typename T>
    std::ostream &
    operator<<(std::ostream &out, const MyPrecond3<T> &A);

} // namespace flens

#endif // IRIS_MYIO_H
