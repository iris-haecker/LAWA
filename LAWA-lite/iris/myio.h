#ifndef IRIS_MYIO_H
#define IRIS_MYIO_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <iostream>

namespace flens {

template <typename T>
    std::ostream &
    operator<<(std::ostream &out,
               const SparseGeMatrix<CRS<T, CRS_General> > &A);

} // namespace flens

#endif // IRIS_MYIO_H