#ifndef IRIS_MYABSSORT_H
#define IRIS_MYABSSORT_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

template <typename VX, typename VI>
    void
    myAbsSort(const DenseVector<VX> &x, DenseVector<VI> &xi);

} // namespace lawa

#endif // IRIS_MYABSSORT_H