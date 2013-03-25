#ifndef IRIS_MYNNZ_H
#define IRIS_MYNNZ_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename VV, typename VS>
    int
    myNNZ(const DenseVector<VV> &v, const DenseVector<VS> &vSorted);

} // namespace lawa

#endif // IRIS_MYNNZ_H