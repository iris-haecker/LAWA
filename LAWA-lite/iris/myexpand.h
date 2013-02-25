#ifndef IRIS_MYEXPAND_H
#define IRIS_MYEXPAND_H 1

#include <lawa/flensforlawa.h>
#include <lawa/lawa.h>

namespace lawa {

template <typename VX, typename TI, typename VY>
    void
    myExpandAdd(const DenseVector<VX>  &x,
                const IndexSet<TI>     &Lambda,
                DenseVector<VY>        &y);

} // namespace lawa

#endif // IRIS_MYEXPAND_H