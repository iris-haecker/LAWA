#ifndef IRIS_MYEXPAND_TCC
#define IRIS_MYEXPAND_TCC 1

#include <lawa/flensforlawa.h>
#include <lawa/lawa.h>

namespace lawa {

template <typename VX, typename TI, typename VY>
void
myExpandAdd(const DenseVector<VX>  &x,
            const IndexSet<TI>     &Lambda,
            DenseVector<VY>        &y)
{
    typedef IndexSet<int>::const_iterator  iterator;

    const int N = Lambda.size();
    assert(x.length()==N);

    int K=1;
    for (iterator k=Lambda.begin(); k!=Lambda.end(); ++k, ++K) {
        y(*k) += x(K);
    }
}

} // namespace lawa

#endif // IRIS_MYEXPAND_TCC