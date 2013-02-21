#ifndef IRIS_MYSUPPORT_TCC
#define IRIS_MYSUPPORT_TCC 1

#include <lawa/flensforlawa.h>

namespace lawa {

template <typename VX, typename TI>
void
mySupport(const DenseVector<VX> &x, IndexSet<TI> &Lambda)
{
    using std::abs;

    typedef typename DenseVector<VX>::ElementType ElementType;
    typedef typename DenseVector<VX>::IndexType   IndexType;

    assert(x.firstIndex()==IndexType(1));

    const ElementType Zero(0);

    Lambda.clear();
    for (IndexType k=1; k<=x.length(); ++k) {
        if (abs(x(k))>Zero) {
            Lambda.insert(k);
        }
    }
}

template <typename IndexType> 
bool
lt<Lexicographical, IndexType>::operator()(const IndexType &left,
                                           const IndexType &right) const
{
    return (left < right);
}

} // namespace lawa

#endif // IRIS_MYSUPPORT_TCC