#ifndef IRIS_MYSUPPORT_H
#define IRIS_MYSUPPORT_H 1

#include <lawa/flensforlawa.h>
#include <lawa/lawa.h>

namespace lawa {

template <typename VX, typename TI>
    void
    mySupport(const DenseVector<VX> &x, IndexSet<TI> &Lambda);

template <typename IndexType> 
struct lt<Lexicographical, IndexType>
{
    bool
    operator()(const IndexType &left, const IndexType &right) const;
};

} // namespace lawa

#endif // IRIS_MYSUPPORT_H