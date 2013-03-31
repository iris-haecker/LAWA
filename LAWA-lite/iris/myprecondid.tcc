#ifndef IRIS_MYPRECONDID_TCC
#define IRIS_MYPRECONDID_TCC 1

#include <iris/mybasis.h>

namespace lawa {

template <typename T>
MyPrecondId<T>::MyPrecondId()
{
}

template <typename T>
T
MyPrecondId<T>::operator()(int absoluteIndex) const
{
    return T(1);
}

} // namespace lawa

#endif // IRIS_MYPRECONDID_TCC
