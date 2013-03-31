#ifndef IRIS_MYPRECOND3_TCC
#define IRIS_MYPRECOND3_TCC 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>

namespace lawa {

template <typename T>
MyPrecond3<T>::MyPrecond3()
{
}

template <typename T>
void
MyPrecond3<T>::resize(int length)
{
    values.engine().resize(length, 1);
}

template <typename T>
T &
MyPrecond3<T>::operator()(int absoluteIndex)
{
    return values(absoluteIndex);
}

} // namespace lawa

#endif // IRIS_MYPRECOND3_TCC
