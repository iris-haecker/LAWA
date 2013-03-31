#ifndef IRIS_MYPRECONDID_H
#define IRIS_MYPRECONDID_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename T>
struct MyPrecondId
{
    MyPrecondId();

    T
    operator()(int absoluteIndex) const;
};

} // namespace lawa

#endif // IRIS_MYPRECONDID_H
