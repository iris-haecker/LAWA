#ifndef IRIS_MYPRECOND3_H
#define IRIS_MYPRECOND3_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>

namespace lawa {

template <typename T>
struct MyPrecond3
{
    typedef flens::DenseVector<Array<T> >      RealDenseVector;

    MyPrecond3();

    void
    resize(int length);

    T &
    operator()(int absoluteIndex);

    RealDenseVector  values;
};

} // namespace lawa

#endif // IRIS_MYPRECOND3_H
