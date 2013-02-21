#ifndef IRIS_MYPRECOND_H
#define IRIS_MYPRECOND_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename Operator>
struct MyPrecond
{
    typedef flens::DenseVector<Array<double> >      RealDenseVector;
    
    MyPrecond(const Operator &A);

    double
    operator()(int absoluteIndex) const;

    const Operator  &A;
    RealDenseVector p;
};

} // namespace lawa

#endif // IRIS_MYPRECOND_H