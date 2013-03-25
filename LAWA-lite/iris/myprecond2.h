#ifndef IRIS_MYPRECOND2_H
#define IRIS_MYPRECOND2_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename Operator>
struct MyPrecond2
{
    typedef flens::DenseVector<Array<double> >      RealDenseVector;
    
    MyPrecond2(const Operator &A);

    double
    operator()(int absoluteIndex) const;

    const Operator  &A;
};

} // namespace lawa

#endif // IRIS_MYPRECOND2_H