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

    template <typename VU>
        void
        apply(DenseVector<VU> &u) const;

    template <typename IS, typename VU>
        void
        apply(const IndexSet<IS> &Lambda, DenseVector<VU> &u) const;

    const Operator  &A;
    RealDenseVector p;
};

} // namespace lawa

#endif // IRIS_MYPRECOND_H