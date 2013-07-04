#ifndef IRIS_MYPRECOND_TCC
#define IRIS_MYPRECOND_TCC 1

#include <iris/mybasis.h>

namespace lawa {

template <typename Operator>
MyPrecond<Operator>::MyPrecond(const Operator &_A)
    : A(_A)
{
    using std::sqrt;

    MyIntegral<double> myIntegral(A.U, A.U);

    p.engine().resize(A.numCols());

    for (int k=1; k<=A.numCols(); ++k) {
        p(k) = double(1)/sqrt(myIntegral(k,0,k,0)+myIntegral(k,1,k,1));
    }
}

template <typename Operator>
double
MyPrecond<Operator>::operator()(int absoluteIndex) const
{
    return p(absoluteIndex);
}

template <typename Operator>
template <typename VU>
void
MyPrecond<Operator>::apply(DenseVector<VU> &u) const
{
    for (int k=u.firstIndex(); k<=u.lastIndex(); ++k) {
        u(k) *= p(k);
    }
}

template <typename Operator>
template <typename IS, typename VU>
void
MyPrecond<Operator>::apply(const IndexSet<IS> &Lambda, DenseVector<VU> &u) const
{
    typedef IndexSet<int>::const_iterator  iterator;

    for (iterator k=Lambda.begin(); k!=Lambda.end(); ++k) {
        u(k) *= p(k);
    }
}

} // namespace lawa

#endif // IRIS_MYPRECOND_TCC