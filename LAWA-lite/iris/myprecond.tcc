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
        p(k) = 1/sqrt(myIntegral(k,0,k,0)+myIntegral(k,1,k,1));
    }
}

template <typename Operator>
double
MyPrecond<Operator>::operator()(int absoluteIndex) const
{
    return p(absoluteIndex);
}

} // namespace lawa

#endif // IRIS_MYPRECOND_TCC