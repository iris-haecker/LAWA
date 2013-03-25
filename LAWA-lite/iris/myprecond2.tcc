#ifndef IRIS_MYPRECOND2_TCC
#define IRIS_MYPRECOND2_TCC 1

#include <iris/mybasis.h>

namespace lawa {

template <typename Operator>
MyPrecond2<Operator>::MyPrecond2(const Operator &_A)
    : A(_A)
{
    using std::sqrt;
}

template <typename Operator>
double
MyPrecond2<Operator>::operator()(int absoluteIndex) const
{
    return 1;
}

} // namespace lawa

#endif // IRIS_MYPRECOND2_TCC