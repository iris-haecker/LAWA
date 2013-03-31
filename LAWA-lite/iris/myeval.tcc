#ifndef IRIS_MYEVAL_TCC
#define IRIS_MYEVAL_TCC 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <fstream>

namespace lawa {

template <typename T>
MyEval<T>::MyEval(const MyBasis<T> &_U, const CoeffVector &_u)
    : U(_U), u(_u)
{
    assert(u.firstIndex()==1);
}

template <typename T>
T
MyEval<T>::operator()(T x) const
{
    T y = 0;
    for (int p=1; p<=u.length(); ++p) {
        y += u(p)*U(x, p, 0);
    }
    return y;
}

template <typename T>
void
MyEval<T>::dump(int N, const char *file) const
{
    std::fstream   out(file, std::fstream::out);

    for (int i=0; i<=N; ++i) {
        const double x = double(i)/N;
        out << x << " " << operator()(x) << std::endl;
    }
}


} // namespace lawa

#endif // IRIS_MYEVAL_TCC
