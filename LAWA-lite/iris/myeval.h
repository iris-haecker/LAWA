#ifndef IRIS_MYEVAL_H
#define IRIS_MYEVAL_H 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <vector>

namespace lawa {

template <typename T>
struct MyEval
{
    typedef DenseVector<Array<T> >  CoeffVector;

    MyEval(const MyBasis<T> &U, const CoeffVector &u);

    T
    operator()(T x) const;

    void
    dump(int N, const char *file) const;

    void
    dump(int N, const Function<T> &refSol, const char *file) const;

    T
    diff_L1(int N, const Function<T> &f) const;

    T
    diff_L2(int N, const Function<T> &f) const;

    T
    diff_LInf(int N, const Function<T> &f) const;

    const CoeffVector  &u;
    const MyBasis<T>   &U;
};


} // namespace lawa

#endif // IRIS_MYEVAL_H