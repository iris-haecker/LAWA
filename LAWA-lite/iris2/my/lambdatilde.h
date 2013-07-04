#ifndef IRIS2_MY_LAMBDATILDE_H
#define IRIS2_MY_LAMBDATILDE_H 1

#include <iris2/iris2.h>

namespace lawa {

template <typename T>
    IndexSet<Index1D>
    lambdaTilde1d(const Index1D &lambda, const LaplaceOperator1D<T> &A, 
                  int sTilde, int jMin, int jMax);

} // namespace lawa

#endif // IRIS2_MY_LAMBDATILDE_H