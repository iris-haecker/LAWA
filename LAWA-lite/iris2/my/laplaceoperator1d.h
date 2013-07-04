#ifndef IRIS2_MY_LAPLACEOPERATOR_H
#define IRIS2_MY_LAPLACEOPERATOR_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename T>
struct LaplaceOperator1D
{

    LaplaceOperator1D(int d, int d_);

    LaplaceOperator1D(int d, int d_, int jMin);

    LaplaceOperator1D(const MyLaplaceOperator1D &A);

    T
    operator()(XType e1, int j1, int k1,
               XType e2, int j2, int k2) const;

    T
    operator()(const Index1D &row_index, const Index1D &col_index) const;

    CompoundBasis<T>  U, V;
    int               j0;
};

} // namespace lawa

#endif // IRIS2_MY_LAPLACEOPERATOR_H
