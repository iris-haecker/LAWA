#ifndef IRIS2_MY_LAPLACE1D_H
#define IRIS2_MY_LAPLACE1D_H 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T>
struct Laplace1D
{

    Laplace1D(int d, int d_);

    Laplace1D(int d, int d_, int jMin);

    Laplace1D(const Laplace1D &A);

    T
    operator()(XType e1, int j1, int k1,
               XType e2, int j2, int k2) const;

    T
    pH1(XType e1, int j1, int k1,
        XType e2, int j2, int k2) const;

    int
    minK2(int j1, int k1, XType e1, int j2, XType e2) const;

    int
    maxK2(int j1, int k1, XType e1, int j2, XType e2) const;

    const CompoundBasis<double> &
    basisRow() const;

    const CompoundBasis<double> &
    basisCol() const;

    T
    operator()(const Index1D &row_index, const Index1D &col_index) const;

    T
    pH1(const Index1D &row_index, const Index1D &col_index) const;

    int
    d() const;

    CompoundBasis<double>  U, V;
    int                    j0;
};

} // namespace lawa

#endif // IRIS2_MY_LAPLACE1D_H
