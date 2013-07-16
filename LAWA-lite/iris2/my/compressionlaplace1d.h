#ifndef IRIS2_MY_COMPRESSIONLAPLACE1D_H
#define IRIS2_MY_COMPRESSIONLAPLACE1D_H 1

#include <iris2/my/laplace1d.h>

namespace lawa {

template <typename T>
struct CompressionLaplace1D
{
    const Laplace1D<T>  &A;
    short               s_tilde, jmin, jmax;

    CompressionLaplace1D(const Laplace1D<T> &A);

    void
    setParameters(const IndexSet<Index1D> &LambdaRow);

    IndexSet<Index1D>
    SparsityPattern(const Index1D &lambda_col,
                    const IndexSet<Index1D> &LambdaRow,
                    int J=-1);
};

} // namespace lawa

#endif // IRIS2_MY_COMPRESSIONLAPLACE1D_H