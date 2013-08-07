#ifndef IRIS2_MY_COMPRESSIONLAPLACE1D_TCC
#define IRIS2_MY_COMPRESSIONLAPLACE1D_TCC 1

#include <iris2/iris2.h>

namespace lawa {

template <typename T>
CompressionLaplace1D<T>::CompressionLaplace1D(const Laplace1D<T> &_A)
    : A(_A), s_tilde(-1), jmin(100), jmax(-30)
{
}

template <typename T>
void
CompressionLaplace1D<T>::setParameters(const IndexSet<Index1D> &LambdaRow)
{
    using std::max;
    using std::min;

    typedef typename IndexSet<Index1D>::const_iterator SetIt;

    s_tilde = -1;
    jmin = 100;
    jmax = -30;
    for (SetIt lambda=LambdaRow.begin(); lambda!=LambdaRow.end(); ++lambda) {
        jmin = min(jmin,(*lambda).j);
        jmax = max(jmax,(*lambda).j);
    }
    s_tilde = jmax-jmin;
}

template <typename T>
IndexSet<Index1D>
CompressionLaplace1D<T>::SparsityPattern(const Index1D &lambda_col,
                                         const IndexSet<Index1D> &LambdaRow,
                                         int J)
{
    using std::max;
    using std::min;

    typedef typename IndexSet<Index1D>::const_iterator SetIt;

    IndexSet<Index1D> LambdaRowSparse;

    int s = max(abs(lambda_col.j-jmin),abs(lambda_col.j-jmax));
    s = max(s,int(s_tilde));

    if (J!=-1) {
        s = min(s,J);
    }

    //Compression level J>s_tilde as indices corresponding to level differences
    //larger than s_tilde cannot appear in LambdaRow

    IndexSet<Index1D> Lambda_x = lambdaTilde1d(lambda_col, A, s, jmin, jmax);

    for (SetIt lambda=Lambda_x.begin(); lambda!=Lambda_x.end(); ++lambda) {
        if (LambdaRow.count(*lambda)>0) {
            LambdaRowSparse.insert(*lambda);
        }
    }
    return LambdaRowSparse;
}

} // namespace lawa

#endif // IRIS2_MY_COMPRESSIONLAPLACE1D_TCC