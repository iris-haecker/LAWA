#ifndef IRIS_MYNNZ_TCC
#define IRIS_MYNNZ_TCC 1

#include <cassert>

namespace lawa {

template <typename VV, typename VS>
int
myNNZ(const DenseVector<VV> &v, const DenseVector<VS> &vSorted)
{
    assert(v.firstIndex()==1);
    assert(vSorted.firstIndex()==1);

    int nnz = 0;
    for (int k=1; k<=v.length(); ++k) {
        if (v(vSorted(k))==0) {
            break;
        }
        ++nnz;
    }
    return nnz;
}

} // namespace lawa

#endif // IRIS_MYNNZ_TCC