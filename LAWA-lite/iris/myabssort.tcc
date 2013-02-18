#ifndef IRIS_MYABSSORT_TCC
#define IRIS_MYABSSORT_TCC 1

#include <algorithm>

namespace lawa {

template <typename V>
struct _MyAbsCompare
{
    typedef typename V::IndexType  IndexType;


    _MyAbsCompare(const V &_x)
        : x(_x)
    {
    }

    bool
    operator() (IndexType i, IndexType j)
    {
        using std::abs;

        return (abs(x(i))>abs(x(j)));
    }

    const V &x;
};


template <typename VX, typename VI>
void
myAbsSort(const DenseVector<VX> &x, DenseVector<VI> &xi)
{
    typedef typename DenseVector<VX>::IndexType  IndexType;

    xi.engine().resize(x.length(), x.firstIndex());

    for (IndexType k=x.firstIndex(); k<=x.lastIndex(); ++k) {
        xi(k) = k;
    }

    _MyAbsCompare<DenseVector<VX> >  cmp(x);

    sort(xi.engine().data(), xi.engine().data()+xi.length(), cmp);
}

} // namespace lawa

#endif // IRIS_MYABSSORT_TCC