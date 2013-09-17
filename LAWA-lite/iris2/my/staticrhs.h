#ifndef IRIS2_MY_STATICRHS_H
#define IRIS2_MY_STATICRHS_H 1

#include <lawa/methods/adaptive/algorithms/thresh.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename Index, typename Preconditioner>
class StaticRHS
{
    public:
        const Preconditioner &P;
        Coefficients<Lexicographical,T,Index> rhs_data;
        Coefficients<AbsoluteValue,T,Index>   rhs_abs_data;

        StaticRHS(const Preconditioner                        &P,
                  const Coefficients<Lexicographical,T,Index> &_rhs_data);

        T
        operator()(const Index &lambda);

        Coefficients<Lexicographical,T,Index>
        operator()(const IndexSet<Index> &Lambda);

        Coefficients<Lexicographical,T,Index>
        operator()(T tol);
};

}    //namespace lawa

#endif // IRIS2_MY_STATICRHS_H
