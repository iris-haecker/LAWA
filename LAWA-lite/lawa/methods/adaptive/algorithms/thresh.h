#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {

template <typename T, typename Index>
    Coefficients<Lexicographical,T,Index >
    THRESH(const Coefficients<Lexicographical,T,Index > &v, T eta);

template <typename T, typename Index>
    Coefficients<Lexicographical,T,Index >
    THRESH(const Coefficients<AbsoluteValue,T,Index > &v, T eta);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/thresh.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H

