#ifndef IRIS2_MY_COEFFICENTTOOLS_H
#define IRIS2_MY_COEFFICENTTOOLS_H 1

#include <map>

namespace lawa {

template <typename T>
using CoefficientVec = Coefficients<Lexicographical,T,Index1D>;

template <typename T>
using CoefficientRows = std::map<Index1D,
                                 CoefficientVec<T>,
                                 lt<Lexicographical, Index1D> >;

template <typename T>
using CoefficientCols = std::map<Index1D,
                                 CoefficientVec<T>,
                                 lt<Lexicographical, Index1D> >;

template <typename T>
    CoefficientCols<T>
    splitCols(const Coefficients<Lexicographical,T,Index2D> &v);

template <typename T>
    CoefficientRows<T>
    splitRows(const CoefficientCols<T> &colMap);

template <typename T>
    CoefficientRows<T>
    splitRows(const Coefficients<Lexicographical,T,Index2D> &v);

template <typename T>
    Coefficients<Lexicographical,T,Index2D>
    joinRows(const CoefficientRows<T> &rowMap);

template <typename T>
    void
    join(Coefficients<Lexicographical,T,Index2D> &v,
         Index1D row, const CoefficientVec<T> &x);

template <typename T>
    void
    join(Coefficients<Lexicographical,T,Index2D> &v,
         const CoefficientVec<T> &x, Index1D col);

} // namespace lawa

#endif // IRIS2_MY_COEFFICENTTOOLS_H