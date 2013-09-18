#ifndef IRIS2_MY_COEFFICENTTOOLS_TCC
#define IRIS2_MY_COEFFICENTTOOLS_TCC 1

#include <iris2/iris2.h>

namespace lawa {

template <typename T>
CoefficientCols<T>
splitCols(const Coefficients<Lexicographical,T,Index2D> &v)
{
    CoefficientCols<T> colMap;

    for (auto it=v.begin(); it!=v.end(); ++it) {
        Index2D index = it->first;
        T       value = it->second;

        Index1D row = index.index1;
        Index1D col = index.index2;

        colMap[col][row] = value;
    }

    return colMap;
}

template <typename T>
CoefficientRows<T>
splitRows(const CoefficientCols<T> &colMap)
{
    CoefficientRows<T> rowMap;

    for (auto colIt=colMap.begin(); colIt!=colMap.end(); ++colIt) {
        Index1D col = colIt->first;

        for (auto rowIt=colIt->second.begin(); rowIt!=colIt->second.end(); ++rowIt) {
            Index1D row   = rowIt->first;
            T       value = rowIt->second;

            rowMap[row][col] = value;
        }
    }
    return rowMap;
}

template <typename T>
CoefficientRows<T>
splitRows(const Coefficients<Lexicographical,T,Index2D> &v)
{
    CoefficientRows<T> rowMap;

    for (auto it=v.begin(); it!=v.end(); ++it) {
        Index2D index = it->first;
        T       value = it->second;

        Index1D row = index.index1;
        Index1D col = index.index2;

        rowMap[row][col] = value;
    }

    return rowMap;
}

template <typename T>
Coefficients<Lexicographical,T,Index2D>
joinRows(const CoefficientRows<T> &rowMap)
{
    Coefficients<Lexicographical,T,Index2D>     v;

    for (auto rowIt=rowMap.begin(); rowIt!=rowMap.end(); ++rowIt) {
        Index1D row = rowIt->first;

        for (auto colIt=rowIt->second.begin(); colIt!=rowIt->second.end(); ++colIt) {
            Index1D col   = colIt->first;
            T       value = colIt->second;

            Index2D index(row, col);
            v[index] = value;
        }
    }
    return v;
}

template <typename T>
void
join(Coefficients<Lexicographical,T,Index2D> &v,
     Index1D row, const CoefficientVec<T> &x)
{
    for (auto itX=x.begin(); itX!=x.end(); ++itX) {
        Index1D col   = itX->first;
        T       value = itX->second;

        Index2D index(row, col);
        v[index] = value;
    }
}

template <typename T>
void
join(Coefficients<Lexicographical,T,Index2D> &v,
     const CoefficientVec<T> &x, Index1D col)
{
    for (auto itX=x.begin(); itX!=x.end(); ++itX) {
        Index1D row   = itX->first;
        T       value = itX->second;

        Index2D index(row, col);
        v[index] = value;
    }
}

} // namespace lawa

#endif // IRIS2_MY_COEFFICENTTOOLS_TCC