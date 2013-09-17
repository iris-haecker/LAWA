#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

template <typename T>
using MappedLaplace1D = MapMatrix<T,
                                  Index1D,
                                  Laplace1D<T>,
                                  CompressionLaplace1D<T>,
                                  PreconditionerLaplace1D<T> >;

template <typename T>
using Rhs = StaticRHS<T, Index1D, PreconditionerLaplace1D<T> >;

template <typename T>
using SymApply1D = SymmetricApply1D<T, MappedLaplace1D<T> >;

template <typename T>
using GHS_Adwav =  GHS_ADWAV1D<T,
                               CompoundBasis<T>,
                               SymApply1D<T>,
                               Rhs<T> >;

template <typename T>
    Coefficients<Lexicographical,T,Index2D>
    initRHS2D(const CompoundBasis<T>                    &V,
              const RHSIntegral2D<T>                    &rhsIntegral,
              const PreconditionerLaplace1D<T>          &P);

template <typename T>
using CoefficientVec = Coefficients<Lexicographical,T,Index1D>;

template <typename T>
using CoefficientRows = map<Index1D,
                            CoefficientVec<T>,
                            lt<Lexicographical, Index1D> >;

template <typename T>
using CoefficientCols = map<Index1D,
                            CoefficientVec<T>,
                            lt<Lexicographical, Index1D> >;

template <typename T>
    CoefficientCols<T>
    splitCols(const Coefficients<Lexicographical,T,Index2D> &v);

template <typename T>
    CoefficientRows<T>
    splitRows(const CoefficientCols<T> &colMap);

template <typename T>
    Coefficients<Lexicographical,T,Index2D>
    joinRows(const CoefficientRows<T> &rowMap);

template <typename T>
    Coefficients<Lexicographical,T,Index1D>
    solveLaplace1D(const PreconditionerLaplace1D<T>                 &P, 
                   const Laplace1D<T>                               &A,
                   const Coefficients<Lexicographical,T,Index1D>    &f,
                   T                                                eps);

int
main()
{
    typedef double                          T;
    typedef SolLaplace1D_M2<T>              RefSol;

    RefSol::setExample(1, T(1));

    const int                               d = 4;
    const int                               d_ = 6;

    Laplace1D<T>                            A(d, d_);
    PreconditionerLaplace1D<T>              P(A);

    Function2D<T>                           rhsFunc(RefSol::rhs);
    RHSIntegral2D<T>                        rhsIntegral(A.V, A.V, rhsFunc, 3);

    Coefficients<Lexicographical,T,Index2D> F = initRHS2D(A.V, rhsIntegral, P);


    CoefficientCols<T> colsF = splitCols(F);
    CoefficientCols<T> colsV;
    CoefficientRows<T> rowsU;

    for (auto it=colsF.begin(); it!=colsF.end(); ++it) {
        Index1D                  col = it->first;
        const CoefficientVec<T>  &f  = it->second;

        cout << "col = " << col << endl;

        colsV[col] = solveLaplace1D(P, A, f, 1e-7);
    }
    
    cout << endl << endl << endl << endl << endl << endl << endl
         << "--------------------------------------------" << endl;

    CoefficientRows<T> rowsV = splitRows(colsV);
    for (auto it=rowsV.begin(); it!=rowsV.end(); ++it) {
        Index1D                  row = it->first;
        const CoefficientVec<T>  &vt = it->second;

        cout << "row = " << row << endl;
        cout << "cols = " << supp(vt) << endl;

        rowsU[row] = solveLaplace1D(P, A, vt, 1e-5);
    }

    Coefficients<Lexicographical,T,Index2D>  U = joinRows(rowsU);

    cout << "U = " << U << endl;
}

template <typename T>
void
computeRHSLambda2D_(IndexSet<Index2D>          &ret,
                    int                        j_,
                    int                        k_,
                    XType                      e_,
                    const CompoundBasis<T>     &V,
                    int                        jMax)
{
    const int j0 = V.j0();
    const int k1 = V.minK(j0, XBSpline, T(0));
    const int k2 = V.maxK(j0, XBSpline, T(1));

    for (int k=k1; k<=k2; ++k) {
        ret.insert(Index2D(Index1D(j_, k_, e_), Index1D(j0, k, XBSpline)));
    }

    for (int j=j0; j<=jMax; ++j) {
        const int k1 = V.minK(j, XWavelet, T(0));
        const int k2 = V.maxK(j, XWavelet, T(1));

        for (int k=k1; k<=k2; ++k) {
            ret.insert(Index2D(Index1D(j_, k_, e_), Index1D(j, k, XWavelet)));
        }
    }
}

template <typename T>
IndexSet<Index2D>
computeRHSLambda2D(const CompoundBasis<T>     &V,
                   int                        jMax)
{
    IndexSet<Index2D> ret;

    const int j0 = V.j0();
    const int k1 = V.minK(j0, XBSpline, T(0));
    const int k2 = V.maxK(j0, XBSpline, T(1));

    for (int k=k1; k<=k2; ++k) {
        computeRHSLambda2D_(ret, j0, k, XBSpline, V, jMax);
    }

    for (int j=j0; j<=jMax; ++j) {
        const int k1 = V.minK(j, XWavelet, T(0));
        const int k2 = V.maxK(j, XWavelet, T(1));

        for (int k=k1; k<=k2; ++k) {
            computeRHSLambda2D_(ret, j, k, XWavelet, V, jMax);
        }
    }

    return ret;
}

template <typename T>
Coefficients<Lexicographical,T,Index2D>
initRHS2D(const CompoundBasis<T>                &V,
          const RHSIntegral2D<T>                &rhsIntegral,
          const PreconditionerLaplace1D<T>      &P)
{
    typedef IndexSet<Index2D>::const_iterator   const_set_it;

    Coefficients<Lexicographical,T,Index2D>     f;

    IndexSet<Index2D> RHSLambda = computeRHSLambda2D(V, V.j0());

    for (const_set_it it=RHSLambda.begin(); it!=RHSLambda.end(); ++it) {
        f[*it] = P(it->index1)*P(it->index2)*rhsIntegral(*it);
    }

    return f;
}

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
Coefficients<Lexicographical,T,Index1D>
solveLaplace1D(const PreconditionerLaplace1D<T>                 &P, 
               const Laplace1D<T>                               &A,
               const Coefficients<Lexicographical,T,Index1D>    &f,
               T                                                eps)
{
    typedef SolLaplace1D<T>                     RefSol;
    RefSol::setExample(1, T(1));

    CompressionLaplace1D<T>                     Compr(A);
    ParametersLaplace1D<T>                      parameters(A);
    MappedLaplace1D<T>                          MA(A,P,Compr);

    Rhs<T>                                      F(P, f);

    SymmetricApply1D<T, MappedLaplace1D<T> >    Apply(parameters, MA);

    GHS_Adwav<T>                                ghs_adwav(A.V, Apply, F);

    int                                         maxNumOfIterations = 100;

    cout << "ADWAV started." << endl;
    ghs_adwav.SOLVE(f.norm(2.), eps, maxNumOfIterations, RefSol::H1norm());
    cout << "ADWAV finished." << endl;

    int numOfIterations = ghs_adwav.solutions.size();

    Coefficients<Lexicographical,T,Index1D>  solution;

    solution = ghs_adwav.solutions[numOfIterations-1];

    return solution;
}

