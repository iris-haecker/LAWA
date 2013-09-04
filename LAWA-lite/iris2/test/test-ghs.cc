#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

typedef double                                      T;

typedef MapMatrix<T, Index1D,
                  Laplace1D<T>,
                  CompressionLaplace1D<T>,
                  PreconditionerLaplace1D<T> >      MappedLaplace1D;

typedef SymmetricApply1D<T, MappedLaplace1D>        SymApply1D;

typedef RHS<T, Index1D,
            RHSIntegral1D<T>,
            PreconditionerLaplace1D<T> >            Rhs;

typedef GHS_ADWAV1D<T,
                    CompoundBasis<T>,
                    SymApply1D,
                    Rhs>                            GHS_Adwav;


template <typename T>
    Coefficients<Lexicographical,T,Index1D>
    initRHS(const CompoundBasis<T>              &V,
            const RHSIntegral1D<T>              &rhsIntegral,
            const PreconditionerLaplace1D<T>    &P);

int
main()
{
    typedef SolLaplace1D<T>                 RefSol;

    RefSol::setExample(1, T(1));

    Laplace1D<T>                            A(2,4);
    PreconditionerLaplace1D<T>              P(A);
    CompressionLaplace1D<T>                 Compr(A);
    ParametersLaplace1D<T>                  parameters(A);
    MappedLaplace1D                         MA(A,P,Compr);

    SymmetricApply1D<T, MappedLaplace1D>    Apply(parameters, MA);

    Function<T>                             rhsFunc(RefSol::rhs);
    RHSIntegral1D<T>                        rhsIntegral(A.V, rhsFunc, 120);
    Coefficients<Lexicographical,T,Index1D> f = initRHS(A.V, rhsIntegral, P);
    Rhs                                     F(rhsIntegral, P, f);

    GHS_Adwav                               ghs_adwav(A.V, Apply, F);

    int                                     maxNumOfIterations = 100;
    T                                       eps                = 1e-5;

    cout << "ADWAV started." << endl;
    ghs_adwav.SOLVE(f.norm(2.), eps, maxNumOfIterations, RefSol::H1norm());
    cout << "ADWAV finished." << endl;
}

template <typename T>
IndexSet<Index1D>
computeRHSLambda(const CompoundBasis<T>     &V,
                 int                        jMax)
{
    IndexSet<Index1D> ret;

    const int j0 = V.j0();
    const int k1 = V.minK(j0, XBSpline, T(0));
    const int k2 = V.maxK(j0, XBSpline, T(1));

    for (int k=k1; k<=k2; ++k) {
        ret.insert(Index1D(j0, k, XBSpline));
    }

    for (int j=j0; j<=jMax; ++j) {
        const int k1 = V.minK(j, XWavelet, T(0));
        const int k2 = V.maxK(j, XWavelet, T(1));

        for (int k=k1; k<=k2; ++k) {
            ret.insert(Index1D(j, k, XWavelet));
        }
    }

    return ret;
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
initRHS(const CompoundBasis<T>              &V,
        const RHSIntegral1D<T>              &rhsIntegral,
        const PreconditionerLaplace1D<T>    &P)
{
    typedef IndexSet<Index1D>::const_iterator const_set_it;

    Coefficients<Lexicographical,T,Index1D>     f;
    IndexSet<Index1D> RHSLambda = computeRHSLambda(V, V.j0()+1);

    for (const_set_it it=RHSLambda.begin(); it!=RHSLambda.end(); ++it) {
        f[*it] = P(*it)*rhsIntegral(*it);
    }

    return f;
}
