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

    Function<T>                             rhs_func(RefSol::rhs);
    RHSIntegral1D<T>                        rhsIntegral(A.V, rhs_func, 120);
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
Coefficients<Lexicographical,T,Index1D>
initRHS(const CompoundBasis<T>              &V,
        const RHSIntegral1D<T>              &rhsIntegral,
        const PreconditionerLaplace1D<T>    &P)
{
    Coefficients<Lexicographical,T,Index1D>     f;
    
    return f;
}
