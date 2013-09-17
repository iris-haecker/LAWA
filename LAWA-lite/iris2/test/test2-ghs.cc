#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

typedef double                                      T;

typedef MapMatrix<T, Index1D,
                  Laplace1D<T>,
                  CompressionLaplace1D<T>,
                  PreconditionerLaplace1D<T> >              MappedLaplace1D;

typedef SymmetricApply1D<T, MappedLaplace1D>                SymApply1D;

typedef StaticRHS<T, Index1D, PreconditionerLaplace1D<T> >  Rhs;

typedef GHS_ADWAV1D<T,
                    CompoundBasis<T>,
                    SymApply1D,
                    Rhs>                                    GHS_Adwav;


template <typename T>
    Coefficients<Lexicographical,T,Index1D>
    initRHS(const CompoundBasis<T>                      &V,
            const RHSIntegral1D<T>                      &rhsIntegral,
            const PreconditionerLaplace1D<T>            &P);

template <typename T, typename Preconditioner>
    void
    plot(const Preconditioner                           &P,
         const CompoundBasis<T>                         &U,
         const Coefficients<Lexicographical,T,Index1D>  &u,
         T                                              a,
         T                                              b,
         T                                              h,
         const char*                                    filename);

template <typename T>
     Coefficients<Lexicographical,T,Index1D>
     solveLaplace1D(const PreconditionerLaplace1D<T>                 &P, 
                    const Laplace1D<T>                               &A,
                    const Coefficients<Lexicographical,T,Index1D>    &f);

int
main()
{
    typedef SolLaplace1D<T>                 RefSol;

    RefSol::setExample(1, T(1));

    const int                               d = 4;
    const int                               d_ = 6;

    Laplace1D<T>                            A(d, d_);
    PreconditionerLaplace1D<T>              P(A);

    Function<T>                             rhsFunc(RefSol::rhs);
    RHSIntegral1D<T>                        rhsIntegral(A.V, rhsFunc, 120);

    Coefficients<Lexicographical,T,Index1D> f, u;

    f = initRHS(A.V, rhsIntegral, P);
    u = solveLaplace1D(P, A, f);

    stringstream plot_filename;
    plot_filename << "test2-ghs-plot"
                  << "_" << d << "_" << d_ << "_" << A.j0 << ".dat";
    cout << "Plot of solution started." << endl;
    plot(P, A.U, u, T(0), T(1), pow2i<T>(-10), plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;

    cout << "u.size() = " << u.size() << endl;
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
    IndexSet<Index1D> RHSLambda = computeRHSLambda(V, V.j0()+6);

    for (const_set_it it=RHSLambda.begin(); it!=RHSLambda.end(); ++it) {
        f[*it] = P(*it)*rhsIntegral(*it);
    }

    return f;
}

template <typename T, typename Preconditioner>
void
plot(const Preconditioner                           &P,
     const CompoundBasis<T>                         &U,
     const Coefficients<Lexicographical,T,Index1D>  &u,
     T                                              a,
     T                                              b,
     T                                              h,
     const char*                                    filename)
{
    typedef Coefficients<Lexicographical,T,Index1D >    Coeff;
    typedef typename Coeff::const_iterator              coeff_it;

    std::ofstream plotfile(filename);
    for (T x=a; x<=b; x+=h) {
        T appr=T(0);
        for (coeff_it it = u.begin(); it != u.end(); ++it) {
            int j = (*it).first.j;
            int k = (*it).first.k;

            T coeff = (*it).second;
            T prec = P((*it).first);

            appr   += prec * coeff * U((*it).first.xtype, j, k, x, 0);
        }
        plotfile << x << " " << appr << std::endl;
    }
    plotfile.close();
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
solveLaplace1D(const PreconditionerLaplace1D<T>                 &P, 
               const Laplace1D<T>                               &A,
               const Coefficients<Lexicographical,T,Index1D>    &f)
{
    typedef SolLaplace1D<T>                 RefSol;

    CompressionLaplace1D<T>                 Compr(A);
    ParametersLaplace1D<T>                  parameters(A);
    MappedLaplace1D                         MA(A,P,Compr);

    Rhs                                     F(P, f);

    SymmetricApply1D<T, MappedLaplace1D>    Apply(parameters, MA);

    GHS_Adwav                               ghs_adwav(A.V, Apply, F);

    int                                     maxNumOfIterations = 100;
    T                                       eps                = 1e-5;

    cout << "ADWAV started." << endl;
    ghs_adwav.SOLVE(f.norm(2.), eps, maxNumOfIterations, RefSol::H1norm());
    cout << "ADWAV finished." << endl;

    int numOfIterations = ghs_adwav.solutions.size();

    Coefficients<Lexicographical,T,Index1D>  solution;

    solution = ghs_adwav.solutions[numOfIterations-1];

    return solution;
}

