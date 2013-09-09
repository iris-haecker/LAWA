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
    initRHS(const CompoundBasis<T>                      &V,
            const RHSIntegral1D<T>                      &rhsIntegral,
            const PreconditionerLaplace1D<T>            &P);

template <typename T, typename SOLVER, typename MA_H, typename RHS_H>
    void
    postprocessing_H1(SOLVER                            &Solver,
                      MA_H                              &A_H1,
                      RHS_H                             &F_H1,
                      T                                 H1norm,
                      const char                        *filename);

template <typename T, typename Preconditioner>
    void
    plot(const CompoundBasis<T>                         &U,
         const Coefficients<Lexicographical,T,Index1D>  coeff,
         const Preconditioner                           &P,
         T                                              (*u)(T),
         T                                              (*du)(T),
         T                                              a,
         T                                              b,
         T                                              h,
         const char*                                    filename);

int
main()
{
    typedef SolLaplace1D<T>                 RefSol;

    RefSol::setExample(1, T(1));

    const int                               d = 4;
    const int                               d_ = 6;

    Laplace1D<T>                            A(d, d_);
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


    cout << "Postprocessing started." << endl;
    stringstream filename;
    filename << "test-ghs-W-XBSpline"
             << d << "_" << d_ << "_" << A.j0 << ".dat";
    ofstream file(filename.str().c_str());
    postprocessing_H1(ghs_adwav,
                      MA,
                      F,
                      RefSol::H1norm(),
                      filename.str().c_str());
    cout << "Postprocessing finished." << endl;


    stringstream plot_filename;
    plot_filename << "test-ghs-W-XBSpline-plot"
                  << "_" << d << "_" << d_ << "_" << A.j0 << ".dat";
    cout << "Plot of solution started." << endl;

    int numOfIterations = ghs_adwav.solutions.size();

    plot(A.U,
         ghs_adwav.solutions[numOfIterations-1],
         P,
         RefSol::u,
         RefSol::d_u,
         T(0),
         T(1),
         pow2i<T>(-5),
         plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;


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

template <typename T, typename SOLVER, typename MA_H, typename RHS_H>
void
postprocessing_H1(SOLVER& Solver, MA_H &A_H1, RHS_H &F_H1, T H1norm, const char* filename)
{

    Coefficients<Lexicographical,T,Index1D> u;
    std::ofstream file(filename);

    for (int i=0; i<int(Solver.solutions.size()); ++i) {
        u = Solver.solutions[i];
        T ErrorH1Norm = computeErrorInH1Norm(A_H1, F_H1, u, H1norm);
        file      << supp(u).size() << " " << Solver.linsolve_iterations[i] << " "
                  << Solver.times[i] << " " << Solver.residuals[i] << " "
                  << ErrorH1Norm << std::endl;
        std::cerr << supp(u).size() << " " << Solver.linsolve_iterations[i] << " "
                  << Solver.times[i] << " " << Solver.residuals[i] << " "
                  << ErrorH1Norm << std::endl;
    }
    file.close();
}


template <typename T, typename Preconditioner>
void
plot(const CompoundBasis<T>                         &U,
     const Coefficients<Lexicographical,T,Index1D>  coeff,
     const Preconditioner                           &P,
     T                                              (*u)(T),
     T                                              (*du)(T),
     T                                              a,
     T                                              b,
     T                                              h,
     const char*                                    filename)
{
    std::cerr << "hello" << std::endl;
    
    typedef Coefficients<Lexicographical,T,Index1D >    Coeff;
    typedef typename Coeff::const_iterator              coeff_it;

    std::ofstream plotfile(filename);
    for (T x=a; x<=b; x+=h) {
        T appr=0., d_appr = 0.0;
        T exact= u(x);
        T d_exact= du(x);
        for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int j = (*it).first.j, k = (*it).first.k;
            T coeff = (*it).second, prec = P((*it).first);

            appr   += prec * coeff * U((*it).first.xtype, j, k, x, 0);
            d_appr += prec * coeff * U((*it).first.xtype, j, k, x, 1);

        }
        plotfile << x << " " << exact << " " << d_exact << " " << appr << " " << d_appr << std::endl;
    }
    plotfile.close();
}

