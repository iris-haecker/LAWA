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
using SymApply1D = SymmetricApply1D<T, MappedLaplace1D<T> >;

template <typename T>
    Coefficients<Lexicographical,T,Index2D>
    initRHS2D(const CompoundBasis<T>                    &V,
              const RHSIntegral2D<T>                    &rhsIntegral,
              const PreconditionerLaplace1D<T>          &P);

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

    const int                                   d = 4;
    const int                                   d_ = 6;

//
//  Setup A1
//
    Laplace1D<T>                                A1(d, d_);
    PreconditionerLaplace1D<T>                  P1(A1);
    CompressionLaplace1D<T>                     Compr1(A1);
    ParametersLaplace1D<T>                      parameters1(A1);
    MappedLaplace1D<T>                          MA1(A1,P1,Compr1);

    SymmetricApply1D<T, MappedLaplace1D<T> >    Apply1(parameters1, MA1);

//
//  Setup A2
//
    Laplace1D<T>                                A2(d, d_);
    PreconditionerLaplace1D<T>                  P2(A2);
    CompressionLaplace1D<T>                     Compr2(A2);
    ParametersLaplace1D<T>                      parameters2(A2);
    MappedLaplace1D<T>                          MA2(A2,P2,Compr2);

    SymmetricApply1D<T, MappedLaplace1D<T> >    Apply2(parameters2, MA2);

//
//  Setup Apply for the tensor Product  A = A1(x)A2
//
    typedef SymmetricApply1D<T, MappedLaplace1D<T> >        SyApply1;
    typedef SymmetricApply1D<T, MappedLaplace1D<T> >        SyApply2;
    
    typedef SymmetricTensorApply<T, SyApply1, SyApply2>     SyTensorApply;


    SyTensorApply                           Apply(Apply1, Apply2);

//
//  Setup right hand side
//
    Function2D<T>                           rhsFunc(RefSol::rhs);
    RHSIntegral2D<T>                        rhsIntegral(A1.V, A2.V, rhsFunc, 3);

    // TODO: initRHS2D should also use A2.V and P2
    Coefficients<Lexicographical,T,Index2D> F = initRHS2D(A1.V, rhsIntegral, P1);

    cout << "F = " << F << endl;

    Coefficients<Lexicographical,T,Index2D> V = Apply(F, 0.0001);

    cout << "V = A*F = " << V << endl;
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


