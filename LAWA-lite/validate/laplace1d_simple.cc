#include <cmath>
#include <fstream>
#include <iostream>

#include <lawa/lawa.h>

#ifdef NEW_LAWA
#include <lawa/lawa.tcc>
#endif

using namespace lawa;
using namespace std;

typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
typedef flens::DenseVector<Array<double> >               RealDenseVector;
typedef flens::DenseVector<Array<int> >                  IntegerDenseVector;


typedef Basis<double,Primal,Interval,Dijkema>            PrimalBasis;
typedef Basis<double,Dual,Interval,Dijkema>              DualBasis;


double
f(double x)
{
    return M_PI*M_PI*sin(M_PI*x);
}


int
main()
{
    const unsigned int d     = 2;
    const unsigned int d_    = 4;

    PrimalBasis   U(d, d_);
    PrimalBasis   V(d, d_);

#ifdef NEW_LAWA
    U.enforceBoundaryConditions<NoBC, DirichletBC>();
    V.enforceBoundaryConditions<NoBC, DirichletBC>();
#else
    U.enforceBoundaryCondition<DirichletBC>();
    V.enforceBoundaryCondition<DirichletBC>();
#endif// NEW_LAWA

    const int j0 = std::max(U.j0, V.j0);
    //const int J  = j0+1;

    cerr << "j0 = " << j0 << endl;

    Integral<Gauss,PrimalBasis,PrimalBasis>  integral(U, V);

    const int sk0 = U.mra.rangeI(j0).firstIndex();
    const int sk1 = U.mra.rangeI(j0).lastIndex();
    const int sk_ = U.mra.rangeI(j0).length();

    const int wk0 = U.rangeJ(j0).firstIndex();
    const int wk1 = U.rangeJ(j0).lastIndex();
    const int wk_ = U.rangeJ(j0).length();

    const int SK0 = V.mra.rangeI(j0).firstIndex();
    const int SK1 = V.mra.rangeI(j0).lastIndex();
    const int SK_ = V.mra.rangeI(j0).length();

    const int WK0 = V.rangeJ(j0).firstIndex();
    const int WK1 = V.rangeJ(j0).lastIndex();
    const int WK_ = V.rangeJ(j0).length();

    RealGeMatrix     A(SK_+WK_, sk_+wk_);

    // (Spline, Spline)
    for (int i=1; i<=SK_; ++i) {
        for (int j=1; j<=sk_; ++j) {
            const int k = sk0+j-1;
            const int K = SK0+i-1;

            double value = integral(j0, k, XBSpline, 1,
                                    j0, K, XBSpline, 1);
            A(i, j) = value;
        }
    }
    cerr << "1) A = " << A << endl;

    // (Wavelet, Spline)
    for (int i=1; i<=SK_; ++i) {
        for (int j=1; j<=wk_; ++j) {
            const int k = wk0+j-1;
            const int K = SK0+i-1;

            double value = integral(j0, k, XWavelet, 1,
                                    j0, K, XBSpline, 1);
            A(i, sk_+j) = value;
        }
    }
    cerr << "2) A = " << A << endl;

    // (Spline, Wavelet)
    for (int i=1; i<=WK_; ++i) {
        for (int j=1; j<=sk_; ++j) {
            const int k = sk0+j-1;
            const int K = WK0+i-1;

            double value = integral(j0, k, XBSpline, 1,
                                    j0, K, XWavelet, 1);
            A(SK_+i, j) = value;
        }
    }
    cerr << "3) A = " << A << endl;

    // (Wavelet, Wavelet)
    for (int i=1; i<=WK_; ++i) {
        for (int j=1; j<=wk_; ++j) {
            const int k = wk0+j-1;
            const int K = WK0+i-1;

            double value = integral(j0, k, XWavelet, 1,
                                    j0, K, XWavelet, 1);
            A(SK_+i, sk_+j) = value;
        }
    }
    cerr << "4) A = " << A << endl;

    cerr << "A = " << A << endl;

    Function<double>                fFunc(f);
    IntegralF<Gauss, PrimalBasis>   rhsIntegral(fFunc, V);
    RealDenseVector                 b(SK_+WK_);

    for (int i=1; i<=SK_; ++i) {
        const int K = SK0+i-1;
        b(i) = rhsIntegral(j0, K, XBSpline, 0);
    }

    for (int i=1; i<=WK_; ++i) {
        const int K = WK0+i-1;

        b(SK_+i) = rhsIntegral(j0, K, XWavelet, 0);
    }

    cerr << "b = " << b << endl;

    const int mn = std::min(A.numRows(), A.numCols());
    IntegerDenseVector      piv(mn);

    flens::sv(A, piv, b);
    cerr << "x = " << b << endl;

    const int N = 3000;
    for (int p=0; p<=N; ++p) {
        const double x = double(p)/N;

        double y = 0;
        double y2 = 0;

        for (int i=1; i<=sk_; ++i) {
            const int k = sk0+i-1;

            y += b(i)*U.mra.phi(x, j0, k, 0);
        }

        for (int i=1; i<=wk_; ++i) {
            const int k = wk0+i-1;

            y2 += b(i+wk_-1)*U.psi(x, j0, k, 0);
        }
        
        cout << x << " " << y << " " << y2 << " " << y+y2 << endl;
    }

}
