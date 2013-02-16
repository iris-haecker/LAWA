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
    return 100*M_PI*M_PI*sin(10*M_PI*x);
}


double
sol(double x)
{
    return sin(10*M_PI*x);
}


int
main()
{
    const unsigned int d     = 2;
    const unsigned int d_    = 4;

    PrimalBasis   U(d, d_);
    PrimalBasis   V(d, d_);

    U.enforceBoundaryCondition<DirichletBC>();
    V.enforceBoundaryCondition<DirichletBC>();


    const int j0 = std::max(U.j0, V.j0);
    const int J  = 0;
    const int j1 = j0 + J;
    
    //cerr << "j0 = " << j0 << endl;

    Integral<Gauss,PrimalBasis,PrimalBasis>  integral(U, V);

    IntegerDenseVector  ku0(_(j0-1,j0+J)), kv0(_(j0-1,j0+J));
    IntegerDenseVector  KU(_(j0-1,j0+J)), KV(_(j0-1,j0+J));
    IntegerDenseVector  KU0(_(j0-1,j0+J)), KV0(_(j0-1,j0+J));
    IntegerDenseVector  KU1(_(j0-1,j0+J)), KV1(_(j0-1,j0+J));

    for (int _ju=j0-1; _ju<=j0+J; ++_ju) {
        ku0(_ju) = (_ju<j0) ? U.mra.rangeI(j0).firstIndex()
                            : U.rangeJ(_ju).firstIndex();
        KU(_ju)  = (_ju<j0) ? U.mra.rangeI(j0).length()
                            : U.rangeJ(_ju).length();
        KU0(_ju) = (_ju<j0) ? 1
                            : KU0(_ju-1) + KU(_ju-1);
        KU1(_ju) = KU0(_ju) + KU(_ju) - 1;
    }

    /*
    cerr << "ku0 = " << ku0 << endl;
    cerr << "KU  = " << KU << endl;
    cerr << "KU0 = " << KU0 << endl;
    cerr << "KU1 = " << KU1 << endl;
    */

    for (int _jv=j0-1; _jv<=j0+J; ++_jv) {
        kv0(_jv) = (_jv<j0) ? V.mra.rangeI(j0).firstIndex()
                            : V.rangeJ(_jv).firstIndex();
        KV(_jv)  = (_jv<j0) ? V.mra.rangeI(j0).length()
                            : V.rangeJ(_jv).length();
        KV0(_jv) = (_jv<j0) ? 1
                            : KV0(_jv-1) + KV(_jv-1);
        KV1(_jv) = KV0(_jv) + KV(_jv) - 1;
    }

    /*
    cerr << "kv0 = " << kv0 << endl;
    cerr << "KV  = " << KV << endl;
    cerr << "KV0 = " << KV0 << endl;
    cerr << "KV1 = " << KV1 << endl;
    */

    RealGeMatrix     A(KV1(j0+J), KU1(j0+J));

    //cerr << "A = " << A << endl;

    
    for (int _jv=j0-1; _jv<=j0+J; ++_jv) {
        int    jv     = (_jv<j0) ? j0 : _jv;
        XType  vType  = (_jv<j0) ? XBSpline : XWavelet;

        for (int _ju=j0-1; _ju<=j0+J; ++_ju) {
            int    ju     = (_ju<j0) ? j0 : _ju;
            XType  uType  = (_ju<j0) ? XBSpline : XWavelet;

            // (uType, vType)
            for (int i=0; i<KV(_jv); ++i) {
                for (int j=0; j<KU(_ju); ++j) {
                    const int ku = ku0(_ju) + j;
                    const int kv = kv0(_jv) + i;

                    double value = integral(ju, ku, uType, 1,
                                            jv, kv, vType, 1);
                    A(KV0(_jv)+i, KU0(_ju)+j) = value;
                }
            }
            //cerr << "A = " << A << endl;
            
        }
    }
    //cerr << "A = " << A << endl;


    Function<double>                fFunc(f);
    IntegralF<Gauss, PrimalBasis>   rhsIntegral(fFunc, V);
    RealDenseVector                 b(KV1(j0+J));

    for (int _jv=j0-1; _jv<=j0+J; ++_jv) {
        int    jv     = (_jv<j0) ? j0 : _jv;
        XType  vType  = (_jv<j0) ? XBSpline : XWavelet;

        for (int i=0; i<KV(_jv); ++i) {
            const int kv = kv0(_jv) + i;

            b(KV0(_jv)+i) = rhsIntegral(jv, kv, vType, 0);
        }
    }


    //cerr << "b = " << b << endl;

    const int mn = std::min(A.numRows(), A.numCols());
    IntegerDenseVector      piv(mn);
    IntegerDenseVector      xSorted(b.length());

    flens::sv(A, piv, b);

    
    // cerr << "x = " << b << endl;

//
//  Bubble sort the abs-values of b
//
    for (int i=1; i<=b.length(); ++i) {
        xSorted(i) = i;
    }

    bool swapped;

    do {
        swapped = false;
        for (int i=1; i<=b.length()-1; ++i) {
            if (std::abs(b(xSorted(i)))<std::abs(b(xSorted(i+1)))) {
                swap(xSorted(i), xSorted(i+1));
                swapped = true;
            }
        }
    } while (swapped);

    /*
    for (int i=1; i<=b.length(); ++i) {
        cerr << b(xSorted(i)) << " ";
    }
    cerr << endl;
    */


    RealDenseVector    c(KV1(j0+J));

    const int numPoints = 3000;

    for (int N=1; N<=KV1(j0+J); ++N) {
//
//      Copy the N largest coefficiencts from b into c
//
        c = 0;
        for (int i=1; i<=N; ++i) {
            c(xSorted(i)) = b(xSorted(i));
        }

        double error;
        double errorMax = 0;
        double errorSqrSum = 0;

        for (int p=0; p<=numPoints; ++p) {
            const double x = double(p)/numPoints;

            double y = 0;
            double y2 = 0;

            for (int _ju=j0-1; _ju<=j0+J; ++_ju) {
                int    ju     = (_ju<j0) ? j0 : _ju;
            
                if (_ju<j0) {
//
//                  Scaling function
//
                    for (int i=0; i<KU(_ju); ++i) {
                        const int ku = ku0(_ju) + i;
                        y += c(KU0(_ju)+i) * U.mra.phi(x, ju, ku, 0);
                    }
                } else {
//
//                  Wavelet
//
                    for (int i=0; i<KU(_ju); ++i) {
                        const int ku = ku0(_ju) + i;
                        y += c(KU0(_ju)+i) * U.psi(x, ju, ku, 0);
                    }
                }
            }
            error = sol(x) - y;
            
            if (std::abs(error)>errorMax) {
                errorMax = std::abs(error);
            }

            errorSqrSum += std::pow(error, 2);

            cout << x << " " << y << " " << error << endl;
        }
        cout << endl << endl;

        /*
        cerr << "N = " << N
             << ", errorMax = " << errorMax
             << ", errorL2 = " << sqrt(errorSqrSum/numPoints)
             << endl;
        */
        cerr << N << " "
             << errorMax << " "
             << sqrt(errorSqrSum/numPoints)
             << endl;
    }

}
