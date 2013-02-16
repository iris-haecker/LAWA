#include <cmath>
#include <fstream>
#include <iostream>

#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
typedef flens::DenseVector<Array<double> >               RealDenseVector;
typedef flens::DenseVector<Array<int> >                  IntegerDenseVector;


typedef Basis<double,Primal,Interval,Dijkema>            PrimalBasis;


const int    k  = 1;
const double u0 = 0;

double
sol(double t)
{
    //return 1+sin(2*M_PI*t);
    return sin(2*M_PI*k*t) + u0;
}

double
g(double t)
{
    //return 1 + sin(2*M_PI*t) + 2*M_PI*cos(2*M_PI*t);
    return sin(2*M_PI*k*t) + 2*M_PI*k*cos(2*M_PI*k*t) + u0;

    /*
    if (t<1.0/3) {
        return 1;
    }
    return 20;
    */
}

int
main()
{
    const unsigned int d     = 3;
    const unsigned int d_    = 5;

    PrimalBasis   U0(d, d_);
    PrimalBasis   U1(d, d_);

    PrimalBasis   V0(d, d_);
    PrimalBasis   V1(d, d_);

    //U0.enforceBoundaryCondition<DirichletBC>();
    V0.enforceBoundaryCondition<DirichletBC>();

    const int j0 = std::max(std::max(std::max(U0.j0, V0.j0), U1.j0), V1.j0);

    const int J  = 1;
    const int j1 = j0 + J;

    //cerr << "j0-1 = " << (j0-1) << endl;
    //cerr << "j1   = " << j1 << endl;

    Integral<Gauss,PrimalBasis,PrimalBasis>  integral00(U0, V0);
    Integral<Gauss,PrimalBasis,PrimalBasis>  integral01(U0, V1);
    Integral<Gauss,PrimalBasis,PrimalBasis>  integral10(U1, V0);
    Integral<Gauss,PrimalBasis,PrimalBasis>  integral11(U1, V1);


    IntegerDenseVector  ku0(_(j0-1,j0+J)), kv0(_(j0-1,j0+J));
    IntegerDenseVector  ku1(_(j0-1,j0+J)), kv1(_(j0-1,j0+J));
    IntegerDenseVector  kum(_(j0-1,j0+J)), kvm(_(j0-1,j0+J));

    IntegerDenseVector  KU(_(j0-1,j0+J)), KV(_(j0-1,j0+J));
    IntegerDenseVector  KU0(_(j0-1,j0+J)), KV0(_(j0-1,j0+J));
    IntegerDenseVector  KU1(_(j0-1,j0+J)), KV1(_(j0-1,j0+J));

    for (int _ju=j0-1; _ju<=j0+J; ++_ju) {
        ku0(_ju) = (_ju<j0) ? U0.mra.rangeI(j0).firstIndex()
                            : U0.rangeJ(_ju).firstIndex();
        ku1(_ju) = (_ju<j0) ? U1.mra.rangeI(j0).lastIndex()
                            : U1.rangeJ(_ju).lastIndex();
        kum(_ju) = (ku1(_ju)+ku0(_ju))/2;

        KU(_ju)  = ku1(_ju)-ku0(_ju)+1;
        KU0(_ju) = (_ju<j0) ? 1
                            : KU0(_ju-1) + KU(_ju-1);
        KU1(_ju) = KU0(_ju) + KU(_ju) - 1;
    }

    cerr << "ku0 = " << ku0 << endl;
    cerr << "ku1 = " << ku1 << endl;
    cerr << "KU  = " << KU << endl;
    cerr << "KU0 = " << KU0 << endl;
    cerr << "KU1 = " << KU1 << endl;

    for (int _jv=j0-1; _jv<=j0+J; ++_jv) {
        kv0(_jv) = (_jv<j0) ? V1.mra.rangeI(j0).firstIndex()
                            : V1.rangeJ(_jv).firstIndex();
        kv1(_jv) = (_jv<j0) ? V0.mra.rangeI(j0).lastIndex()
                            : V0.rangeJ(_jv).lastIndex();
        kvm(_jv) = (kv0(_jv)+kv1(_jv))/2;

        KV(_jv)  = kv1(_jv)-kv0(_jv)+1;
        KV0(_jv) = (_jv<j0) ? 1
                            : KV0(_jv-1) + KV(_jv-1);
        KV1(_jv) = KV0(_jv) + KV(_jv) - 1;
    }

    /*
    cerr << "kv0 = " << kv0 << endl;
    cerr << "kv1 = " << kv1 << endl;
    cerr << "kvm = " << kvm << endl;
    cerr << "KV  = " << KV << endl;
    cerr << "KV0 = " << KV0 << endl;
    cerr << "KV1 = " << KV1 << endl;
    */
//
//  Setup stiffness matrix
//

    RealGeMatrix     A(KV1(j0+J), KU1(j0+J));

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

                    double value;

                    if (ku<kum(_ju)) {
                        if (kv<kvm(_jv)) {
                            value = -integral01(ju, ku, uType, 0,
                                                jv, kv, vType, 1)
                                    +integral01(ju, ku, uType, 0,
                                                jv, kv, vType, 0);
                        } else {
                            value = -integral00(ju, ku, uType, 0,
                                                jv, kv, vType, 1)
                                    +integral00(ju, ku, uType, 0,
                                                jv, kv, vType, 0);
                        }
                    } else {
                        if (kv<kvm(_jv)) {
                            value = -integral11(ju, ku, uType, 0,
                                                jv, kv, vType, 1)
                                    +integral11(ju, ku, uType, 0,
                                                jv, kv, vType, 0);
                        } else {
                            value = -integral10(ju, ku, uType, 0,
                                                jv, kv, vType, 1)
                                    +integral10(ju, ku, uType, 0,
                                                jv, kv, vType, 0);
                        }
                    }
                    A(KV0(_jv)+i, KU0(_ju)+j) = value;
                }
            }
        }
    }
    cerr << "A = " << A << endl;

//
//  Setup right-hand side
//
    Function<double>                gFunc(g);
    IntegralF<Gauss, PrimalBasis>   rhsIntegral0(gFunc, V0);
    IntegralF<Gauss, PrimalBasis>   rhsIntegral1(gFunc, V1);
    RealDenseVector                 b(KV1(j0+J));

    for (int _jv=j0-1; _jv<=j0+J; ++_jv) {
        int    jv     = (_jv<j0) ? j0 : _jv;
        XType  vType  = (_jv<j0) ? XBSpline : XWavelet;

        for (int i=0; i<KV(_jv); ++i) {
            const int kv = kv0(_jv) + i;

            if (kv<=kvm(_jv)) {
                b(KV0(_jv)+i) = rhsIntegral1(jv, kv, vType, 0);
            } else {
                b(KV0(_jv)+i) = rhsIntegral0(jv, kv, vType, 0);
            }
        }
    }
    cerr << "-> b = " << b << endl;

//
//  Solve least square problem Ax=b (by solving A'*A x = A'*b)
//

    RealGeMatrix  AtA;
    blas::mm(Trans, NoTrans, 1.0, A, A, 0.0, AtA);


    RealGeMatrix  _AtA = AtA;


    RealDenseVector Atb;
    blas::mv(Trans, 1.0, A, b, 0.0, Atb);
    Atb.engine().changeIndexBase(1);

    cerr << "AtA = " << AtA << endl;
    cerr << "Atb = " << Atb << endl;

    RealDenseVector  _Atb = Atb;

    IntegerDenseVector      piv(AtA.numRows());

    flens::sv(AtA, piv, Atb);
    cerr << "AtA = " << AtA << endl;
    cerr << "piv = " << piv << endl;
    cerr << "x = " << Atb << endl;

    RealDenseVector  x(Atb.length());
    x = Atb;

    assert(x.firstIndex()==1);
/*
//
//  Bubble sort the abs-values of x
//
    RealDenseVector xSorted(x.length());
    for (int i=1; i<=x.length(); ++i) {
        xSorted(i) = i;
    }

    bool swapped;

    do {
        swapped = false;
        for (int i=1; i<=b.length()-1; ++i) {
            if (std::abs(x(xSorted(i)))<std::abs(x(xSorted(i+1)))) {
                swap(xSorted(i), xSorted(i+1));
                swapped = true;
            }
        }
    } while (swapped);


    RealDenseVector  c(x.length());
    
    int maxN = std::min(140, x.length());
    
    //cerr << "maxN = " << maxN << endl;


//
//  Plot coefficients
//

    for (int N=1; N<=maxN; ++N) {
//
//      Copy the N largest coefficiencts from x into c
//
        c = 0;
        for (int i=1; i<=N; ++i) {
            c(xSorted(i)) = x(xSorted(i));
        }

        cerr << "N = " << N << endl;
        cerr << "c = " << c << endl;

        Coefficients<AbsoluteValue, double, Index1D >  coeff;
        Index1D                                        index;
        double                                         coeffValue;

        for (int _ju=j0-1; _ju<=j0+J; ++_ju) {
            int ju     = (_ju<j0) ? j0 : _ju;

            if (_ju<j0) {
//
//              Scaling function
//
                for (int i=0; i<KU(_ju); ++i) {
                    const int ku = KU0(_ju) + i;
                    const int k  = ku0(_ju)+i;

                    if (std::abs(c(ku))>0) {
                        Index1D     index(ju, ku0(_ju)+i, XBSpline);

                        coeff.insert(pair<double, Index1D>(c(ku), index));

                        cerr << "insert: "
                             << "N = " << N
                             << ", ku = " << ku
                             << ", k = " << k
                             << ", j = " << index.j
                             << ", k = " << index.k
                             << ", type = " << index.xtype
                             << ", coeffValue = " << c(ku)
                             << endl;
                    }

                }
            } else {
//
//              Wavelet
//
                for (int i=0; i<KU(_ju); ++i) {
                    const int ku = KU0(_ju) + i;
                    const int k  = ku0(_ju)+i;

                    if (std::abs(c(ku))>0) {
                        Index1D     index(ju, ku0(_ju)+i, XWavelet);

                        coeff.insert(pair<double, Index1D>(c(ku), index));

                        cerr << "insert: "
                             << "N = " << N
                             << ", ku = " << ku
                             << ", k = " << k
                             << ", j = " << index.j
                             << ", k = " << index.k
                             << ", type = " << index.xtype
                             << ", coeffValue = " << c(ku)
                             << endl;
                    }
                }
            }

        }

        stringstream filename;
        filename << "ode_nterm_coeff" << setw(3) << setfill('0') << N;

        plotCoeff(coeff, U0, j0, J, filename.str().c_str());
    }



//
//  Plotting functions
//
    std::ofstream errorData("ode_nterm_error.dat");

    std::ofstream errorPlot("ode_nterm_error.gps");
    errorPlot << "reset" << std::endl;
    errorPlot << "set terminal png; set output 'ode_nterm_error.png'"
              << std::endl;
    errorPlot << "set palette color; set colorbox vertical" << std::endl;
    errorPlot << "plot \"ode_nterm_error.dat\" "
              << "using 1:2 with lines, "
              << "\"ode_nterm_error.dat\" "
              << "using 1:3 with lines";
    errorPlot.close();


    for (int N=1; N<=maxN; ++N) {
        stringstream filenameData, filenameGps, filenamePng;

        filenameData << "ode_nterm_plot"
                     << setw(3) << setfill('0') << N
                     << ".dat";

        filenameGps << "ode_nterm_plot"
                    << setw(3) << setfill('0') << N
                    << ".gps";

        filenamePng << "ode_nterm_plot"
                    << setw(3) << setfill('0') << N
                    << ".png";


        std::ofstream plotData(filenameData.str().c_str());
        std::ofstream plotGps(filenameGps.str().c_str());

        plotGps << "reset" << std::endl;
        plotGps << "set terminal png; set output '"
                << filenamePng.str()
                << "'" << std::endl;
        plotGps << "set palette color; set colorbox vertical" << std::endl;
        plotGps << "set xrange[0:1]" << std::endl;
        plotGps << "set yrange[-2:2]" << std::endl;
        plotGps << "plot \"" << filenameData.str() << "\" "
                << "using 1:2 with lines, "
                << "\"" << filenameData.str() << "\" "
                << "using 1:4 with lines";
        plotGps.close();

//
//      Copy the N largest coefficiencts from x into c
//
        c = 0;
        for (int i=1; i<=N; ++i) {
            c(xSorted(i)) = x(xSorted(i));
        }
        
        cerr << "N = " << N << endl;
        cerr << "c = " << c << endl;

        Coefficients<AbsoluteValue, double, Index1D >  coeff;
        Index1D                                        index;
        double                                         coeffValue;

        double solution;
        double error;
        double errorMax = 0;
        double errorSqrSum = 0;

        const int numPoints = 5000;
        for (int p=0; p<=numPoints; ++p) {
            const double x = double(p)/numPoints;

            double y = 0;

            for (int _ju=j0-1; _ju<=j0+J; ++_ju) {
                int    ju     = (_ju<j0) ? j0 : _ju;

                if (_ju<j0) {
//
//                  Scaling function
//
                    for (int i=0; i<KU(_ju); ++i) {
                        const int ku = ku0(_ju) + i;


                        double value = 0;
                        if (ku<=kum(_ju)) {
                            value = c(KU0(_ju)+i) * U0.mra.phi(x, ju, ku, 0);
                        } else {
                            value = c(KU0(_ju)+i) * U1.mra.phi(x, ju, ku, 0);
                        }
                        y += value;
                    }
                } else {
//
//                  Wavelet
//
                    for (int i=0; i<KU(_ju); ++i) {
                        const int ku = ku0(_ju) + i;

                        double value = 0;
                        if (ku<=kum(_ju)) {
                            value = c(KU0(_ju)+i) * U0.psi(x, ju, ku, 0);
                        } else {
                            value = c(KU0(_ju)+i) * U1.psi(x, ju, ku, 0);
                        }
                        y += value;
                        coeffValue = c(KU0(_ju)+i);
                    }
                }
                
            }
            solution = sol(x);
            error    = solution - y;
            
            if (std::abs(error)>errorMax) {
                errorMax = std::abs(error);
            }

            errorSqrSum += std::pow(error, 2);

            plotData << x << " "
                     << y << " "
                     << error << " "
                     << solution << endl;
         }

         errorData << N << " "
                   << errorMax << " "
                   << sqrt(errorSqrSum/numPoints)
                   << endl;
         cerr << N << " "
              << errorMax << " "
              << sqrt(errorSqrSum/numPoints)
              << endl;


    }
*/
}