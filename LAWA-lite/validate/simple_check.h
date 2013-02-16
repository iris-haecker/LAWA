
double
amax(const RealGeMatrix &A)
{
    double result = 0;
    for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
        for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
            if (std::abs(A(i,j))>result) {
                result = std::abs(A(i,j));
            }
        }
    }
    return result;
}

double
amax(const RealDenseVector &x)
{
    double result = 0;
    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        if (std::abs(x(i))>result) {
            result = std::abs(x(i));
        }
    }
    return result;
}



void
simpleCheck(const RealGeMatrix    &_A,
            const RealGeMatrix    &_AtA,
            const RealDenseVector &_b,
            const RealDenseVector &_Atb,
            const RealDenseVector &_c)
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

        cerr << "j0 = " << j0 << endl;

        Integral<Gauss,PrimalBasis,PrimalBasis>  integral00(U0, V0);
        Integral<Gauss,PrimalBasis,PrimalBasis>  integral01(U0, V1);
        Integral<Gauss,PrimalBasis,PrimalBasis>  integral10(U1, V0);
        Integral<Gauss,PrimalBasis,PrimalBasis>  integral11(U1, V1);

        const int sk0 = U0.mra.rangeI(j0).firstIndex();
        const int sk1 = U1.mra.rangeI(j0).lastIndex();
        const int skm = (sk1+sk0)/2;
        const int sk_ = sk1-sk0+1;

        const int wk0 = U0.rangeJ(j0).firstIndex();
        const int wk1 = U1.rangeJ(j0).lastIndex();
        const int wkm = (wk1+wk0)/2;
        const int wk_ = wk1-wk0+1;

        const int SK0 = V1.mra.rangeI(j0).firstIndex();
        const int SK1 = V0.mra.rangeI(j0).lastIndex();
        const int SKM = (SK1+SK0)/2;
        const int SK_ = SK1-SK0+1;

        const int WK0 = V1.rangeJ(j0).firstIndex();
        const int WK1 = V0.rangeJ(j0).lastIndex();
        const int WKM = (WK1+WK0)/2;
        const int WK_ = WK1-WK0+1;

        cerr << "U1.mra.rangeI(j0) = " << U1.mra.rangeI(j0) << endl;
        cerr << "U0.mra.rangeI(j0) = " << U0.mra.rangeI(j0) << endl;

        cerr << "V1.mra.rangeI(j0) = " << V1.mra.rangeI(j0) << endl;
        cerr << "V0.mra.rangeI(j0) = " << V0.mra.rangeI(j0) << endl;

        RealGeMatrix     A(SK_+WK_, sk_+wk_);
        RealGeMatrix     A1(SK_+WK_, sk_+wk_);
        RealGeMatrix     A2(SK_+WK_, sk_+wk_);
        RealGeMatrix     A3(SK_+WK_, sk_+wk_);
        RealGeMatrix     A4(SK_+WK_, sk_+wk_);

        // (Spline, Spline)
        for (int i=1; i<=SK_; ++i) {
            for (int j=1; j<=sk_; ++j) {
                const int k = sk0+j-1;
                const int K = SK0+i-1;

                cerr << "k = " << k << ", K = " << K << endl;

                double value;
                if (k<skm) {
                    if (K<SKM) {
                        value = -integral01(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 1)
                                +integral01(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 0);
                    } else {
                        value = -integral00(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 1)
                                +integral00(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 0);
                    }
                } else {
                    if (K<SKM) {
                        value = -integral11(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 1)
                                +integral11(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 0);
                    } else {
                        value = -integral10(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 1)
                                +integral10(j0, k, XBSpline, 0,
                                            j0, K, XBSpline, 0);
                    }
                }
                A1(i, j) = value;
            }
        }
        cerr << "1) A1 = " << A1 << endl;

        // (Wavelet, Spline)
        for (int i=1; i<=SK_; ++i) {
            for (int j=1; j<=wk_; ++j) {
                const int k = wk0+j-1;
                const int K = SK0+i-1;

                double value;
                if (k<wkm) {
                    if (K<SKM) {
                        value = -integral01(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 1)
                                +integral01(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 0);
                    } else {
                        value = -integral00(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 1)
                                +integral00(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 0);
                    }
                } else {
                    if (K<SKM) {
                        value = -integral11(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 1)
                                +integral11(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 0);
                    } else {
                        value = -integral10(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 1)
                                +integral10(j0, k, XWavelet, 0,
                                            j0, K, XBSpline, 0);
                    }
                }
                A2(i, sk_+j) = value;
            }
        }
        cerr << "2) A2 = " << A2 << endl;

        // (Spline, Wavelet)
        for (int i=1; i<=WK_; ++i) {
            for (int j=1; j<=sk_; ++j) {
                const int k = sk0+j-1;
                const int K = WK0+i-1;

                double value;
                if (k<skm) {
                    if (K<WKM) {
                        value = -integral01(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 1)
                                +integral01(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 0);
                    } else {
                        value = -integral00(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 1)
                                +integral00(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 0);
                    }
                } else {
                    if (K<WKM) {
                        value = -integral11(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 1)
                                +integral11(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 0);
                    } else {
                        value = -integral10(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 1)
                                +integral10(j0, k, XBSpline, 0,
                                            j0, K, XWavelet, 0);
                    }
                }
                A3(SK_+i, j) = value;
            }
        }
        cerr << "3) A3 = " << A3 << endl;

        // (Wavelet, Wavelet)
        for (int i=1; i<=WK_; ++i) {
            for (int j=1; j<=wk_; ++j) {
                const int k = wk0+j-1;
                const int K = WK0+i-1;

                double value;
                if (k<wkm) {
                    if (K<WKM) {
                        value = -integral01(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 1)
                                +integral01(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 0);
                    } else {
                        value = -integral00(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 1)
                                +integral00(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 0);
                    }
                } else {
                    if (K<WKM) {
                        value = -integral11(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 1)
                                +integral11(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 0);
                    } else {
                        value = -integral10(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 1)
                                +integral10(j0, k, XWavelet, 0,
                                            j0, K, XWavelet, 0);
                    }
                }
                A4(SK_+i, sk_+j) = value;
            }
        }
        cerr << "4) A4 = " << A4 << endl;

        A += A1;
        A += A2;
        A += A3;
        A += A4;

        cerr << "A = " << A << endl;



        Function<double>                gFunc(g);
        IntegralF<Gauss, PrimalBasis>   rhsIntegral0(gFunc, V0);
        IntegralF<Gauss, PrimalBasis>   rhsIntegral1(gFunc, V1);
        RealDenseVector                 b(SK_+WK_);

        for (int i=1; i<=SK_; ++i) {
            const int K = SK0+i-1;
            if (K>SKM) {
                b(i) = rhsIntegral0(j0, K, XBSpline, 0);
                     //+ u0*V0.mra.phi(0.0, j0, K, 0);
            } else {
                b(i) = rhsIntegral1(j0, K, XBSpline, 0);
                     //+ u0*V1.mra.phi(0.0, j0, K, 0);
            }
        }

        for (int i=1; i<=WK_; ++i) {
            const int K = WK0+i-1;
            if (K>WKM) {
                b(SK_+i) = rhsIntegral0(j0, K, XWavelet, 0);
                         //+ u0*V0.psi(0.0, j0, K, 0);
            } else {
                b(SK_+i) = rhsIntegral1(j0, K, XWavelet, 0);
                         //+ u0*V1.psi(0.0, j0, K, 0);
            }
        }

        cerr << "b = " << b << endl;

        RealGeMatrix    ADiff = A;
        RealDenseVector bDiff = b;


        ADiff -= _A;
        bDiff -= _b;

        cerr.precision(30);

        cerr << "ADiff = " << ADiff << ", amax(ADiff) = " << amax(ADiff) << endl;
        cerr << "bDiff = " << bDiff << ", amax(bDiff) = " << amax(bDiff) << endl;


        for (int i=bDiff.firstIndex(); i<=bDiff.lastIndex(); ++i) {
            cerr << "i = " << i << ", bDiff(i) = " << bDiff(i) << endl;
        }

        RealDenseVector  c;

        RealGeMatrix  AtA;
        blas::mm(Trans, NoTrans, 1.0, A, A, 0.0, AtA);


        RealGeMatrix     AtADiff = AtA;
        
        AtADiff -= _AtA;

        
        cerr << "AtADiff = " << AtADiff << endl;

        cerr << "amax(AtADiff) = " << amax(AtADiff) << endl;


        RealDenseVector Atb;
        blas::mv(Trans, 1.0, A, b, 0.0, Atb);
        Atb.engine().changeIndexBase(1);


        RealDenseVector  AtbDiff = Atb;
        
        AtbDiff -= _Atb;
        cerr << "AtbDiff = " << AtbDiff << ", amax(AtbDiff) = " << amax(AtbDiff) << endl;

        for (int i=AtbDiff.firstIndex(); i<=AtbDiff.lastIndex(); ++i) {
            cerr << "i = " << i << ", AtbDiff(i) = " << AtbDiff(i) << endl;
        }


        //cerr << "AtA = " << AtA << endl;
        //cerr << "Atb = " << Atb << endl;

        IntegerDenseVector      piv(AtA.numRows());

        flens::sv(AtA, piv, Atb);
        //cerr << "x = " << Atb << endl;
        c = Atb;

        /*
        const int mn = std::min(A.numRows(), A.numCols());
        IntegerDenseVector      piv(mn);

        flens::sv(A, piv, b);
        cerr << "x = " << b << endl;
        c = b;
        */


        c.engine().changeIndexBase(1);

        cerr << "c =  " << c << endl;
        cerr << "_c = " << _c << endl;

        c -= _c;
        cerr << "diff: c = " << c << endl;
}


