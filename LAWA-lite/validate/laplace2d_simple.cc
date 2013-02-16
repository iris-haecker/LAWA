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


const int time_d  = 2;
const int time_d_ = 4;

const int space_d  = 2;
const int space_d_ = 4;

const int j0t = 3;
const int j0x = 3;

const int Jt = 1;
const int Jx = 1;


double
gx(double x)
{
//    return sin(2*M_PI*x);
    return 1;
}

double
gt(double t)
{
    return 1;
//    return (1+4*M_PI*M_PI*t);
}

void
spaceStiffnessMatrix(bool deriv, RealGeMatrix &A)
{
    const unsigned int d  = space_d;
    const unsigned int d_ = space_d_;

    PrimalBasis   U0(d, d_);
    PrimalBasis   U1(d, d_);

    PrimalBasis   V0(d, d_);
    PrimalBasis   V1(d, d_);

    U0.enforceBoundaryCondition<DirichletBC>();
    U1.enforceBoundaryCondition<DirichletBC>();
    V0.enforceBoundaryCondition<DirichletBC>();
    V1.enforceBoundaryCondition<DirichletBC>();

    const int _j0 = std::max(std::max(std::max(U0.j0, V0.j0), U1.j0), V1.j0);
    assert(j0x>=_j0);

    const int j1 = j0x + Jx;

    //cerr << "j0-1 = " << (j0-1) << endl;
    //cerr << "j1   = " << j1 << endl;

    Integral<Gauss,PrimalBasis,PrimalBasis>  integral00(U0, V0);
    Integral<Gauss,PrimalBasis,PrimalBasis>  integral01(U0, V1);
    Integral<Gauss,PrimalBasis,PrimalBasis>  integral10(U1, V0);
    Integral<Gauss,PrimalBasis,PrimalBasis>  integral11(U1, V1);


    IntegerDenseVector  ku0(_(j0x-1,j0x+Jx)), kv0(_(j0x-1,j0x+Jx));
    IntegerDenseVector  ku1(_(j0x-1,j0x+Jx)), kv1(_(j0x-1,j0x+Jx));
    IntegerDenseVector  kum(_(j0x-1,j0x+Jx)), kvm(_(j0x-1,j0x+Jx));

    IntegerDenseVector  KU(_(j0x-1,j0x+Jx)), KV(_(j0x-1,j0x+Jx));
    IntegerDenseVector  KU0(_(j0x-1,j0x+Jx)), KV0(_(j0x-1,j0x+Jx));
    IntegerDenseVector  KU1(_(j0x-1,j0x+Jx)), KV1(_(j0x-1,j0x+Jx));

    for (int _ju=j0x-1; _ju<=j0x+Jx; ++_ju) {
        ku0(_ju) = (_ju<j0x) ? U0.mra.rangeI(j0x).firstIndex()
                             : U0.rangeJ(_ju).firstIndex();
        ku1(_ju) = (_ju<j0x) ? U1.mra.rangeI(j0x).lastIndex()
                            : U1.rangeJ(_ju).lastIndex();
        kum(_ju) = (ku1(_ju)+ku0(_ju))/2;
        KU(_ju)  = ku1(_ju)-ku0(_ju)+1;
        KU0(_ju) = (_ju<j0x) ? 1
                             : KU0(_ju-1) + KU(_ju-1);
        KU1(_ju) = KU0(_ju) + KU(_ju) - 1;
    }

    /*
    cerr << "ku0 = " << ku0 << endl;
    cerr << "ku1 = " << ku1 << endl;
    cerr << "KU  = " << KU << endl;
    cerr << "KU0 = " << KU0 << endl;
    cerr << "KU1 = " << KU1 << endl;
    */

    for (int _jv=j0x-1; _jv<=j0x+Jx; ++_jv) {
        kv0(_jv) = (_jv<j0x) ? V1.mra.rangeI(j0x).firstIndex()
                             : V1.rangeJ(_jv).firstIndex();
        kv1(_jv) = (_jv<j0x) ? V0.mra.rangeI(j0x).lastIndex()
                             : V0.rangeJ(_jv).lastIndex();
        kvm(_jv) = (kv0(_jv)+kv1(_jv))/2;

        KV(_jv)  = kv1(_jv)-kv0(_jv)+1;
        KV0(_jv) = (_jv<j0x) ? 1
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

    A.engine().resize(KV1(j0x+Jx), KU1(j0x+Jx));

    for (int _jv=j0x-1; _jv<=j0x+Jx; ++_jv) {
        int    jv     = (_jv<j0x) ? j0x : _jv;
        XType  vType  = (_jv<j0x) ? XBSpline : XWavelet;

        for (int _ju=j0x-1; _ju<=j0x+Jx; ++_ju) {
            int    ju     = (_ju<j0x) ? j0x : _ju;
            XType  uType  = (_ju<j0x) ? XBSpline : XWavelet;

            // (uType, vType)
            for (int i=0; i<KV(_jv); ++i) {
                for (int j=0; j<KU(_ju); ++j) {
                    const int ku = ku0(_ju) + j;
                    const int kv = kv0(_jv) + i;

                    double value;

                    if (ku<kum(_ju)) {
                        if (kv<kvm(_jv)) {
                            if (deriv) {
                                value = integral01(ju, ku, uType, 1,
                                                   jv, kv, vType, 1);
                            } else {
                                value = integral01(ju, ku, uType, 0,
                                                   jv, kv, vType, 0);
                            }
                        } else {
                            if (deriv) {
                                value = integral00(ju, ku, uType, 1,
                                                   jv, kv, vType, 1);
                            } else {
                                value = integral00(ju, ku, uType, 0,
                                                   jv, kv, vType, 0);
                            }
                        }
                    } else {
                        if (kv<kvm(_jv)) {
                            if (deriv) {
                                value = integral11(ju, ku, uType, 1,
                                                   jv, kv, vType, 1);
                            } else {
                                value = integral11(ju, ku, uType, 0,
                                                   jv, kv, vType, 0);
                            }
                        } else {
                            if (deriv) {
                                value = integral10(ju, ku, uType, 1,
                                                   jv, kv, vType, 1);
                            } else {
                                value = integral10(ju, ku, uType, 0,
                                                   jv, kv, vType, 0);
                            }
                        }
                    }
                    A(KV0(_jv)+i, KU0(_ju)+j) = value;
                }
            }
        }
    }
//  cerr << "A = " << A << endl;
}

void
setupRhs(RealDenseVector &b1, RealDenseVector &b2)
{
    PrimalBasis   V0x(space_d, space_d_);
    PrimalBasis   V1x(space_d, space_d_);

    V0x.enforceBoundaryCondition<DirichletBC>();
    V1x.enforceBoundaryCondition<DirichletBC>();

    const int _j0x = std::max(V0x.j0, V1x.j0);
    assert(j0x>=_j0x);

    IntegerDenseVector  kv0x(_(j0x-1,j0x+Jx));
    IntegerDenseVector  kv1x(_(j0x-1,j0x+Jx));
    IntegerDenseVector  kvmx(_(j0x-1,j0x+Jx));

    IntegerDenseVector  KVx(_(j0x-1,j0x+Jx));
    IntegerDenseVector  KV0x(_(j0x-1,j0x+Jx));
    IntegerDenseVector  KV1x(_(j0x-1,j0x+Jx));

    for (int _jvx=j0x-1; _jvx<=j0x+Jx; ++_jvx) {
        kv0x(_jvx) = (_jvx<j0x) ? V1x.mra.rangeI(j0x).firstIndex()
                                : V1x.rangeJ(_jvx).firstIndex();
        kv1x(_jvx) = (_jvx<j0x) ? V0x.mra.rangeI(j0x).lastIndex()
                                : V0x.rangeJ(_jvx).lastIndex();
        kvmx(_jvx) = (kv0x(_jvx)+kv1x(_jvx))/2;

        KVx(_jvx)  = kv1x(_jvx)-kv0x(_jvx)+1;
        KV0x(_jvx) = (_jvx<j0x) ? 1
                                : KV0x(_jvx-1) + KVx(_jvx-1);
        KV1x(_jvx) = KV0x(_jvx) + KVx(_jvx) - 1;
    }

    /*
    cerr << "kv0x = " << kv0x << endl;
    cerr << "kv1x = " << kv1x << endl;
    cerr << "kvmx = " << kvmx << endl;
    cerr << "KVx  = " << KVx << endl;
    cerr << "KV0x = " << KV0x << endl;
    cerr << "KV1x = " << KV1x << endl;
    */

    PrimalBasis   V0t(time_d, time_d_);
    PrimalBasis   V1t(time_d, time_d_);

    V0t.enforceBoundaryCondition<DirichletBC>();
    V1t.enforceBoundaryCondition<DirichletBC>();

    const int _j0t = std::max(V0t.j0, V1t.j0);
    assert(j0t>=_j0t);


    IntegerDenseVector  kv0t(_(j0t-1,j0t+Jt));
    IntegerDenseVector  kv1t(_(j0t-1,j0t+Jt));
    IntegerDenseVector  kvmt(_(j0t-1,j0t+Jt));

    IntegerDenseVector  KVt(_(j0t-1,j0t+Jt));
    IntegerDenseVector  KV0t(_(j0t-1,j0t+Jt));
    IntegerDenseVector  KV1t(_(j0t-1,j0t+Jt));

    for (int _jvt=j0t-1; _jvt<=j0t+Jt; ++_jvt) {
        kv0t(_jvt) = (_jvt<j0t) ? V1t.mra.rangeI(j0t).firstIndex()
                                : V1t.rangeJ(_jvt).firstIndex();
        kv1t(_jvt) = (_jvt<j0t) ? V0t.mra.rangeI(j0t).lastIndex()
                                : V0t.rangeJ(_jvt).lastIndex();
        kvmt(_jvt) = (kv0t(_jvt)+kv1t(_jvt))/2;

        KVt(_jvt)  = kv1t(_jvt)-kv0t(_jvt)+1;
        KV0t(_jvt) = (_jvt<j0t) ? 1
                                : KV0t(_jvt-1) + KVt(_jvt-1);
        KV1t(_jvt) = KV0t(_jvt) + KVt(_jvt) - 1;
    }

    /*
    cerr << "kv0t = " << kv0t << endl;
    cerr << "kv1t = " << kv1t << endl;
    cerr << "kvmt = " << kvmt << endl;
    cerr << "KVt  = " << KVt << endl;
    cerr << "KV0t = " << KV0t << endl;
    cerr << "KV1t = " << KV1t << endl;
    */

    b1.engine().resize(KV1x(j0x+Jx));
    b2.engine().resize(KV1t(j0t+Jt));
    
    Function<double>                gxFunc(gx);
    Function<double>                gtFunc(gt);

    IntegralF<Gauss, PrimalBasis>   rhsIntegralX0(gxFunc, V0x);
    IntegralF<Gauss, PrimalBasis>   rhsIntegralX1(gxFunc, V1x);

    IntegralF<Gauss, PrimalBasis>   rhsIntegralT0(gtFunc, V0t);
    IntegralF<Gauss, PrimalBasis>   rhsIntegralT1(gtFunc, V1t);

    for (int _jvx=j0x-1; _jvx<=j0x+Jx; ++_jvx) {
        int    jvx     = (_jvx<j0x) ? j0x : _jvx;
        XType  vxType  = (_jvx<j0x) ? XBSpline : XWavelet;

        for (int ix=0; ix<KVx(_jvx); ++ix) {
            const int kvx = kv0x(_jvx) + ix;

            double valueX = 0;
            if (kvx<=kvmx(_jvx)) {
                valueX = rhsIntegralX0(jvx, kvx, vxType, 0);
            } else {
                valueX = rhsIntegralX1(jvx, kvx, vxType, 0);
            }
            b1(KV0x(_jvx)+ix) = valueX;
        }
    }

    for (int _jvt=j0t-1; _jvt<=j0t+Jt; ++_jvt) {
        int    jvt     = (_jvt<j0t) ? j0t : _jvt;
        XType  vtType  = (_jvt<j0t) ? XBSpline : XWavelet;

        for (int it=0; it<KVt(_jvt); ++it) {
            const int kvt = kv0t(_jvt) + it;

            double valueT = 0;
            if (kvt<=kvmt(_jvt)) {
                valueT = rhsIntegralT1(jvt, kvt, vtType, 0);
            } else {
                valueT = rhsIntegralT0(jvt, kvt, vtType, 0);
            }
            b2(KV0t(_jvt)+it) = valueT;
        }
    }

}
/*
void
setupStiffnessMatrix(const RealGeMatrix &A1,
                     const RealGeMatrix &A1x,
                     const RealGeMatrix &A2,
                     const RealGeMatrix &A2t,
                     RealGeMatrix &A)
{
    const int m = A1.numRows();
    const int n = A1.numCols();
    const int M = A2.numRows();
    const int N = A2.numCols();

    assert(A1x.numRows()==m);
    assert(A1x.numCols()==n);

    assert(A2t.numRows()==M);
    assert(A2t.numCols()==N);

    RealGeMatrix B(A2.numRows(), A2.numCols());

    A.engine().resize(m*M, n*N, 1, 1);
    for (int i=0; i<m; ++i) {
        Range<int> rowBlock = _(1+i*M, (i+1)*M);

        for (int j=0; j<n; ++j) {
            Range<int> colBlock = _(1+j*N, (j+1)*N);

            A(rowBlock, colBlock) = A2t;
            A(rowBlock, colBlock) *= -A1(i+1,j+1);

            B = A2;
            B *= A1x(i+1,j+1);

            A(rowBlock, colBlock) += B;
        }
    }
}
*/

void
setupStiffnessMatrix(const RealGeMatrix &A1,
                     const RealGeMatrix &A1x,
                     RealGeMatrix &A)
{
    const int m = A1.numRows();
    const int n = A1.numCols();
    const int M = A1x.numRows();
    const int N = A1x.numCols();


    RealGeMatrix B(M, N);

    A.engine().resize(m*M, n*N, 1, 1);

    for (int i=0; i<m; ++i) {
        Range<int> rowBlock = _(1+i*M, (i+1)*M);

        for (int j=0; j<n; ++j) {
            Range<int> colBlock = _(1+j*N, (j+1)*N);

            A(rowBlock, colBlock) = A1;
            A(rowBlock, colBlock) *= A1x(i+1,j+1);

            B = A1x;
            B *= A1(i+1,j+1);
            A(rowBlock, colBlock) += B;
        }
    }
}

void
assembleRhs(const RealDenseVector &b1,
            const RealDenseVector &b2,
            RealDenseVector &b)
{
    const int n1 = b1.length();
    const int n2 = b2.length();
    
    b.engine().resize(n1*n2);

    for (int i=0; i<n1; ++i) {
        Range<int> range = _(1+i*n2, (i+1)*n2);
        b(range) = b2;
        b(range) *= b1(i+1);
        //cout << "range = " << range << endl;
    }
}


double
eval(const RealDenseVector &coeff, double x, double t)
{
    PrimalBasis   U0x(space_d, space_d_);
    PrimalBasis   U1x(space_d, space_d_);

    U0x.enforceBoundaryCondition<DirichletBC>();
    U1x.enforceBoundaryCondition<DirichletBC>();

    const int _j0x = std::max(U0x.j0, U1x.j0);
    assert(j0x>=_j0x);

    IntegerDenseVector  ku0x(_(j0x-1,j0x+Jx));
    IntegerDenseVector  ku1x(_(j0x-1,j0x+Jx));
    IntegerDenseVector  kumx(_(j0x-1,j0x+Jx));

    IntegerDenseVector  KUx(_(j0x-1,j0x+Jx));
    IntegerDenseVector  KU0x(_(j0x-1,j0x+Jx));
    IntegerDenseVector  KU1x(_(j0x-1,j0x+Jx));

    for (int _jux=j0x-1; _jux<=j0x+Jx; ++_jux) {
        ku0x(_jux) = (_jux<j0x) ? U0x.mra.rangeI(j0x).firstIndex()
                                : U0x.rangeJ(_jux).firstIndex();
        ku1x(_jux) = (_jux<j0x) ? U1x.mra.rangeI(j0x).lastIndex()
                                : U1x.rangeJ(_jux).lastIndex();
        kumx(_jux) = (ku1x(_jux)+ku0x(_jux))/2;
        KUx(_jux)  = ku1x(_jux)-ku0x(_jux)+1;
        KU0x(_jux) = (_jux<j0x) ? 1
                                : KU0x(_jux-1) + KUx(_jux-1);
        KU1x(_jux) = KU0x(_jux) + KUx(_jux) - 1;
    }

    /*
    cerr << "ku0x = " << ku0x << endl;
    cerr << "ku1x = " << ku1x << endl;
    cerr << "kumx = " << kumx << endl;
    cerr << "KUx  = " << KUx << endl;
    cerr << "KU0x = " << KU0x << endl;
    cerr << "KU1x = " << KU1x << endl;
    */

    PrimalBasis   U0t(time_d, time_d_);
    PrimalBasis   U1t(time_d, time_d_);

    U0t.enforceBoundaryCondition<DirichletBC>();
    U1t.enforceBoundaryCondition<DirichletBC>();

    const int _j0t = std::max(U0t.j0, U1t.j0);
    assert(j0t>=_j0t);

    IntegerDenseVector  ku0t(_(j0t-1,j0t+Jt));
    IntegerDenseVector  ku1t(_(j0t-1,j0t+Jt));
    IntegerDenseVector  kumt(_(j0t-1,j0t+Jt));

    IntegerDenseVector  KUt(_(j0t-1,j0t+Jt));
    IntegerDenseVector  KU0t(_(j0t-1,j0t+Jt));
    IntegerDenseVector  KU1t(_(j0t-1,j0t+Jt));

    for (int _jut=j0t-1; _jut<=j0t+Jt; ++_jut) {
        ku0t(_jut) = (_jut<j0t) ? U0t.mra.rangeI(j0t).firstIndex()
                                : U0t.rangeJ(_jut).firstIndex();
        ku1t(_jut) = (_jut<j0t) ? U1t.mra.rangeI(j0t).lastIndex()
                                : U1t.rangeJ(_jut).lastIndex();
        kumt(_jut) = (ku1t(_jut)+ku0t(_jut))/2;
        KUt(_jut)  = ku1t(_jut)-ku0t(_jut)+1;
        KU0t(_jut) = (_jut<j0t) ? 1
                                : KU0t(_jut-1) + KUt(_jut-1);
        KU1t(_jut) = KU0t(_jut) + KUt(_jut) - 1;
    }

    /*
    cerr << "ku0t = " << ku0t << endl;
    cerr << "ku1t = " << ku1t << endl;
    cerr << "kumt = " << kumt << endl;
    cerr << "KUt  = " << KUt << endl;
    cerr << "KU0t = " << KU0t << endl;
    cerr << "KU1t = " << KU1t << endl;
    */

    double value = 0;
    int col = 1;

    for (int _jux=j0x-1; _jux<=j0x+Jx; ++_jux) {
        int jux = (_jux<j0x) ? j0x : _jux;

        for (int ix=0; ix<KUx(_jux); ++ix) {
            const int kux = ku0x(_jux) + ix;

            double valueX = 0;
            if (kux<=kumx(_jux)) {
                if (_jux<j0x) {
//
//                  Scaling function
//
                    valueX = U0x.mra.phi(x, jux, kux, 0);
                } else {
//
//                  Wavelet
//
                    //valueX = U0x.psi(x, jux, kux, 0);
                }
            } else {
                if (_jux<j0x) {
//
//                  Scaling function
//
                    valueX = U1x.mra.phi(x, jux, kux, 0);
                } else {
//
//                  Wavelet
//
                    //valueX = U1x.psi(x, jux, kux, 0);
                }
            }

            /*
            cerr << "jux = " << jux
                 << ", kux = " << kux
                 << ", valueX = " << valueX
                 << ", typeX = " << int(_jux<j0x)
                 << endl;
            */

            for (int _jut=j0t-1; _jut<=j0t+Jt; ++_jut) {
                int jut = (_jut<j0t) ? j0t : _jut;

                for (int it=0; it<KUt(_jut); ++it, ++col) {
                    const int kut = ku0t(_jut) + it;

                    double valueT = 0;
                    if (kut<=kumt(_jut)) {
                        if (_jut<j0t) {
//
//                          Scaling function
//
                            valueT = U0t.mra.phi(t, jut, kut, 0);
                        } else {
//
//                          Wavelet
//
                            //valueT = U0t.psi(t, jut, kut, 0);
                        }
                    } else {
                        if (_jut<j0t) {
//
//                          Scaling function
//
                            valueT = U1t.mra.phi(t, jut, kut, 0);
                        } else {
//
//                          Wavelet
//
                            //valueT = U1t.psi(x, jut, kut, 0);
                        }
                    }
                    value += coeff(col)*valueX*valueT;

                    /*
                    cerr << "jut = " << jut
                              << ", kut = " << kut
                              << ", valueT = " << valueT
                              << ", typeT = " << int(_jut<j0t)
                              << ", col = " << col
                              << endl;
                    */
                }
            }


        }
    }

    return value;
}

int
main()
{
    RealGeMatrix     A1, A1x, A2, A2t, A;
    RealDenseVector  b1, b2, b;


    spaceStiffnessMatrix(false, A1);
    spaceStiffnessMatrix(true, A1x);

    //timeStiffnessMatrix(false, A2);
    //timeStiffnessMatrix(true, A2t);

    cerr << "A1 = " << A1 << endl;
    cerr << "A1x = " << A1x << endl;

/*
    cerr << "A2 = " << A2 << endl;
    cerr << "A2t = " << A2t << endl;
*/
    setupRhs(b1, b2);
    assembleRhs(b1, b2, b);

    cerr << "b1 = " << b1 << endl;
    cerr << "b2 = " << b2 << endl;

    setupStiffnessMatrix(A1, A1x, A);
    
    
    IntegerDenseVector   piv(A.numRows());
    flens::sv(A, piv, b);

    cerr << "x = " << b << endl;

/*
    setupStiffnessMatrix(A1, A1x, A2, A2t, A);

//
//  Solve the least square problem
//
    RealGeMatrix  AtA;
    blas::mm(Trans, NoTrans, 1.0, A, A, 0.0, AtA);

    RealDenseVector Atb;
    blas::mv(Trans, 1.0, A, b, 0.0, Atb);

    AtA.engine().changeIndexBase(1,1);
    Atb.engine().changeIndexBase(1);

    IntegerDenseVector   piv(AtA.numRows());


    flens::sv(AtA, piv, Atb);
    cerr << "Atb = " << Atb << endl;
*/


    const int numPointsT = 100;
    const int numPointsX = 100;

    for (int pt=0; pt<=numPointsT; ++pt) {
        const double t = double(pt)/numPointsT;

        for (int px=0; px<=numPointsX; ++px) {
            const double x = double(px)/numPointsX;

            cout << x << " " << t << " " << eval(b, x, t) << endl;
        }
        cout << endl;
    }
}
