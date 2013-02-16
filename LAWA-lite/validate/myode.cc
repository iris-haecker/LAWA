#include <cmath>
#include <fstream>
#include <iostream>

#include <lawa/lawa.h>

#include <validate/mybasis.h>
#include <validate/mybasis.tcc>

#include <validate/myintegral.h>
#include <validate/myintegral.tcc>

#include <validate/myrhsintegral.h>
#include <validate/myrhsintegral.tcc>


using namespace lawa;
using namespace std;

typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
typedef flens::DenseVector<Array<double> >               RealDenseVector;
typedef flens::DenseVector<Array<int> >                  IntegerDenseVector;

const int    k  = 1;
const double u0 = 1;

double
g(double t)
{
    return u0;
//    return sin(2*M_PI*k*t) + 2*M_PI*k*cos(2*M_PI*k*t) + u0;
}

int
main()
{
    MyBasis<double> U(3, 5);
    MyBasis<double> V(3, 5);

    V.basisRight.enforceBoundaryCondition<DirichletBC>();

    const int j0 = min(U.j0(), V.j0());
    const int J = j0+2;
    const int jMax = J+3;

    MyIntegral<double>  myIntegral(U, V);
        
    Range<int>    KU(U.getFirstAbsoluteIndex(j0-1), U.getLastAbsoluteIndex(J));
    Range<int>    KV(V.getFirstAbsoluteIndex(j0-1), V.getLastAbsoluteIndex(J));

//
//  A'*A berechnen
//

    RealGeMatrix     AtA(KU, KU);

    const double tol = 0.0001;

    for (int row=KU.firstIndex(); row<=KU.lastIndex(); ++row) {
        
        std::cerr << "row = " << row << std::endl;
        
        for (int col=KU.firstIndex(); col<=KU.lastIndex(); ++col) {

            AtA(row,col) = 0;

            for (int j=j0-1; j<=jMax; ++j) {
                double value = 0;

                int k0 = V.getFirstAbsoluteIndex(j);
                int k1 = V.getLastAbsoluteIndex(j);

                for (int k=k0; k<=k1; ++k) {

                    value += (-myIntegral(row,0,k,1) + myIntegral(row,0,k,0))
                           * (-myIntegral(col,0,k,1) + myIntegral(col,0,k,0));
                }

                AtA(row,col) += value;

                /*
                if ((j>=J) && (value>=tol*(k1-k0+1))) {
                    std::cerr << "j = " << j << std::endl;
                    std::cerr << "value = " << value << std::endl;
                    std::cerr << "k1-k0+1 = " << k1-k0+1 << std::endl;
                    std::cerr << "tol*(k1-k0+1) = " << tol*(k1-k0+1)
                              << std::endl
                              << std::endl;
                }

                if ((j>jMax) || (j>=J) && (value<tol*(k1-k0+1))) {
                    cerr << "row = " << row
                         << ", col = " << col
                         << ", j = " << j << endl;
                    break;
                }
                */

            }
        }
    }

    cerr << "AtA = " << AtA << endl;

//
//  A berechnen
//
    RealGeMatrix     A(KV, KU);

    for (int row=KV.firstIndex(); row<=KV.lastIndex(); ++row) {
        for (int col=KU.firstIndex(); col<=KU.lastIndex(); ++col) {

            A(row,col) = -myIntegral(col,0,row,1) + myIntegral(col,0,row,0);

        }
    }

    
//    blas::mm(Trans, NoTrans, 1.0, A, A, 0.0, AtA);

//    cerr << "AtA = " << AtA << endl;
    cerr << "A = " << A << endl;


//
//  b berechnen
//
    RealDenseVector  b(KV);

    MyRhsIntegral<double>  myRhsIntegral(Function<double>(g), V);

    for (int kv=KV.firstIndex(); kv<=KV.lastIndex(); ++kv) {
        b(kv) = myRhsIntegral(kv, 0);
    }
    cerr << "b = " << b << endl;

//
//  A'*b berechnen
//
    RealDenseVector  Atb(KU);

    blas::mv(Trans, 1.0, A, b, 0.0, Atb);
    cerr << "Atb = " << Atb << endl;

//
//  Solve least square problem Ax=b (by solving A'*A x = A'*b)
//

    IntegerDenseVector      piv(AtA.numRows());
    Atb.engine().changeIndexBase(1);

    flens::sv(AtA, piv, Atb);
    cerr << "AtA = " << AtA << endl;
    cerr << "piv = " << piv << endl;
    cerr << "x = " << Atb << endl;

//
//  Evaluate solution
//
    const int numPoints = 5000;
    for (int p=0; p<=numPoints; ++p) {
        const double x = double(p)/numPoints;

        double y = 0;
        for (int absIndex=1; absIndex<=Atb.length(); ++absIndex) {
            y += Atb(absIndex) * U(x, absIndex, 0);
        }
        cout << x << " " << y << endl;
    }
}
