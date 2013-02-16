#include <iostream>
#define SOLVER_DEBUG

#include <iris/iris.cxx>
#include <iris/richardson.h>
#include <iris/richardson.tcc>


using namespace lawa;
using namespace std;

const int    l  = 1;
const double u0 = 0;

double
g(double t)
{
    return sin(2*M_PI*l*t) + 2*M_PI*l*cos(2*M_PI*l*t) + u0;
}


int
main()
{
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;

    typedef MyOperator<double>                               Operator;

    const int d  = 2;
    const int d_ = 4;

    const int jMax = 11;
    const double eps = 0.00001;
    const double tol = 0.00001;

    
    cerr << "Init:" << endl;
    
    Operator                         OperatorA(d, d_, jMax);
    MyApply<Operator>                A(OperatorA, eps, tol);
    MyApplyAtA<MyApply<Operator> >   AtA(A);

    RealDenseVector     Atb(OperatorA.numCols());
    RealDenseVector     x(OperatorA.numCols());

    RealDenseVector     b(OperatorA.numRows());


    cerr << "Init RHS:" << endl;

    MyRhsIntegral<double>  myRhsIntegral(Function<double>(g), OperatorA.V);

    for (int p=1; p<=OperatorA.V.getLastAbsoluteIndex(10); ++p) {
        b(p) = myRhsIntegral(p, 0);
    }

    cerr << "num rows = "<< A.numRows() << endl;

    for (int p=1; p<=x.length(); ++p) {
        x(p) = x.length() - p;
        // x(p) = p;
    }
    
    /*
    for (int p=1; p<=100; ++p) {
        cerr << "p = " << p << endl;
        b = OperatorA*x;
    }
    */

    
    /*
    cg(A, x, b, 0.001, 4000);
    */

    cerr << "Compute  A^T * b:" << endl;

    //Atb = transpose(A)*b;
    //cg(AtA, x, Atb, 0.001, 40);

    // richardson(A, b, 0.001, x, 0.001, 10);
    RealGeMatrix        MA1, MA2;
    OperatorA.densify(MA1, OperatorA.j1, true);
    OperatorA.densify(MA2, OperatorA.j1, false);

    /*
    cerr << "MA1 = " << MA1 << endl;
    cerr << "MA2 = " << MA2 << endl;
    */

    double norm = 0;
    for (int i=1; i<=MA2.numRows(); ++i) {
        for (int j=1; j<=MA2.numCols(); ++j) {
            double tmp = MA2(i,j) - MA1(i,j);
            norm += tmp*tmp;
        }
    }
    cerr << "norm = " << sqrt(norm) << endl;

/*
    cerr << "A.densify(jMax, MA):" << endl;

    A.densify(jMax, MA);

    RealDenseVector     vb(MA.numRows()), vx(MA.numCols());

    for (int p=1; p<=vb.length(); ++p) {
        vb(p) = myRhsIntegral(p, 0);
    }


    cerr << "MA = " << MA << endl;

    cerr << "richardson" << endl;

    //richardson(MA, vb, 0.0000000000001, vx, 0.00000001, 10);
    
    cerr << "cg" << endl;
    vx = 0;
    cg(A, vx, vb, 0.00000001, 40);
*/
}