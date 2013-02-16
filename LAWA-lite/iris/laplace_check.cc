#include <iostream>
//#define SOLVER_DEBUG

#include <iris/iris.cxx>

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
    typedef DenseVector<Array<double> >                RealDenseVector;
    typedef GeMatrix<FullStorage<double, ColMajor> >   RealGeMatrix;
    typedef MyOperator<double>                         Operator;

    const int d  = 2;
    const int d_ = 4;

    const int jMax = 10;
    const double eps = 0.0001;
    const double tol = 0.00000000000000000001;

    
    cerr << "Init:" << endl;
    
    Operator            OperatorA(d, d_, jMax);
    MyApply<Operator>   A(OperatorA, eps, tol);

    RealDenseVector     x(A.numCols());
    RealDenseVector     b(A.numRows()), b2;

    cerr << "OperatorA.numRows() = " << OperatorA.numRows() << std::endl;
    cerr << "OperatorA.numCols() = " << OperatorA.numCols() << std::endl;

    cerr << "A.numRows() = " << A.numRows() << std::endl;
    cerr << "A.numCols() = " << A.numCols() << std::endl;

    cerr << "x.length() = " << x.length() << std::endl;
    cerr << "b.length() = " << b.length() << std::endl;


    cerr << "Init RHS:" << endl;

    MyRhsIntegral<double>  myRhsIntegral(Function<double>(g), OperatorA.V);

    for (int p=1; p<=OperatorA.V.getLastAbsoluteIndex(min(8,jMax)); ++p) {
        b(p) = myRhsIntegral(p, 0);
    }

    //cg(A, x, b, 0.001, 4000);

    cerr << "Densify Operator A:" << endl;
    RealGeMatrix   MA;
    A.densify(jMax, MA);

    //std::cerr << "MA = " << MA << std::endl;

    RealDenseVector p;
    A.precond(jMax, p);
    DiagonalMatrix<double>   P(p);

    cerr << "pcg:" << endl;

    //std::cerr << "P = " << P << std::endl;

    //RealGeMatrix   PA(P.numRows(), A.numCols());

    // PA = P*A;
    //blas::mm(cxxblas::NoTrans, cxxblas::NoTrans, 1., P, MA, 0., PA);

    //std::cerr << "PA = " << PA << std::endl;



    x = 0;
    pcg(P, MA, x, b, 0.00001, 40000);

    x = 0;
    //cg(MA, x, b, 0.001, 40000);


/*
    cerr << "Densify Operator A:" << endl;
    RealGeMatrix   MA;
    A.densify(jMax, MA);

    cerr << "Test:" << endl;
    x(1)  = 1;
    x(2)  = 0.5;
    x(3)  = 0.25;
    x(4)  = 0.125;
    x(5)  = 0.0064;
    x(6)  = 0.0032;
    x(7)  = 0.0016;
    x(8)  = 0.00016;
    x(9)  = 0.0008;
    x(10) = 0.00008;
    b2  = A*x;
    b   = A*x;
    b  -= MA*x;
    
    double bNormSquare = b*b;
    
    cerr << "b = " << b << endl;
    cerr << "b2 = " << b2 << endl;
    cerr << "bNormSquare = " << bNormSquare << endl;
    cerr << "bNorm = " << sqrt(bNormSquare) << endl;
    */
    
    
    //cerr << "MA = " << MA << endl;
    /*
    cg(MA, x, b, 0.001, 40000);
    */
}