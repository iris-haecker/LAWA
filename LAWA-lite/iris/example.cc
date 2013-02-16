#include <iostream>
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
    typedef flens::DenseVector<Array<double> >      RealDenseVector;
    typedef MyOperator<double>                      Operator;

    const int d = 4;
    const int d_ = 6;

    const int jMax = 16;
    const double eps = 0.00001;
    const double tol = 0.001;

    
    cerr << "Init:" << endl;
    
    Operator            _A(d, d_, jMax);
    MyApply<Operator>   A(_A, eps, tol);
    MyApply_AtA<MyApply<Operator> > AtA(A);

    cerr << "[" << _A.firstRow() << "," << _A.lastRow() << "]x";
    cerr << "[" << _A.firstCol() << "," << _A.lastCol() << "];" << endl;

    RealDenseVector     x(_A.numCols(), _A.firstCol()),
                        Atb(_A.numCols(), _A.firstCol());
    RealDenseVector     b(_A.numRows(), _A.firstRow());


    cerr << "Init RHS:" << endl;

    MyRhsIntegral<double>  myRhsIntegral(Function<double>(g), _A.V);

    for (int k=b.firstIndex(); k<=b.lastIndex(); ++k) {
        b(k) = myRhsIntegral(k, 0);
    }
    Atb = transpose(A)*b;
    

    cerr << "Solve:  A^T A x = A^T b" << endl;

    cg(AtA, x, Atb, 0.001, 40);

}