#include <iostream>
#define SOLVER_DEBUG

#include <iris/iris.cxx>


using namespace lawa;
using namespace std;

double
g(double t)
{
    const int l  = 1;
    const int u0 = 0;

    // return sin(2*M_PI*l*t) + 2*M_PI*l*cos(2*M_PI*l*t) + u0;
    // return t;
    return 2;
}

int
main()
{
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::DenseVector<Array<int> >                  IntDenseVector;
    typedef MyOperator<double>                               Operator;
    typedef MyRhs<double>                                    Rhs;

    const int d  = 3;
    const int d_ = 5;

    const int jMax = 3;
    const double eps = 0.00001;


    Operator             OperatorA(d, d_, jMax);
    Rhs                  rhs(Function<double>(g), OperatorA);

    MyApply<Operator>    A(OperatorA, eps);

    RealDenseVector      x(A.numCols());
    RealDenseVector      b;

    rhs.filter(0.2, b);
    cerr << "b = " << b << endl;

    rhs.filter(0.1, b);
    cerr << "b = " << b << endl;

    rhs.filter(0, b);
    cerr << "b = " << b << endl;
}