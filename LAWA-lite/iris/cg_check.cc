#include <iostream>
#define SOLVER_DEBUG

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
    typedef GeMatrix<FullStorage<double, ColMajor> >    RealGeMatrix;
    typedef DenseVector<Array<double> >                 RealDenseVector;

    const int n = 4000;

    RealGeMatrix      A(n,n);
    RealDenseVector   x(n), b(n);
    
    A.diag( 0) =  2;
    A.diag( 1) = -1;
    A.diag(-1) = -1;

    b = 1;

    cg(A, x, b, 0.001, 40000);
    
    //cout << "x = " << x << endl;
}