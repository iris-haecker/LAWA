#include <iostream>
#define SOLVER_DEBUG

#include <iris/iris.cxx>


using namespace lawa;
using namespace std;

double
g(double t)
{
    const int l  = 100;
    const int u0 = 0;

    return sin(2*M_PI*l*t) + 2*M_PI*l*cos(2*M_PI*l*t) + u0;
}

int
main()
{
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::DenseVector<Array<int> >                  IntDenseVector;
    typedef MyOperator<double>                               Operator;
    typedef MyPrecond<Operator>                              Precond;
    typedef MyRhs<double, Precond>                           Rhs;
    typedef MyGHS<Operator, Rhs, Precond>                    GHS;

    const int d  = 2;
    const int d_ = 4;

    const int    jMax = 3;
    const double eps = 0.00001;

    Operator              operatorA(d, d_, jMax);
    Precond               P(operatorA);
    Rhs                   rhs(Function<double>(g), operatorA, P);

    RealDenseVector       x(operatorA.numCols());

    GHS                   ghs(operatorA, rhs, P, 0.4, 0.012618, 2./7);

    double                nu;
    IndexSet<int>         Lambda;

    ghs.grow(x, rhs.norm, eps, nu, Lambda);

    cerr << "nu = " << nu << endl;
    cerr << "Lambda = " << Lambda << endl;
    
    
    SparseGeMatrix<CRS<double ,CRS_General> >  B;
    myRestrict(operatorA, P, Lambda, B);
    
    //cerr << "B = " << B << endl;
    
}