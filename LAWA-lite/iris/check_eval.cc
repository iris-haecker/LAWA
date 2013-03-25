#include <iostream>
#define SOLVER_DEBUG


#include <iris/iris.cxx>

using namespace flens;
using namespace lawa;
using namespace std;

int
main()
{
    typedef flens::DenseVector<Array<double> >              RealDenseVector;
    typedef MyOperator<double>                              Operator;
    typedef MyEval<double>                                  Eval;

    const int d  = 2;
    const int d_ = 4;

    const int    jMax = 8;
    const double eps = 0.0000001;

    Operator              A(d, d_, jMax);
    RealDenseVector       u(A.numCols());

    u(1) = 1;
    u(2) = 0.7;

    const int N = 1000;

    Eval sol(A.U, u);

    sol.dump(N, "test.dat");

}