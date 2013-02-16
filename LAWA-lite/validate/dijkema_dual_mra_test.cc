#include <fstream>
#include <iostream>

#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

int
main()
{
    const int d  = 3;
    const int d_ = 3;

    MRA<double,Dual,Interval,Dijkema>  mra_(d, d_);
    //mra_.enforceBoundaryCondition<DirichletBC>();

    const int m = mra_.M0_.numRows();
    const int n = mra_.M0_.numCols();
    cerr << m << " x " << n << endl;
    GeMatrix<FullStorage<double, ColMajor> > D(m,n);
    DenseVector<Array<double> > e(n);
    for (int i=1; i<=n; ++i) {
        e(i) = 1.;
        D(_,i) = mra_.M0_*e;
        e(i) = 0.;
    }
    cout << D << endl;


    return 0;
}
