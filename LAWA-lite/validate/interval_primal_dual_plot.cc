#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
//#include <extensions/extensions.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    
    //const int          J     = 1;
    const unsigned int d     = 2;
    const unsigned int d_    = 4;
    //const unsigned int deriv = 0;

    typedef Basis<double,Primal,Interval,Dijkema>  PrimalBasis;
    typedef Basis<double,Dual,Interval,Dijkema>    DualBasis;

    PrimalBasis primal(d, d_);
    DualBasis   dual(d, d_);
    
    //primal.enforceBoundaryCondition<DirichletBC>();
    
    const int j0 = std::min(primal.j0, dual.j0);
    
    cerr << "j0 = " << j0 << endl;
    cerr << "primal.M1.numCols() = " << primal.M1.numCols() << endl;
    cerr << "primal.mra.rangeI(j0) = " << primal.rangeJ(j0) << endl;

    const int k0 = primal.rangeJ(j0).firstIndex();
    const int k1 = primal.rangeJ(j0).lastIndex();

    const int N = 1000;
    for (int k=k0; k<=k1; ++k) {
        cerr << "k = " << k << endl;
        for (int i=0; i<=N; ++i) {
            const double x = i/double(N);
            //cout << x << " " << primal.mra.phi(x, j0, k, 0) << endl;
            cout << x << " " << primal.psi(x, j0, k, 0) << endl;
        }
        cout << endl << endl;
    }

    return 0;
}
