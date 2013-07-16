#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

int
main()
{
    Laplace1D<double>  A(2,4);

    int j1 = A.j0;

    int minK1 = A.V.rangeI(j1).firstIndex();
    int maxK1 = A.V.rangeI(j1).lastIndex();

    for (int k1=minK1; k1<=maxK1; ++k1) {
        typedef typename IndexSet<Index1D>::iterator const_set_it;

        Index1D            nu(j1, k1, XBSpline);
        IndexSet<Index1D>  lambdas = lambdaTilde1d(nu, A, 3, j1, j1+3);

        for (const_set_it it=lambdas.begin(); it!=lambdas.end(); ++it) {
            cerr << "*it = " << *it << "  " << A(nu, *it) << "  " << endl;
            
            // cout << A(nu, *it) << "  ";
        }
        cout << endl;
    }
}