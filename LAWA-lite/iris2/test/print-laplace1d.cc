#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

int
main()
{
    Laplace1D<double>  A(2,4);

    int j1 = A.j0;
    int j2 = A.j0;

    int minK1 = A.V.rangeI(j1).firstIndex();
    int maxK1 = A.V.rangeI(j1).lastIndex();

    for (int k1=minK1; k1<=maxK1; ++k1) {

        Index1D nu(j1, k1, XBSpline);

        int minK2 = A.minK2(j1, k1, XBSpline, j2, XBSpline);
        int maxK2 = A.maxK2(j1, k1, XBSpline, j2, XBSpline);

        for (int k2=minK2; k2<=maxK2; ++k2) {
            Index1D lambda(j2, k2, XBSpline);

            cout << A(nu, lambda) << "  ";
        }
        cout << endl;
    }
}