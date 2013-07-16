#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

int
main()
{
    Laplace1D<double>                   A(4,6);
    PreconditionerLaplace1D<double>     P(A);

    const int j0 = A.j0;
    const int j1 = j0 + 1;

    const int k1 = A.V.rangeI(j0).firstIndex();
    const int k2 = A.V.rangeI(j0).lastIndex();

    for (int k=k1; k<=k2; ++k) {
        Index1D  lambda(j0, k, XBSpline);

        cout << lambda << ": " << P(lambda) << endl;
    }

    for (int j=j0; j<=j1; ++j) {
        const int k1 = A.V.rangeJ(j).firstIndex();
        const int k2 = A.V.rangeJ(j).lastIndex();

        for (int k=k1; k<=k2; ++k) {
            Index1D  lambda(j, k, XWavelet);

            cout << lambda << ": " << P(lambda) << endl;
        }
    }

}