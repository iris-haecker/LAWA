#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

int
main()
{
    typedef double              T;
    typedef SolLaplace1D<T>     RefSol;

    const int           d = 2;
    const int           d_ = 4;

    Function<T>         rhs_func(RefSol::rhs);
    CompoundBasis<T>    V(d, d_);

    V.enforceBoundaryCondition<NoBC, NoBC>();

    RHSIntegral1D<T>    rhsIntegral(V, rhs_func, 120);

    RefSol::setExample(1, T(1));

    const int j0 = V.j0();
    const int j1 = j0 + 1;

    const int k1 = V.rangeI(j0).firstIndex();
    const int k2 = V.rangeI(j0).lastIndex();

    for (int k=k1; k<=k2; ++k) {
        Index1D  lambda(j0, k, XBSpline);

        cout << lambda << ": " << rhsIntegral(lambda) << endl;
    }

    for (int j=j0; j<=j1; ++j) {
        const int k1 = V.rangeJ(j).firstIndex();
        const int k2 = V.rangeJ(j).lastIndex();

        for (int k=k1; k<=k2; ++k) {
            Index1D  lambda(j, k, XWavelet);

            cout << lambda << ": " << rhsIntegral(lambda) << endl;
        }
    }
}