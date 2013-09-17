#include <iris2/iris2.cxx>

using namespace lawa;
using namespace std;

typedef double                                      T;

int
main()
{

    CompoundBasis<T>   U(2,4);
    
    int j0 = U.j0();

    int kMin, kMax;

    cout << "Spline" << endl;

    kMin = U.minK(j0, XBSpline, T(0));
    kMax = U.maxK(j0, XBSpline, T(1));

    for (int k=kMin; k<=kMax; ++k) {
        cout << "k = " << k << endl;

        cout << "supp = " << U.support(j0, k, XBSpline) << endl;
        cout << "singularSupport = "
            << U.singularSupport(j0, k, XBSpline)
            << endl;
    }


    cout << endl << endl << "Wavelet" << endl;

    kMin = U.minK(j0+1, XWavelet, T(0));
    kMax = U.maxK(j0+1, XWavelet, T(1));

    for (int k=kMin; k<=kMax; ++k) {
        cout << "k = " << k << endl;

        cout << "supp = " << U.support(j0+1, k, XWavelet) << endl;
        cout << "singularSupport = "
            << U.singularSupport(j0+1, k, XWavelet)
            << endl;
    }
}
