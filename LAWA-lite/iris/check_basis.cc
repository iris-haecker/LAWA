#include <iostream>
#include <iris/iris.cxx>

using namespace lawa;
using namespace std;

int
main()
{
    int   d  = 3;
    int   d_ = 5;
    
    MyBasis<double>  U(d, d_);
    const int        j0 = U.j0();

    int k1 = U.basisLeft.mra.rangeI(j0).firstIndex();
    int k2 = U.basisLeft.mra.rangeI(j0).lastIndex();

    double x = 0.1875;
    cerr << "x = " << x << endl;

    cerr << "mink:" << U.basisLeft.generator(XBSpline).minK(j0, x) << endl;
    cerr << "maxk:" << U.basisRight.generator(XBSpline).maxK(j0, x) << endl;

    cerr << "XBSpline:" << endl;
    for (int k=k1; k<=k2; ++k) {
        cerr << k << " "
             << U.basisLeft.generator(XBSpline).support(j0, k) << endl;
    }


    k1 = U.basisLeft.rangeJ(j0).firstIndex();
    k2 = U.basisLeft.rangeJ(j0).lastIndex();

    cerr << "XWavelet:" << endl;
    for (int k=k1; k<=k2; ++k) {
        if (k<=U.basisLeft.rangeJL(j0).lastIndex()) {
            cerr << "LEFT   ";
        } else if (k<=U.basisLeft.rangeJI(j0).lastIndex()) {
            cerr << "INNER  ";
        } else {
            cerr << "RIGHT  ";
        }
        cerr << k << " "
             << U.basisLeft.generator(XWavelet).support(j0, k) << endl;
    }

    cerr << "x = " << x << endl;
    cerr << "mink:" << U.basisLeft.generator(XWavelet).minK(j0, x) << endl;
    cerr << "maxk:" << U.basisRight.generator(XWavelet).maxK(j0, x) << endl;

    x = 0.4;
    cerr << "x = " << x << endl;
    cerr << "mink:" << U.basisLeft.generator(XWavelet).minK(j0, x) << endl;
    cerr << "maxk:" << U.basisRight.generator(XWavelet).maxK(j0, x) << endl;

    x = 0.6;
    cerr << "x = " << x << endl;
    cerr << "mink:" << U.basisLeft.generator(XWavelet).minK(j0, x) << endl;
    cerr << "maxk:" << U.basisRight.generator(XWavelet).maxK(j0, x) << endl;

}

