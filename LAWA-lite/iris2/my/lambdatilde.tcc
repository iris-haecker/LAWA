#ifndef IRIS2_MY_LAMBDATILDE_TCC
#define IRIS2_MY_LAMBDATILDE_TCC 1

#include <iris2/my/lambdatilde.h>

namespace lawa {

template <typename T>
IndexSet<Index1D>
lambdaTilde1d(const Index1D &lambda, const Laplace1D<T> &A, 
              int sTilde, int jMin, int jMax)
{
    using std::abs;
    using std::max;
    using std::min;

    IndexSet<Index1D> ret;

    int j1 = lambda.j;
    int k1 = lambda.k;

    // TODO: remove
    jMax = std::min(16, jMax);


    if (lambda.xtype == XBSpline) {

        assert(j1==jMin);

        //
        //  (XBSpline,XBSpline)
        //
        int minK2 = A.minK2(j1, k1, XBSpline, j1, XBSpline);
        int maxK2 = A.maxK2(j1, k1, XBSpline, j1, XBSpline);

        for (int k2=minK2; k2<=maxK2; ++k2) {
            ret.insert(Index1D(j1,k2,XBSpline));
        }

        //
        //  (XBSpline, XWavelet)
        //
        for (int j2=j1; j2<=min(j1+sTilde, jMax); ++j2) {
            int minK2 = A.minK2(j1, k1, XBSpline, j2, XWavelet);
            int maxK2 = A.maxK2(j1, k1, XBSpline, j2, XWavelet);

            for (int k2=minK2; k2<=maxK2; ++k2) {
                ret.insert(Index1D(j2,k2,XWavelet));
            }
        }

    } else if (lambda.xtype == XWavelet) {

        if (j1-sTilde <= jMin) {

        //
        //  (XWavelet, XBSpline)
        //
            int minK2 = A.minK2(j1, k1, XWavelet, jMin, XBSpline);
            int maxK2 = A.maxK2(j1, k1, XWavelet, jMin, XBSpline);

            for (int k2=minK2; k2<=maxK2; ++k2) {
                ret.insert(Index1D(jMin,k2,XBSpline));
            }
        }

        //
        //  (XWavelet, XWavelet)
        //
        for (int j2=max(j1-sTilde,jMin); j2<=min(j2+sTilde,jMax); ++j2) {
            int minK2 = A.minK2(j1, k1, XBSpline, j2, XWavelet);
            int maxK2 = A.maxK2(j1, k1, XBSpline, j2, XWavelet);

            for (int k2=minK2; k2<=maxK2; ++k2) {
                ret.insert(Index1D(j2,k2,XWavelet));
            }
        }
    }
    return ret;
}

} // namespace lawa

#endif // IRIS2_MY_LAMBDATILDE_TCC