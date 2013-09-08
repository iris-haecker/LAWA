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

    typedef DenseVector<Array<T> >   Vector;

    IndexSet<Index1D> ret;

    int j1 = lambda.j;
    int k1 = lambda.k;

    // TODO: remove
    //jMax = std::min(16, jMax);


    if (lambda.xtype == XBSpline) {

        assert(j1==jMin);

        //
        //  (XBSpline,XBSpline)
        //
        Support<T> supp1 = A.basisRow().support(j1, k1, XBSpline);

        int minK2 = A.minK2(j1, k1, XBSpline, j1, XBSpline);
        int maxK2 = A.maxK2(j1, k1, XBSpline, j1, XBSpline);

        for (int k2=minK2; k2<=maxK2; ++k2) {
            if (overlap(supp1, A.basisCol().support(j1,k2,XBSpline)) > 0) {
                ret.insert(Index1D(j1,k2,XBSpline));
            }
        }

        //
        //  (XBSpline, XWavelet)
        //
        for (int j2=j1; j2<=min(j1+sTilde, jMax); ++j2) {

            Vector singularSupp = A.basisRow().singularSupport(j1,k1,XBSpline);

            const int i0 = singularSupp.firstIndex();
            const int i1 = singularSupp.lastIndex();

            for (int i=i0; i<=i1; ++i) {
                const T x = singularSupp(i);
                
                int minK2 = A.basisCol().minK(j2, XWavelet, x);
                int maxK2 = A.basisCol().maxK(j2, XWavelet, x);

                for (int k2=minK2; k2<=maxK2; ++k2) {
                    Support<T> supp2 = A.basisCol().support(j2,k2,XWavelet);

                    if (overlap(supp1, supp2) > 0
                     && distance(singularSupp, supp2) <= 0 )
                    {
                        ret.insert(Index1D(j2,k2,XWavelet));
                    }
                }
            }
        }

    } else if (lambda.xtype == XWavelet) {

        Support<T> supp1 = A.basisRow().support(j1, k1, XWavelet);

        if (j1-sTilde <= jMin) {

        //
        //  (XWavelet, XBSpline)
        //
            int minK2 = A.minK2(j1, k1, XWavelet, jMin, XBSpline);
            int maxK2 = A.maxK2(j1, k1, XWavelet, jMin, XBSpline);

            for (int k2=minK2; k2<=maxK2; ++k2) {
                Support<T> supp2 = A.basisCol().support(jMin,k2,XBSpline);
                if (overlap(supp1, supp2) > 0) {
                    ret.insert(Index1D(jMin,k2,XBSpline));
                }
            }
        }

        //
        //  (XWavelet, XWavelet)
        //
        for (int j2=max(j1-sTilde,jMin); j2<=min(j2+sTilde,jMax); ++j2) {

            Vector singularSupp = A.basisRow().singularSupport(j1,k1,XWavelet);

            const int i0 = singularSupp.firstIndex();
            const int i1 = singularSupp.lastIndex();

            for (int i=i0; i<=i1; ++i) {
                const T x = singularSupp(i);
                
                int minK2 = A.basisCol().minK(j2, XWavelet, x);
                int maxK2 = A.basisCol().maxK(j2, XWavelet, x);

                for (int k2=minK2; k2<=maxK2; ++k2) {
                    Support<T> supp2 = A.basisCol().support(j2,k2,XWavelet);

                    if (overlap(supp1, supp2) > 0
                     && distance(singularSupp, supp2) <= 0 )
                    {
                        ret.insert(Index1D(j2,k2,XWavelet));
                    }
                }

            }
        }
    }
    return ret;
}

} // namespace lawa

#endif // IRIS2_MY_LAMBDATILDE_TCC