#ifndef IRIS2_MY_SYMMETRICAPPLY1D_TCC
#define IRIS2_MY_SYMMETRICAPPLY1D_TCC 1

#include <iris2/my/compoundbasis.h>

namespace lawa {

template <typename T, typename MA>
SymmetricApply1D<T,MA>::SymmetricApply1D(const ParametersLaplace1D<T> &_params,
                                         MA    &_A)
    : parameters(_params), A(_A)
{
}

template <typename T, typename MA>
typename SymmetricApply1D<T,MA>::CoefficientsLex
SymmetricApply1D<T,MA>::operator()(const CoefficientsLex &v, int k)
{
    CoefficientsLex ret;

    int d = A.a.d();

    if (v.size()==0) {
        return ret;
    }

    CoefficientsAbs tmp;
    tmp = v;

    int s     = 0;
    int count = 0;

    for (CoeffAbsIt it=tmp.begin(); (it!=tmp.end()) && (s<=k); ++it) {
        IndexSet<Index1D> Lambda_v;

        Lambda_v = lambdaTilde1d((*it).second,
                                 A.a,               // bilinear form
                                 (k-s),
                                 A.a.j0,
                                 (*it).second.j+(k-s)+1);

        for (IndexSetIt mu = Lambda_v.begin(); mu!=Lambda_v.end(); ++mu) {
            ret[*mu] += A(*mu, (*it).second) * (*it).first;
        }
        ++count;
        s = int(log(T(count))/log(T(2))) + 1;
    }

    return ret;
}

template <typename T, typename MA>
typename SymmetricApply1D<T,MA>::CoefficientsLex
SymmetricApply1D<T,MA>::operator()(const CoefficientsLex &v, T eps)
{
    CoefficientsLex ret;
    CoefficientsAbs v_abs;

    v_abs = v;
    ret   = operator()(v, findK(v_abs, eps));
    return ret;
}

template <typename T, typename MA>
int
SymmetricApply1D<T,MA>::findK(const CoefficientsAbs &v, T eps)
{
    using std::min;
    using std::pow;

    int d = A.a.d();

    if (v.size() == 0) {
        return 1;
    }

    //s = gamma-1, gamma the smoothness index of the wavelet basis
    T s = d - T(1.5);
    T tau = T(1) / (s + T(0.5));

    // here the constant in (7.27) (KU-Wavelet) is estimated with 10
    int k_eps = static_cast<int>(10*log(pow(eps, -1.0/s)*pow(v.wtauNorm(tau), 1.0/s)) / log(2.0));
    DenseVector<Array<T> > normsec = v.norm_sections();
    T ErrorEstimateFactor = T(1);
    //std::cout << "eps = " << eps << ", k_eps = " << k_eps << std::endl;

    for (int k=1; k<=k_eps; ++k) {
        //std::cout << "At k = " << setw(3) << k;

        T R_k = 0.0;
        for (int i=k; i<=normsec.lastIndex()-1; ++i) {
            R_k += normsec(i+1);
        }
        R_k *= parameters.CA;
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k += pow(2.,-k*s) * normsec(1);
        //std::cout << ", R_k = " << setw(10) << R_k;

        for (int l=0; l<=k-1; ++l) {
            if (k-l<=normsec.lastIndex()-1) {
                //R_k += std::pow(l,-1.01)*std::pow(2.,-l*s) * normsec(k-l+1);
                R_k += pow(2.,-l*s) * normsec(k-l+1);
            }
        }
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k *= ErrorEstimateFactor;
        //std::cout << ", R_k = " << setw(10) << R_k << ", eps = " << setw(10) << eps << endl;

        if (R_k<=eps) {
            std::cout << "   findK ==> k = " << k << ", k_eps = " << k_eps << std::endl;
            int maxlevel=22;
            if (d==2)         {    maxlevel=25; }
            else if (d==3)    {   maxlevel=19; }    //for non-singular examples, also lower values are possible
            return min(k,maxlevel);
        }
    }
    return min(k_eps,25);    //higher level differences result in translation indices that cannot be stored in int.
}

}   //namespace lawa

#endif // IRIS2_MY_SYMMETRICAPPLY1D_TCC
