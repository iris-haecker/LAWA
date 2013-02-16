/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <iostream>
#include <vector>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/refsols_pde_realline1d.h>
#include <applications/unbounded_domains/parameters/parameters.h>

typedef double T;
using namespace lawa;
using namespace std;

typedef Basis<T,Primal,R,CDF>   Basis1D;
typedef Wavelet<T,Primal,R,CDF> WaveletR;

typedef IndexSet<Index1D>::const_iterator const_set_it;

//Operator definitions
typedef HelmholtzOperator1D<T,Basis1D>                              HelmholtzBilinearForm1D;
typedef H1NormPreconditioner1D<T,Basis1D>                           Preconditioner1D;

typedef CompressionPDE1D_WO_XBSpline<T,Basis1D>                     Compression1D_WO_XBSpline;

//Righthandsides definitions
typedef RHSWithPeaks1D_WO_XBSpline<T>                               RhsIntegral1D_WO_XBSpline;
typedef SumOfTwoRHSIntegrals<T,Index1D,RhsIntegral1D_WO_XBSpline,
                             RhsIntegral1D_WO_XBSpline>             SumOfRHSIntegral1D_WO_XBSpline;
typedef RHS<T,Index1D, SumOfRHSIntegral1D_WO_XBSpline,
            Preconditioner1D>                                       Rhs_WO_XBSpline;
typedef RHS<T,Index1D, RhsIntegral1D_WO_XBSpline, Preconditioner1D> Rhs_PP_WO_XBSpline;


//Matrix definitions
typedef MapMatrix<T,Index1D,HelmholtzBilinearForm1D,
                  Compression1D_WO_XBSpline, Preconditioner1D>      MA_WO_XBSpline;

//APPLY definitions
typedef Parameters<T, Basis1D, HelmholtzBilinearForm1D,
                   Preconditioner1D>                                Parameters1D;
typedef SYM_APPLY_1D<T,Index1D,Basis1D,Parameters1D,MA_WO_XBSpline> APPLY1D_WO_XBSpline;

//Algorithm definition
typedef GHS_ADWAV1D<T,Basis1D,APPLY1D_WO_XBSpline,Rhs_WO_XBSpline>  GHS_Adwav_WO_XBSpline;

template <typename T>
IndexSet<Index1D>
computeRHSLambda_SingularPart(const WaveletR &psi, const DenseVector<Array<T> > &_f_singularPoints,
                              int J_plus, int J_minus);

template <typename T>
IndexSet<Index1D>
computeRHSLambda_SmoothPart(const WaveletR &psi, T a, T b, int J_plus, int J_minus);

template <typename T>
Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const WaveletR &psi, const RhsIntegral1D_WO_XBSpline &rhsintegral1d_singular,
                    const RhsIntegral1D_WO_XBSpline &rhsintegral1d_smooth,
                    const Preconditioner1D &P, RefSols_PDE_Realline1D<T> &refsol);


int main (int argc, char *argv[]) {
    cout.precision(8);
    if (argc != 5) {
        cout << "usage " << argv[0] << " d d_ max_its example" << endl; exit(1);
    }

    int d=atoi(argv[1]), d_=atoi(argv[2]);

    int NumOfIterations=atoi(argv[3]);

    int example=atoi(argv[4]);
    assert(example>=1); assert(example<=6);
    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,1.,0.,1.);

    int j0=0;

    T eps=1e-5;

    Basis1D basis(d,d_,0);
    WaveletR psi(basis);
    HelmholtzBilinearForm1D Bil(basis,1.);
    Preconditioner1D P(basis);
    Compression1D_WO_XBSpline Compr(basis);
    Parameters1D parameters(basis,Bil,false);


    T left_bound = 0., right_bound = 0.;
    int J_plus_smooth = 0, J_minus_smooth = 0, J_plus_singular = 0, J_minus_singular = 0;
    bool singular_integral=false;
    refsol.getRHS_WO_XBSplineParameters(d, d_, left_bound, right_bound,
                                        J_plus_smooth, J_minus_smooth,
                                        J_plus_singular, J_minus_singular, singular_integral);

    RhsIntegral1D_WO_XBSpline      rhsintegral1d_singular(psi, refsol.rhs, refsol.sing_pts,
                                                          refsol.deltas, left_bound, right_bound,
                                                          1., 100, true, false);

    RhsIntegral1D_WO_XBSpline      rhsintegral1d_smooth(psi, refsol.rhs, refsol.sing_pts,
                                                        refsol.deltas, left_bound, right_bound,
                                                        1., 100, false, true);

    SumOfRHSIntegral1D_WO_XBSpline rhsintegral(rhsintegral1d_smooth, rhsintegral1d_singular);
    Coefficients<Lexicographical,T,Index1D> f;
    f = initializeRHSVector(psi, rhsintegral1d_singular, rhsintegral1d_smooth, P, refsol);
    Rhs_WO_XBSpline F(rhsintegral,P,f);

    MA_WO_XBSpline A(Bil,P,Compr);

    APPLY1D_WO_XBSpline Apply(parameters,basis,A);

    GHS_Adwav_WO_XBSpline ghs_adwav(basis,Apply,F);
    ghs_adwav.SOLVE(f.norm(2.),eps,NumOfIterations,refsol.H1norm());


    RhsIntegral1D_WO_XBSpline rhsintegral1d_pp(psi, refsol.rhs, refsol.sing_pts,
                                               refsol.deltas, left_bound-100, right_bound+100,
                                               1., 300, true, true);
    Rhs_PP_WO_XBSpline F_pp(rhsintegral1d_pp,P);
    cout << "Postprocessing started." << endl;
    stringstream filename;
    filename << "ghs-adwav-realline-helmholtz1d-WO-XBSpline-conv_"
             << example << "_" << d << "_" << d_ << "_" << j0 << ".dat";
    ofstream file(filename.str().c_str());
    postprocessing_H1<T,Index1D, GHS_Adwav_WO_XBSpline,
                      MA_WO_XBSpline, Rhs_PP_WO_XBSpline>(ghs_adwav, A, F_pp, refsol.H1norm(),
                                                          filename.str().c_str());
    cout << "Postprocessing finished." << endl;

    stringstream plot_filename;
    plot_filename << "ghs-adwav-realline-helmholtz1d-WO-XBSpline-plot_" << example
                  << "_" << d << "_" << d_  << ".dat";
    cout << "Plot of solution started." << endl;
    plot<T, Basis1D, Preconditioner1D>(basis, ghs_adwav.solutions[NumOfIterations-1], P, refsol.u,
                                       refsol.d_u, -10., 10., pow2i<T>(-5),
                                       plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;

/*
    stringstream filename;
    filename << "adwav-ghs-realline-helmholtz1d-WO-XBSpline-conv_"
             << example << "_" << d << "_" << d_ << ".dat";
    ofstream file(filename.str().c_str());
    Coefficients<Lexicographical,T,Index1D> u;
    for (int i=0; i<NumOfIterations; ++i) {
        u = ghs_adwav.solutions[i];
        cout << "   Iteration " << i+1 << ", Lambda.size() = " << supp(u).size() << endl;
        T Error_H_energy = estimateError_H_energy(A, F_pp, u, refsol.H1norm());
        cout << "      Error ||u_h-u||  = " << Error_H_energy << endl;
        file << supp(u).size() << " " << ghs_adwav.times[i] << " " <<  ghs_adwav.residuals[i] << " "
             << Error_H_energy << endl;
        if ((u.size() >= 300) && (u.size()<=400)) {
            Coefficients<AbsoluteValue,T,Index1D > u_abs;
            u_abs = u;
            plotCoeff(u_abs, basis, "u_coeff");
        }
        if (i==NumOfIterations-1) {
            T H1norm=0.;
            plot(basis, u, P, refsol.exact, refsol.d_exact, -20.,20., pow2i<T>(-5), H1norm,
                 "adwav-ghs-realline-helmholtz1d-WO-XBSpline");
        }

    }
*/
    return 0;
}

template <typename T>
IndexSet<Index1D>
computeRHSLambda_SingularPart(const WaveletR &psi, const DenseVector<Array<T> > &_f_singularPoints,
                              int J_plus, int J_minus)
{
    T l1 = psi.support(0,0).l1, l2 = psi.support(0,0).l2;
    IndexSet<Index1D> ret;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        for (int j=0; j<=J_plus; ++j) {
            int k_left =  std::floor(float(pow2i<T>(j)*x-l2));
            int k_right = std::ceil(float(pow2i<T>(j)*x-l1));
            for (int k=k_left; k<=k_right; ++k) {
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
        for (int j=-1; j>=J_minus; --j) {
            int k_left =  std::floor(float(pow2i<T>(j)*x-l2));
            int k_right = std::ceil(float(pow2i<T>(j)*x-l1));
            for (int k=k_left; k<=k_right; ++k) {
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
computeRHSLambda_SmoothPart(const WaveletR &psi, T a, T b, int J_plus, int J_minus)
{
    T l1 = psi.support(0,0).l1, l2 = psi.support(0,0).l2;
    IndexSet<Index1D> ret;
    for (int j=0; j<=J_plus; ++j) {
        int k_left =  std::floor(float(pow2i<T>(j)*a-l2));
        int k_right = std::ceil(float(pow2i<T>(j)*b-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(j,k,XWavelet));
        }
    }
    for (int j=-1; j>=J_minus; --j) {
        int k_left =  std::floor(float(pow2i<T>(j)*a-l2));
        int k_right = std::ceil(float(pow2i<T>(j)*b-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(j,k,XWavelet));
        }
    }
    return ret;
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const WaveletR &psi, const RhsIntegral1D_WO_XBSpline &rhsintegral1d_singular,
                    const RhsIntegral1D_WO_XBSpline &rhsintegral1d_smooth,
                    const Preconditioner1D &P, RefSols_PDE_Realline1D<T> &refsol)
{
    T left_bound = 0., right_bound = 0.;
    int J_plus_smooth = 0, J_minus_smooth = 0, J_plus_singular = 0, J_minus_singular = 0;
    bool singular_integral=false;
    refsol.getRHS_WO_XBSplineParameters(psi.d, psi.d_, left_bound, right_bound,
                                        J_plus_smooth, J_minus_smooth,
                                        J_plus_singular, J_minus_singular, singular_integral);

    Coefficients<Lexicographical,T,Index1D> f_singular;
    IndexSet<Index1D> LambdaRHS_singular = computeRHSLambda_SingularPart(psi, refsol.sing_pts,
                                                                         J_plus_singular,
                                                                         J_minus_singular);

    cout << "Initializing singular part of rhs, size of indexset: "
         << LambdaRHS_singular.size() << endl;
    for (const_set_it it=LambdaRHS_singular.begin(); it!=LambdaRHS_singular.end(); ++it) {
        f_singular[*it] = P(*it)*rhsintegral1d_singular(*it);
    }

    Coefficients<Lexicographical,T,Index1D> f_smooth;
    IndexSet<Index1D> LambdaRHS_smooth   = computeRHSLambda_SmoothPart(psi, left_bound,right_bound,
                                                                       J_plus_smooth,
                                                                       J_minus_singular);

    cout << "Initializing smooth part of rhs, size of indexset: " << LambdaRHS_smooth.size() << endl;
    for (const_set_it it=LambdaRHS_smooth.begin(); it!=LambdaRHS_smooth.end(); ++it) {
        f_smooth[*it] = P(*it)*rhsintegral1d_smooth(*it);
    }
    if (singular_integral) {
        for (const_set_it it=LambdaRHS_singular.begin(); it!=LambdaRHS_singular.end(); ++it) {
            if (LambdaRHS_smooth.count(*it)>0) continue;
            f_smooth[*it] = P(*it)*rhsintegral1d_smooth(*it);
        }
    }

    return f_smooth + f_singular;
}
