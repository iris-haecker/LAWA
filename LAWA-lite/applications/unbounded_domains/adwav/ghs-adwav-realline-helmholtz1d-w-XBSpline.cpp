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

typedef CompressionPDE1D<T,Basis1D>                                 Compression1D;

//Righthandsides definitions
typedef RHSWithPeaks1D<T,Basis1D> RhsIntegral1D;
typedef SumOfTwoRHSIntegrals<T,Index1D,RhsIntegral1D,RhsIntegral1D> SumOfRHSIntegral1D;
typedef RHS<T,Index1D, SumOfRHSIntegral1D, Preconditioner1D>        Rhs;
typedef RHS<T,Index1D, RhsIntegral1D, Preconditioner1D>             Rhs_PP;


//Matrix definitions
typedef MapMatrix<T,Index1D,HelmholtzBilinearForm1D,Compression1D,
                  Preconditioner1D>                                 MA;

//APPLY definitions
typedef Parameters<T, Basis1D, HelmholtzBilinearForm1D, Preconditioner1D> Parameters1D;
typedef SYM_APPLY_1D<T,Index1D,Basis1D,Parameters1D,MA>                   APPLY1D;

//Algorithm definition
typedef GHS_ADWAV1D<T,Basis1D,APPLY1D,Rhs>                                GHS_Adwav;

template <typename T>
IndexSet<Index1D>
computeRHSLambda_SingularPart(const Basis1D &basis, const DenseVector<Array<T> > &_f_singularPoints,
                              int J_plus);

template <typename T>
IndexSet<Index1D>
computeRHSLambda_SmoothPart(const Basis1D &basis, T a, T b, int J_plus);

template <typename T>
Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const Basis1D &basis, const RhsIntegral1D &rhsintegral1d_singular,
                    const RhsIntegral1D &rhsintegral1d_smooth, const Preconditioner1D &P,
                    RefSols_PDE_Realline1D<T> &refsol);


int main (int argc, char *argv[]) {
    cout.precision(8);
    if (argc != 5 && argc !=6) {
        cout << "usage " << argv[0] << " d d_ max_its example [j0]" << endl; exit(1);
    }

    int d=atoi(argv[1]), d_=atoi(argv[2]);

    int NumOfIterations=atoi(argv[3]);

    int example=atoi(argv[4]);
    assert(example>=1); assert(example<=6);
    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,1.,0.,1.);

    int jmin=0;
    if (argc==6) {
        jmin=atoi(argv[5]);
    }
    else {
        jmin=refsol.getMinimalLevel(d,d_);
    }

    T eps=1e-5;

    Basis1D basis(d,d_,jmin);
    HelmholtzBilinearForm1D Bil(basis,1.);
    Preconditioner1D P(basis);
    Compression1D Compr(basis);
    Parameters1D  parameters(basis,Bil,true,basis.j0);

    Function<T>        rhs_func(refsol.rhs,refsol.sing_pts);
    RhsIntegral1D      rhsintegral1d_singular(basis, rhs_func, refsol.deltas, 120, true, false);
    RhsIntegral1D      rhsintegral1d_smooth(basis, rhs_func, refsol.deltas, 120, false, true);
    SumOfRHSIntegral1D rhsintegral(rhsintegral1d_smooth, rhsintegral1d_singular);
    Coefficients<Lexicographical,T,Index1D> f;
    f = initializeRHSVector(basis, rhsintegral1d_singular, rhsintegral1d_smooth, P, refsol);
    Rhs F(rhsintegral,P,f);

    MA A(Bil,P,Compr);

    APPLY1D Apply(parameters,basis,A);

    GHS_Adwav ghs_adwav(basis,Apply,F);
    cout << "ADWAV started." << endl;
    ghs_adwav.SOLVE(f.norm(2.),eps,NumOfIterations,refsol.H1norm());
    cout << "ADWAV finished." << endl;


    RhsIntegral1D rhsintegral1d_pp(basis,rhs_func, refsol.deltas,300);
    Rhs_PP F_pp(rhsintegral1d_pp,P);

    cout << "Postprocessing started." << endl;
    stringstream filename;
    filename << "ghs-adwav-realline-helmholtz1d-W-XBSpline-conv_"
             << example << "_" << d << "_" << d_ << "_" << j0 << ".dat";
    ofstream file(filename.str().c_str());
    postprocessing_H1<T,Index1D, GHS_Adwav, MA, Rhs_PP>(ghs_adwav, A, F_pp, refsol.H1norm(),
                                                        filename.str().c_str());
    cout << "Postprocessing finished." << endl;

    stringstream plot_filename;
    plot_filename << "ghs-adwav-realline-helmholtz1d-W-XBSpline-plot_" << example
                  << "_" << d << "_" << d_ << "_" << jmin << ".dat";
    cout << "Plot of solution started." << endl;
    plot<T, Basis1D, Preconditioner1D>(basis, ghs_adwav.solutions[NumOfIterations-1], P, refsol.u,
                                       refsol.d_u, -10., 10., pow2i<T>(-5),
                                       plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;

    return 0;
}

template <typename T>
IndexSet<Index1D>
computeRHSLambda_SingularPart(const Basis1D &basis, const DenseVector<Array<T> > &_f_singularPoints,
                              int J_plus)
{
    IndexSet<Index1D> ret;
    BSpline<T,Primal,R,CDF> phi(basis.mra);
    T l1, l2;
    l1 = phi.support(0,0).l1, l2 = phi.support(0,0).l2;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        int k_left =  std::floor(float(pow2i<T>(basis.j0)*x-l2));
        int k_right = std::ceil(float(pow2i<T>(basis.j0)*x-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(basis.j0,k,XBSpline));
        }
    }

    Wavelet<T,Primal,R,CDF> psi(basis);
    l1 = psi.support(0,0).l1, l2 = psi.support(0,0).l2;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        for (int j=basis.j0; j<=J_plus; ++j) {
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
computeRHSLambda_SmoothPart(const Basis1D &basis, T a, T b, int J_plus)
{
    IndexSet<Index1D> ret;
    BSpline<T,Primal,R,CDF> phi(basis.mra);
    T l1, l2;
    l1 = phi.support(0,0).l1, l2 = phi.support(0,0).l2;
    int k_left =  std::floor(float(pow2i<T>(basis.j0)*a-l2));
    int k_right = std::ceil(float(pow2i<T>(basis.j0)*b-l1));
    for (int k=k_left; k<=k_right; ++k) {
        ret.insert(Index1D(basis.j0,k,XBSpline));
    }
    Wavelet<T,Primal,R,CDF> psi(basis);
    l1 = psi.support(0,0).l1, l2 = psi.support(0,0).l2;
    for (int j=basis.j0; j<=J_plus; ++j) {
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
initializeRHSVector(const Basis1D &basis, const RhsIntegral1D &rhsintegral1d_singular,
                    const RhsIntegral1D &rhsintegral1d_smooth, const Preconditioner1D &P,
                    RefSols_PDE_Realline1D<T> &refsol)
{
    T left_bound = 0., right_bound = 0.;
    int J_plus_smooth = 0, J_plus_singular = 0;
    bool singular_integral=false;
    refsol.getRHS_W_XBSplineParameters(basis.d, basis.d_, left_bound, right_bound, J_plus_smooth,
                                       J_plus_singular, singular_integral);


    Coefficients<Lexicographical,T,Index1D> f_singular;
    IndexSet<Index1D> LambdaRHS_singular = computeRHSLambda_SingularPart(basis,refsol.sing_pts,
                                                                         J_plus_singular);
    cout << "Initializing singular part of rhs, size of indexset: "
         << LambdaRHS_singular.size() << endl;
    for (const_set_it it=LambdaRHS_singular.begin(); it!=LambdaRHS_singular.end(); ++it) {
        f_singular[*it] = P(*it)*rhsintegral1d_singular(*it);
    }


    Coefficients<Lexicographical,T,Index1D> f_smooth;
    IndexSet<Index1D> LambdaRHS_smooth =computeRHSLambda_SmoothPart(basis,left_bound,
                                                                    right_bound,J_plus_smooth);
    cout << "Initializing smooth part of rhs, size of indexset: "
         << LambdaRHS_smooth.size() << endl;
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
