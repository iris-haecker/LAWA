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
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/referencesolutions.h>

#define ROW_SIZE 4*8192
#define COL_SIZE 4*2048

typedef double T;
using namespace lawa;
using namespace std;

const DomainType domain = R;
const Construction cons = CDF;

//Basis definitions
typedef Basis<T,Primal,domain,cons> ReallineBasis;
typedef TensorBasis2D<Adaptive, ReallineBasis,ReallineBasis> Basis2D;

//Operator definitions
typedef HelmholtzOperator2D<T, Basis2D>                             HelmholtzBilinearForm2D;
//typedef DiagonalMatrixPreconditioner2D<T,Basis2D,
//                                       HelmholtzBilinearForm2D >    Preconditioner2D;

typedef H1NormPreconditioner2D<T,Basis2D>                           Preconditioner2D;

typedef AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner2D>   MA;

//Righthandsides definitions (tensor)
typedef SeparableRHS2D<T,Basis2D >                                  SeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                SumOfSeparableRhsIntegral2D;
typedef RHS<T,Index2D, SumOfSeparableRhsIntegral2D,
            Preconditioner2D>                                       SumOfSeparableRhs;

//Righthandsides definitions (non tensor, smooth rhs)
typedef SmoothRHSWithAlignedSing2D<T, Basis2D, FullGridGL>          NonSeparableRhsIntegralFG2D;
typedef RHS<T,Index2D, NonSeparableRhsIntegralFG2D,
            Preconditioner2D>                                       NonSeparableRhsFG2D;

//Righthandsides definitions (non tensor, non smooth rhs)
typedef SmoothRHSWithAlignedSing2D<T, Basis2D, FullGridGL>          NonSeparableRhsIntegralFG2D;
typedef SumOfThreeRHSIntegrals<T, Index2D,
                               NonSeparableRhsIntegralFG2D>         SumOfNonSeparableRhsIntegralFG2D;
typedef RHS<T,Index2D, SumOfNonSeparableRhsIntegralFG2D,
            Preconditioner2D>                                       SumOfNonSeparableRhsFG2D;


//Algorithm definition
typedef S_ADWAV<T,Index2D, Basis2D, MA, SumOfSeparableRhs>          S_Adwav_Tensor;
typedef S_ADWAV<T,Index2D, Basis2D, MA, NonSeparableRhsFG2D>        S_Adwav_NonTensor1;
typedef S_ADWAV<T,Index2D, Basis2D, MA, SumOfNonSeparableRhsFG2D>   S_Adwav_NonTensor2;

void
estimateMinimalLevel(int example, int d, int d_, int &jmin_x, int &jmin_y);

int main (int argc, char *argv[]) {
    if (argc!=5 && argc !=7) {
        cout << "usage " << argv[0] << " d d_ max_its example [jmin_x jmin_y]" << endl; exit(1);
    }
    cout.precision(16);
    T contraction = 0.125;
    T threshTol = 0.4;
    T cgTol = 0.1*threshTol;//1e-12;
    T resTol=1e-4;

    int d=atoi(argv[1]), d_=atoi(argv[2]);
    int NumOfIterations=atoi(argv[3]);
    int example=atoi(argv[4]);
    int jmin_x, jmin_y;
    if (argc==5) {
        estimateMinimalLevel(example, d, d_, jmin_x, jmin_y);
    }
    else {
        jmin_x=atoi(argv[5]);
        jmin_y=atoi(argv[6]);
    }

    ReallineBasis basis_x(d,d_,jmin_x);
    ReallineBasis basis_y(d,d_,jmin_y);

    Basis2D basis2d(basis_x,basis_y);
    HelmholtzBilinearForm2D Bil(basis2d, 1.);
    //Preconditioner2D P(Bil);
    Preconditioner2D P(basis2d);
    MA A(basis2d, 1., P, 1e-12, 8092, 8092);

    IndexSet<Index2D> InitialLambda;
    Index1D index_x(jmin_x,0,XBSpline);
    Index1D index_y(jmin_y,0,XBSpline);
    InitialLambda.insert(Index2D(index_x,index_y));

    //Righthand side construction for tensor solution
    if (example==1 || example==2 || example==3) {
        int order = 127;
        TensorRefSols_PDE_Realline2D<T> refsol;
        refsol.setExample(example, 1.);
        SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                        refsol.exact_y, refsol.sing_pts_y);

        SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                        refsol.rhs_y, refsol.sing_pts_y);
        GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
        SeparableRhsIntegral2D rhsintegral_x(basis2d, SepFunc1, refsol.deltas_x, no_deltas, order);
        SeparableRhsIntegral2D rhsintegral_y(basis2d, SepFunc2, no_deltas, refsol.deltas_y, order);
        SumOfSeparableRhsIntegral2D rhsintegral2d(rhsintegral_x,rhsintegral_y);
        SumOfSeparableRhs F(rhsintegral2d,P);

        S_Adwav_Tensor s_adwav(basis2d, A, F, contraction, threshTol, cgTol, resTol,
                               NumOfIterations, 3, 1e-2);
        Timer time;
        time.start();
        s_adwav.solve_cg(InitialLambda, refsol.H1norm());
        time.stop();
        cout << "S-ADWAV required " << time.elapsed() << " seconds real time" << endl;

        stringstream conv_filename;
        conv_filename << "s-adwav-realline-helmholtz2d-conv_" << example << "_" << d << "_" << d_
                      << "_" << jmin_x << "_" << jmin_y << ".dat";
        cout << "Postprocessing started." << endl;
        postprocessing_H1<T,Index2D, S_Adwav_Tensor, MA, SumOfSeparableRhs>
                         (s_adwav, A, F, refsol.H1norm(), conv_filename.str().c_str());
        cout << "Postprocessing finished." << endl;

        stringstream plot_filename;
        plot_filename << "s-adwav-realline-helmholtz2d-plot_" << example << "_" << d << "_" << d_
                      << "_" << jmin_x << "_" << jmin_y << ".dat";
        cout << "Plot of solution started." << endl;
        plot2D(basis2d, s_adwav.solutions[NumOfIterations-1], P, refsol.exact, -10., 10., -10., 10.,
               pow2i<T>(-3), plot_filename.str().c_str());
        cout << "Plot of solution finished." << endl;
    }
    //Righthand side construction for non tensor solution and smooth rhs
    else if (example==4) {
        int order = 40;
        RefSols_PDE_Realline2D<T> refsol;
        refsol.setExample(1, 1.);
        Function2D<T> Func2d(refsol.rhs, refsol.sing_pts_x, refsol.sing_pts_y);
        NonSeparableRhsIntegralFG2D rhsintegral2d(basis2d, Func2d, order);
        NonSeparableRhsFG2D F(rhsintegral2d,P);

        S_Adwav_NonTensor1 s_adwav(basis2d, A, F, contraction, threshTol, cgTol, resTol,
                                   NumOfIterations, 3, 1e-2);
        Timer time;
        time.start();
        s_adwav.solve_cg(InitialLambda, refsol.H1norm());
        time.stop();
        cout << "S-ADWAV required " << time.elapsed() << " seconds real time" << endl;

        stringstream filename;
        filename << "s-adwav-realline-helmholtz2d-conv_" << example << "_" << d << "_" << d_
                 << "_" << jmin_x << "_" << jmin_y << ".dat";
        postprocessing_H1<T,Index2D, S_Adwav_NonTensor1, MA, NonSeparableRhsFG2D>
                         (s_adwav, A, F, refsol.H1norm(), filename.str().c_str());

        stringstream plot_filename;
        plot_filename << "s-adwav-realline-helmholtz2d-plot_" << example << "_" << d << "_" << d_
                      << "_" << jmin_x << "_" << jmin_y << ".dat";
        cout << "Plot of solution started." << endl;
        plot2D(basis2d, s_adwav.solutions[NumOfIterations-1], P, refsol.exact, -10., 10., -10., 10.,
               pow2i<T>(-3), plot_filename.str().c_str());
        cout << "Plot of solution finished." << endl;
    }
    //Righthand side construction for non tensor solution and sum of non-smooth rhs
    else if (example==5) {

        RefSols_PDE_Realline2D<T> refsol;
        refsol.setExample(2, 1.);
        if (d==2) {
            int order = 20;
            Function2D<T> Func2d(refsol.exact, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_x(refsol.exact_dx, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_y(refsol.exact_dy, refsol.sing_pts_x, refsol.sing_pts_y);
            NonSeparableRhsIntegralFG2D rhsintegral_reaction(basis2d, Func2d, order);
            NonSeparableRhsIntegralFG2D rhsintegral_diffusion_x(basis2d, Func2d_x, order, 1, 0);
            NonSeparableRhsIntegralFG2D rhsintegral_diffusion_y(basis2d, Func2d_y, order, 0, 1);
            SumOfNonSeparableRhsIntegralFG2D rhsintegral2d(rhsintegral_diffusion_x,
                                                           rhsintegral_diffusion_y,
                                                           rhsintegral_reaction);
            SumOfNonSeparableRhsFG2D F(rhsintegral2d,P);

            S_Adwav_NonTensor2 s_adwav(basis2d, A, F, contraction, threshTol, cgTol, resTol,
                                       NumOfIterations, 3, 1e-2);
            Timer time;
            time.start();
            s_adwav.solve_cg(InitialLambda, refsol.H1norm());
            time.stop();
            cout << "S-ADWAV required " << time.elapsed() << " seconds real time" << endl;

            stringstream filename;
            filename << "s-adwav-realline-helmholtz2d-conv_" << example << "_" << d << "_" << d_
                     << "_" << jmin_x << "_" << jmin_y << ".dat";
            postprocessing_H1<T,Index2D, S_Adwav_NonTensor2, MA, SumOfNonSeparableRhsFG2D>
                             (s_adwav, A, F, refsol.H1norm(), filename.str().c_str());

            stringstream plot_filename;
            plot_filename << "s-adwav-realline-helmholtz2d-plot_" << example << "_" << d << "_" << d_
                          << "_" << jmin_x << "_" << jmin_y << ".dat";
            cout << "Plot of solution started." << endl;
            plot2D(basis2d, s_adwav.solutions[NumOfIterations-1], P, refsol.exact, -10., 10., -10., 10.,
                   pow2i<T>(-3), plot_filename.str().c_str());
            cout << "Plot of solution finished." << endl;
        }
        else  {
            int order = 40;
            Function2D<T> Func2d(refsol.exact, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Minus_Func2d(refsol.minus_exact, refsol.sing_pts_x, refsol.sing_pts_y);
            NonSeparableRhsIntegralFG2D rhsintegral_reaction(basis2d, Func2d, order);
            NonSeparableRhsIntegralFG2D rhsintegral_diffusion_x(basis2d, Minus_Func2d, order, 2, 0);
            NonSeparableRhsIntegralFG2D rhsintegral_diffusion_y(basis2d, Minus_Func2d, order, 0, 2);
            SumOfNonSeparableRhsIntegralFG2D rhsintegral2d(rhsintegral_diffusion_x,
                                                           rhsintegral_diffusion_y,
                                                           rhsintegral_reaction);
            SumOfNonSeparableRhsFG2D F(rhsintegral2d,P);

            S_Adwav_NonTensor2 s_adwav(basis2d, A, F, contraction, threshTol, cgTol, resTol,
                                       NumOfIterations, 3, 1e-2);
            Timer time;
            time.start();
            s_adwav.solve_cg(InitialLambda, refsol.H1norm());
            time.stop();
            cout << "S-ADWAV required " << time.elapsed() << " seconds real time" << endl;

            stringstream filename;
            filename << "s-adwav-realline-helmholtz2d-conv_" << example << "_" << d << "_" << d_
                     << "_" << jmin_x << "_" << jmin_y << ".dat";
            postprocessing_H1<T,Index2D, S_Adwav_NonTensor2, MA, SumOfNonSeparableRhsFG2D>
                             (s_adwav, A, F, refsol.H1norm(), filename.str().c_str());

            stringstream plot_filename;
            plot_filename << "s-adwav-realline-helmholtz2d-plot_" << example << "_" << d << "_" << d_
                          << "_" << jmin_x << "_" << jmin_y << ".dat";
            cout << "Plot of solution started." << endl;
            plot2D(basis2d, s_adwav.solutions[NumOfIterations-1], P, refsol.exact, -10., 10., -10., 10.,
                   pow2i<T>(-3), plot_filename.str().c_str());
            cout << "Plot of solution finished." << endl;
        }
    }

    return 0;

}

void
estimateMinimalLevel(int example, int d, int d_, int &jmin_x, int &jmin_y)
{
    ReallineBasis basis_x(d,d_,0);
    ReallineBasis basis_y(d,d_,0);
    Basis2D basis2d(basis_x,basis_y);
    HelmholtzBilinearForm2D Bil(basis2d,1.);
    //Preconditioner2D P(Bil);
    Preconditioner2D P(basis2d);

    Coefficients<Lexicographical,T,Index2D> f_LambdaWavelet, f_LambdaBSpline;
    Coefficients<AbsoluteValue,T,Index2D>   f_LambdaWavelet_abs, f_LambdaBSpline_abs;

    //Righthand side construction for tensor solution
    if (example==1 || example==2 || example==3) {
        TensorRefSols_PDE_Realline2D<T> refsol;
        refsol.setExample(example, 1.);
        SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                        refsol.exact_y, refsol.sing_pts_y);

        SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                        refsol.rhs_y, refsol.sing_pts_y);
        GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
        SeparableRhsIntegral2D rhsintegral_x(basis2d, SepFunc1, refsol.deltas_x, no_deltas, 6);
        SeparableRhsIntegral2D rhsintegral_y(basis2d, SepFunc2, no_deltas, refsol.deltas_y, 6);
        SumOfSeparableRhsIntegral2D rhsintegral2d(rhsintegral_x,rhsintegral_y);

        for (int j_x=0; j_x>=-4; --j_x) {
            for (int j_y=0; j_y>=-4; --j_y) {
                for (int k_x=-10; k_x<=10; ++k_x) {
                    for (int k_y=-10; k_y<=10; ++k_y) {
                        Index2D waveletindex(Index1D(j_x,k_x,XWavelet),Index1D(j_y,k_y,XWavelet));
                        Index2D bsplineindex(Index1D(j_x,k_x,XBSpline),Index1D(j_y,k_y,XBSpline));
                        T wavelet_val = rhsintegral2d(waveletindex)*P(waveletindex);
                        T bspline_val = rhsintegral2d(bsplineindex)*P(bsplineindex);
                        f_LambdaWavelet[waveletindex] = wavelet_val;
                        f_LambdaBSpline[bsplineindex] = bspline_val;
                    }
                }
            }
        }
    }
    else if (example==4) {
        RefSols_PDE_Realline2D<T> refsol;
        refsol.setExample(1, 1.);
        Function2D<T> Func2d(refsol.rhs, refsol.sing_pts_x, refsol.sing_pts_y);
        NonSeparableRhsIntegralFG2D rhsintegral2d(basis2d, Func2d, 4);
        int count = 0;
        for (int j_x=0; j_x>=-4; --j_x) {
            for (int j_y=0; j_y>=-4; --j_y) {
                for (int k_x=-5; k_x<=5; ++k_x) {
                    for (int k_y=-5; k_y<=5; ++k_y) {
                        Index2D waveletindex(Index1D(j_x,k_x,XWavelet),Index1D(j_y,k_y,XWavelet));
                        Index2D bsplineindex(Index1D(j_x,k_x,XBSpline),Index1D(j_y,k_y,XBSpline));
                        T wavelet_val = rhsintegral2d(waveletindex)*P(waveletindex);
                        T bspline_val = rhsintegral2d(bsplineindex)*P(bsplineindex);
                        f_LambdaWavelet[waveletindex] = wavelet_val;
                        f_LambdaBSpline[bsplineindex] = bspline_val;
                        ++count;
                        if (count %100 == 0) cout << "Estimate minimal level: " << count << endl;
                    }
                }
            }
        }
    }
    else if (example==5 || example==6) {
        RefSols_PDE_Realline2D<T> refsol;
        refsol.setExample(1, 1.);
        Function2D<T> Func2d(refsol.exact, refsol.sing_pts_x, refsol.sing_pts_y);
        Function2D<T> Func2d_x(refsol.exact_dx, refsol.sing_pts_x, refsol.sing_pts_y);
        Function2D<T> Func2d_y(refsol.exact_dy, refsol.sing_pts_x, refsol.sing_pts_y);
        NonSeparableRhsIntegralFG2D rhsintegral_reaction(basis2d, Func2d, 7);
        NonSeparableRhsIntegralFG2D rhsintegral_diffusion_x(basis2d, Func2d_x, 7, 1, 0);
        NonSeparableRhsIntegralFG2D rhsintegral_diffusion_y(basis2d, Func2d_y, 7, 0, 1);
        SumOfNonSeparableRhsIntegralFG2D rhsintegral2d(rhsintegral_diffusion_x,
                                                       rhsintegral_diffusion_y,
                                                       rhsintegral_reaction);
        int count = 0;
        for (int j_x=0; j_x>=-4; --j_x) {
            for (int j_y=0; j_y>=-4; --j_y) {
                for (int k_x=-5; k_x<=5; ++k_x) {
                    for (int k_y=-5; k_y<=5; ++k_y) {
                        Index2D waveletindex(Index1D(j_x,k_x,XWavelet),Index1D(j_y,k_y,XWavelet));
                        Index2D bsplineindex(Index1D(j_x,k_x,XBSpline),Index1D(j_y,k_y,XBSpline));
                        T wavelet_val = rhsintegral2d(waveletindex)*P(waveletindex);
                        T bspline_val = rhsintegral2d(bsplineindex)*P(bsplineindex);
                        f_LambdaWavelet[waveletindex] = wavelet_val;
                        f_LambdaBSpline[bsplineindex] = bspline_val;
                        ++count;
                        if (count %100 == 0) cout << "Estimate minimal level: " << count << endl;
                    }
                }
            }
        }
    }

    f_LambdaWavelet_abs = f_LambdaWavelet;
    f_LambdaBSpline_abs = f_LambdaBSpline;

    typedef Coefficients<AbsoluteValue,T,Index2D>::const_iterator const_it;
    const_it wavelet_estim = f_LambdaWavelet_abs.begin();
    const_it bspline_estim = f_LambdaBSpline_abs.begin();
    cout << "Wavelet estimate: " << (*wavelet_estim).second << endl;
    cout << "BSpline estimate: " << (*bspline_estim).second << endl;

    //jmin_x = (*bspline_estim).second.index1.j;
    //jmin_y = (*bspline_estim).second.index2.j;
    jmin_x = (*wavelet_estim).second.index1.j;
    jmin_y = (*wavelet_estim).second.index2.j;


}

