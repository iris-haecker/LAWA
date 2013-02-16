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

typedef double T;
using namespace lawa;
using namespace std;

const DomainType domain = R;
const Construction cons = CDF;

typedef Basis<T,Primal,domain,cons> Basis1D;

//Operator definitions
typedef HelmholtzOperator1D<T,Basis1D>     HelmholtzBilinearForm1D;
typedef CompressionPDE1D<T,Basis1D>        Compression1D;
typedef H1NormPreconditioner1D<T,Basis1D>  Preconditioner1D;

//Righthandsides definitions
typedef RHSWithPeaks1D<T,Basis1D > RhsIntegral1D;
typedef RHS<T,Index1D, RhsIntegral1D, Preconditioner1D> Rhs;

//MapMatrix definition
//typedef MapMatrixWithZeros<T,Index1D,HelmholtzBilinearForm1D,Compression1D,Preconditioner1D> MA;
typedef MapMatrix<T,Index1D,HelmholtzBilinearForm1D,Compression1D,Preconditioner1D> MA;

//Algorithm definition
typedef S_ADWAV<T,Index1D, Basis1D, MA, Rhs> S_Adwav;

int
estimateMinimalLevel(int example, T c, int d, int d_);

int main (int argc, char *argv[]) {
    if (argc != 5 && argc != 6) {
        cout << "usage " << argv[0] << " d d_ max_its example [jmin]" << endl; exit(1);
    }
    T c = 1.;
    T contraction = 0.125;
    T threshTol = 0.1, cgTol = 0.1*threshTol, resTol=1e-4;

    int d=atoi(argv[1]), d_=atoi(argv[2]);
    int NumOfIterations=atoi(argv[3]);
    int example=atoi(argv[4]);
    int jmin=0;
    int rhs_order=10;
    if (argc==5) {
        if (d==2 && d_==2) {
            if (example==1)      jmin=-2;
            else if (example==2) jmin=-4;
            else if (example==3) jmin=-1;
            else if (example==4) jmin=0;
            else if (example==5) jmin=0;
            else if (example==6) jmin=-1;       //better convergence behaviour
        }
        else if (d==3 && d_==3) {
            if (example==1)      jmin=-2;
            else if (example==2) jmin=-4;
            else if (example==3) jmin=-1;
            else if (example==4) jmin=0;
            else if (example==5) jmin=-1;
            else if (example==6) jmin=-1;       //better convergence behaviour
        }
        else if (d==3 && d_==5) {
            if (example==1)      jmin=-2;
            else if (example==2) jmin=-4;
            else if (example==3) jmin=-1;
            else if (example==4) jmin=0;
            else if (example==5) jmin=-1;
            else if (example==6) jmin=-1;       //better convergence behaviour
        }
        //jmin = estimateMinimalLevel(example, c, d,d_);
    }
    else {
        jmin=atoi(argv[5]);
    }
    if (example==6) rhs_order=20;

    cout << "Initializing S-ADWAV, jmin = " << jmin << endl;

    Basis1D basis(d,d_,jmin);
    HelmholtzBilinearForm1D Bil(basis,c);
    Preconditioner1D P(basis);
    Compression1D Compr(basis);
    //MA A(Bil,P,Compr,0,2*8192,2*8192);
    MA A(Bil,P,Compr);
    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example, 1., 0, c);
    Function<T> rhs(refsol.rhs,refsol.sing_pts);
    RhsIntegral1D rhsintegral1d(basis, rhs, refsol.deltas, 25);
    Rhs F(rhsintegral1d,P);

    IndexSet<Index1D> InitialLambda;
    InitialLambda.insert(Index1D(jmin,1,XBSpline));

    S_Adwav s_adwav(basis, A, F, contraction, threshTol, cgTol, resTol, NumOfIterations, 2, 1e-2);
    cout << "... finished." << endl;

    Timer time;
    time.start();
    s_adwav.solve_cg(InitialLambda, refsol.H1norm());
    time.stop();
    cout << "S-ADWAV required " << time.elapsed() << " seconds real time" << endl;

    RhsIntegral1D rhsintegral1d_postproc(basis,rhs, refsol.deltas,60);
    Rhs F_postproc(rhsintegral1d_postproc,P);
    cout << "Postprocessing started." << endl;
    stringstream filename;
    filename << "s-adwav-realline-helmholtz1d-conv_" << example << "_" << d << "_" << d_
             << "_" << jmin << ".dat";
    assert(c==1);
    postprocessing_H1<T,Index1D, S_Adwav, MA, Rhs>(s_adwav, A, F_postproc, refsol.H1norm(),
                                                   filename.str().c_str());
    cout << "Postprocessing finished." << endl;

    stringstream plot_filename;
    plot_filename << "s-adwav-realline-helmholtz1d-plot_" << example << "_" << d << "_" << d_
                  << "_" << jmin << ".dat";
    cout << "Plot of solution started." << endl;
    plot<T, Basis1D, Preconditioner1D>(basis, s_adwav.solutions[NumOfIterations-1], P, refsol.u, refsol.d_u, -10., 10.,
         pow2i<T>(-5), plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;
    return 0;
}


int
estimateMinimalLevel(int example, T c, int d, int d_)
{
    //if (example==6) {
    //  return -3;  //works also with estimated j0, but H1-errors are more stable on a lower level
    //}
    Basis1D basis(d,d_,0);
    HelmholtzBilinearForm1D Bil(basis,c);
    Preconditioner1D P(basis);
    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example, 1., 0., c);
    Function<T> rhs_func(refsol.rhs,refsol.sing_pts);
    RhsIntegral1D rhsintegral(basis, rhs_func, refsol.deltas,25);

    Coefficients<Lexicographical,T,Index1D> f_LambdaWavelet, f_LambdaBSpline;
    Coefficients<AbsoluteValue,T,Index1D>   f_LambdaWavelet_abs, f_LambdaBSpline_abs;
    for (int j=0; j>=-8; --j) {
        for (int k=-40; k<=40; ++k) {
            f_LambdaWavelet[Index1D(j,k,XWavelet)] = rhsintegral(XWavelet,j,k)*P(XWavelet,j,k);
            f_LambdaBSpline[Index1D(j,k,XBSpline)] = rhsintegral(XBSpline,j,k)*P(XBSpline,j,k);
        }
    }
    f_LambdaWavelet_abs = f_LambdaWavelet;
    f_LambdaBSpline_abs = f_LambdaBSpline;

    typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator const_it;
    const_it wavelet_estim = f_LambdaWavelet_abs.begin();
    const_it bspline_estim = f_LambdaBSpline_abs.begin();
    cout << "Wavelet estimate: " << (*wavelet_estim).second.j
         << ", BSpline estimate: " << (*bspline_estim).second.j << endl;
    if (example==3) {
        return (*wavelet_estim).second.j;
    }
    else {
        return (*bspline_estim).second.j;
    }


    return 0;
}

