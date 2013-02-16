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
#include <fstream>
#include <sstream>
#include <vector>
#include <lawa/lawa.h>
#include <tutorials/examples/referencesolutions/spacetimetensorrefsols1d.h>


#define ROW_SIZE 4*4096
#define COL_SIZE 2*4*2048

typedef double T;
using namespace lawa;
using namespace std;

// Basis definitions
typedef Basis<T,Primal,Periodic,CDF>                        PeriodicBasis;
typedef Basis<T,Primal,Interval,Dijkema>                    IntervalBasis;
typedef TensorBasis2D<Adaptive,PeriodicBasis,IntervalBasis> Basis2D;

// Operator definitions
typedef LeftNormPreconditioner2D<T,Basis2D>                             LeftPrec2D;
typedef RightNormPreconditioner2D_c<T,Basis2D>                          RightPrec2D;
typedef AdaptiveSpaceTimePDEOperator1D<T, Basis2D, LeftPrec2D,
                                RightPrec2D, NoInitialCondition>        AdaptSpaceTimePDEOp;

// Righthandsides definitions
typedef SeparableRHS2D<T,Basis2D >                                      SepRhs2D;
typedef SumOfTwoRHSIntegrals<T, Index2D, SepRhs2D, SepRhs2D>            SumOf2SepRhs2D;
typedef RHS<T,Index2D, SumOf2SepRhs2D, LeftPrec2D>                      AdaptRHS;

// Algorithm definition
typedef S_ADWAV<T,Index2D, Basis2D, AdaptSpaceTimePDEOp, AdaptRHS>      S_Adwav;

// FLENS definitions
typedef GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >             FullColMatrixT;

int main (int argc, char *argv[]) {
    
    /* PARAMETERS: */

    if (argc != 7) {
        cerr << "Usage " << argv[0] << " d d_ j0_t j0_x max_its ex" << endl; 
        exit(1);
    }

    int d       = atoi(argv[1]);
    int d_      = atoi(argv[2]);
    int j0_t    = atoi(argv[3]);
    int j0_x    = atoi(argv[4]);
    int numIts  = atoi(argv[5]);
    int example = atoi(argv[6]);
    
    // S-Adwav parameters
    T contraction   = 0.125;
    T threshTol     = 0.1;
    T cgTol         = 1e-6;
    T resTol        = 1e-4;

    // Basis Initialization 
    PeriodicBasis   basis_t(d, d_, j0_t);
    IntervalBasis   basis_x(d,d_,j0_x);
    basis_x.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2D(basis_t, basis_x);

    // Operator Initialization
    T c = 1.;
    LeftPrec2D    P_left(basis2D);
    RightPrec2D   P_right(basis2D);
    AdaptSpaceTimePDEOp  A(basis2D, P_left, P_right, c);

    // Righthandside Initialization
    SpaceTimeTensorRefSols1D<T, Basis2D> refsol;
    refsol.setExample(example, c, 0, 0);
    //      RHS for tensor solutions has structure f = rhs_contrib_t * u_x + u_t * rhs_contrib_x 
    //          (see also spacetimetensorrefsols.h)
    //      Initialization of separable functions for each summand
    SeparableFunction2D<T> SepFunc1(refsol.rhs_contrib_t,   refsol.sing_pts_t, 
                                    refsol.u_x,             refsol.sing_pts_x);  
    SeparableFunction2D<T> SepFunc2(refsol.u_t,             refsol.sing_pts_t, 
                                    refsol.rhs_contrib_x,   refsol.sing_pts_x); 
    //      Initialization of integral over these functions for each summand
    FullColMatrixT no_deltas;
    SepRhs2D RHS_1(basis2D, SepFunc1, no_deltas, no_deltas, 10);
    SepRhs2D RHS_2(basis2D, SepFunc2, no_deltas, no_deltas, 10);
    //      Sum of both integrals
    SumOf2SepRhs2D rhsintegral2d(RHS_1,RHS_2);
    //      Preconditioned integral evaluations as needed by S-Adwav
    AdaptRHS F(rhsintegral2d,P_left);
    
    // S-Adwav
    //      Initialization of index set
    IndexSet<Index2D>   InitialLambda;
    for (int k_t = basis_t.mra.rangeI(j0_t).firstIndex(); k_t <= basis_t.mra.rangeI(j0_t).lastIndex(); ++k_t) {
        for (int k_x = basis_x.mra.rangeI(j0_x).firstIndex(); k_x <= basis_x.mra.rangeI(j0_x).lastIndex(); ++k_x) {
            Index1D index_t(j0_t,k_t,XBSpline);
            Index1D index_x(j0_x,k_x,XBSpline);
            InitialLambda.insert(Index2D(index_t,index_x));
        }
    }
    //      Solving the problem, using gmres due to nonsymmetry
    S_Adwav s_adwav(basis2D, A, F, contraction, threshTol, cgTol, resTol, numIts, 4, 1e-2);
    s_adwav.solve_gmres(InitialLambda);

    // Post-processing
    cout << "Starting post-processing..." << endl;
    
    //      Calculate maximal level used
    Coefficients<Lexicographical,T,Index2D> u_ref = s_adwav.solutions[numIts-1];
    typedef Coefficients<Lexicographical,T,Index2D>::const_iterator const_it;
    typedef Coefficients<Lexicographical,T,Index2D>::value_type val_type;

    int J_x = 0;
    int J_y = 0;
    for(const_it it = u_ref.begin(); it != u_ref.end(); ++it){
        if((*it).first.index1.j > J_x){ J_x = (*it).first.index1.j; }
        if((*it).first.index2.j > J_y){ J_y = (*it).first.index2.j; }
    }
    cout << "J_t = " << J_x << ", J_x = " << J_y << endl;
    
    //      Expand Index Set twice as approximation for an exact solution
    cout << "Calculate exact solution..." << endl;        
    Coefficients<Lexicographical,T,Index2D> u_exact;
    u_exact = s_adwav.solutions[numIts-1];
    IndexSet<Index2D> ExpandedLambda;
    ExpandedLambda = supp(u_exact);
    ExpandedLambda = ExpandedLambda + C(ExpandedLambda, contraction, basis2D);
    ExpandedLambda = ExpandedLambda + C(ExpandedLambda, contraction, basis2D);
    cout << "Size of expanded Lambda: " << ExpandedLambda.size() << endl;
    FillWithZeros(ExpandedLambda,u_exact);
    Coefficients<Lexicographical,T,Index2D> f = F(ExpandedLambda);
    T tmp = 0.0;
    cout << GMRES_Solve(ExpandedLambda, A, u_exact, f, tmp, cgTol) << " Gmres iterations for exact solution" << endl;
    cout << "Residuum of exact solution: " << tmp << endl;
    
    cout << "Size of u_exact: " << supp(u_exact).size() << endl;
    
    A.clear();
    ExpandedLambda.clear();
    
    //      Calculate errors
    cout << "Calculate errors..." << endl;
    Coefficients<Lexicographical,T,Index2D> u;
    string filename, filename_indexset;
    switch(example){
        case 1: filename = "ex1_s-adwav-heat-conv";
                break;
        case 2: filename = "ex2_s-adwav-heat-conv";
                break;
        case 3: filename = "ex3_s-adwav-heat-conv";
            break;
        default: break;
    }

    string convFilename = filename + ".txt";
    ofstream file(convFilename.c_str());
    file << "N TotalTime Residuals   L2(0,T; L2)   L2(0,T; H1)   W(0,T)" << endl;

    for (int i=0; i<numIts; ++i) {
        u = s_adwav.solutions[i];
        cout << "   Iteration " << i+1 << ", Lambda.size() = " << supp(u).size() << endl;
        
        T Error_W0T = estimate_SpaceTimeError_W0T(u, u_exact, P_right);
        T Error_L2_H1 = estimate_SpaceTimeError_L0T_H1(u, u_exact, P_right);
        T Error_L2_L2 = estimate_SpaceTimeError_L0T_L2(u, u_exact, P_right);
                
        file << supp(u).size() << " " << s_adwav.times[i] << " "  << s_adwav.residuals[i] 
             << " " << Error_L2_L2 << " " << Error_L2_H1 << " " << Error_W0T << endl;
        //      Plot solution in last iteration
        if (i==numIts-1) {
            cout << "... plotting" << endl;
            stringstream sizeString;
            sizeString << i;
            string fullFilename = filename + "_Data_Iter" + sizeString.str();
            plot2D(basis2D, u, P_right, refsol.sol, 0., 1., 0., 1., pow2i<T>(-6), fullFilename.c_str());
        }
        
    }

    cout << "... finished" << endl;

    return 0;

}
