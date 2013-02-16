#include <iostream>
#include <fstream>
#include <sstream>
#include <lawa/lawa.h>
#include <tutorials/examples/referencesolutions/referencesolutions.h>

using namespace std;
using namespace lawa;

typedef double T;

// Basis definitions
typedef Basis<T, Primal, Interval, Dijkema>     PrimalBasis;
typedef Basis<T, Dual, Interval, Dijkema>       DualBasis;
typedef Basis<T, Primal, Periodic, CDF>         PeriodicBasis;
typedef PDEConstCoeffOperator1D<T, PrimalBasis> PDE1DOp;

// Righthandsides definitions
//      Timedependent Rhs, as we have to be able to evaluate
//      \int_\Omega f(t,x)*v(x) = f(t) * \int_\Omega f(x)*v(x) for each t
typedef TimedepSeparableRHS1D<T, PrimalBasis>    TimedepRHS1D;
typedef SumOfTimedepRHS1D<T, TimedepRHS1D>       SumOfRHS1D;

// TimeStepping Methods
typedef ThetaScheme1D_LTI<T, PrimalBasis, PDE1DOp, SumOfRHS1D>                  Theta;
typedef TimeStepping<T, Theta>                                                  TimeStepper;
typedef FixedPointSolver<T, TimeStepper>                                        ThetaFPSolver;
typedef MultiGrid_2ndKind_LTI<T, PrimalBasis, DualBasis, PDE1DOp, SumOfRHS1D >  Multigrid;


typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

const double c = 1.;

int main(int argc, char* argv[]){
    
    if(argc != 7){
        cerr << "Usage: " << argv[0] << " j0 J d d_ steps example" << endl;
        return 0;
    }
    
    /* PARAMETERS: minimal level etc */
    
    int j0 = atoi(argv[1]);
    int J = atoi(argv[2]);
    int d = atoi(argv[3]);
    int d_ = atoi(argv[4]);
    int steps = atoi(argv[5]);
    int ex = atoi(argv[6]);
    T theta = 0.5;    
    T timestep = 1./steps;
    
    T c = 1.;   // diffusion constant
    T k = 0.;   // convection constant
    T r = 0.;   // reaction constant
    
    // Basis Initialization
    PrimalBasis basis1d(d,d_,j0);
    DualBasis   dualbasis1d(d,d_,j0);
    basis1d.enforceBoundaryCondition<DirichletBC>();
    dualbasis1d.enforceBoundaryCondition<DirichletBC>();
    cout << "Basis: " << basis1d.mra.cardI(J) << endl;
    
    //  Operator Initialization
    Assembler1D<T, PrimalBasis> assembler(basis1d);
    PDE1DOp A(basis1d, r, k, c);
    
    // Right Hand Side Initialization
    //      Use Reference Solutions
    SpaceTimeTensorRefSols1D<T,TensorBasis2D<Uniform, PeriodicBasis, PrimalBasis> > refsol;
    refsol.setExample(ex, c, k, r);
    SeparableFunction2D<T> rhs_fct_1(refsol.rhs_contrib_t, refsol.sing_pts_t,
                                     refsol.u_x, refsol.sing_pts_x);
    SeparableFunction2D<T> rhs_fct_2(refsol.u_t, refsol.sing_pts_t,
                                     refsol.rhs_contrib_x, refsol.sing_pts_x);
    TimedepRHS1D rhs_1(basis1d, rhs_fct_1, 10);
    TimedepRHS1D rhs_2(basis1d, rhs_fct_2, 10);
    SumOfRHS1D   F(rhs_1, rhs_2);
    
    // Initialize Solvers
    DenseVectorT u_0(basis1d.mra.rangeI(j0));

    //Theta scheme(theta, basis1d, A, F);
    //TimeStepper timestepmethod(scheme, timestep, steps, J);
    //ThetaFPSolver fixedpointsolver(timestepmethod);
    Multigrid mgsolver(basis1d, dualbasis1d, A, F, theta, timestep, steps, j0);
    Timer timer;
    
    timer.start();
    DenseVectorT u = mgsolver.solve(u_0, J);
    timer.stop();    
    
    // Postprocessing
    cout << steps << " " << J << " " << basis1d.mra.cardI(J) << " " << timer.elapsed() << " seconds " << endl;
    print_U(u, basis1d, J, "FP_u_T.txt");

    cout << "L2_error at t = T = " << calculateL2Error(1., u, refsol.sol, basis1d, J, 1./pow2i<T>(J+2)) << endl;
    cout << "H1_error at t = T = " << calculateH1Error(1., u, refsol.sol, refsol.dx_sol, basis1d, J, 1./pow2i<T>(J+2)) << endl;
    
    /*FullColMatrixT& solutions = fixedpointsolver.getSolutions();
    print_U(solutions, basis1d, J, "FP_u_2d.txt", timestep, steps);

    cout << "Error in L2(0,T; L2) = " << calculateL2_L2_Error(solutions, refsol.sol, basis1d, J, timestep, steps, 1./pow2i<T>(J+2)) << endl;
    cout << "Error in L2(0,T; H1) = " << calculateL2_H1_Error(solutions, refsol.sol, refsol.dx_sol, basis1d, J, timestep, steps, 1./pow2i<T>(J+2)) << endl;
    */
    return 0;
}
