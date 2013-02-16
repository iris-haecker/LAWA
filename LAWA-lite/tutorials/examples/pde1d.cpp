/* POISSON PROBLEM 1D
 *
 *  This example calculates a poisson problem with constant forcing f on the
 *  one-dimensional domain [0,1], i.e.
 *          - u'' = f on (0,1) , u(0) = u(1) = 0.
 *  The solution is obtained using a uniform Wavelet-Galerkin method with a
 *  diagonal scaling preconditioner.
 */

/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Typedefs for problem components:
///     Primal Basis over an interval, using Dijkema construction
typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
///     PDE-Operator in 1D, i.e. for $a(v,u) = \int(v_x \cdot u_x) \int(b(x) v \cdot u') + \int(a(x) v \cdot u)$
typedef PDENonConstCoeffOperator1D<T, PrimalBasis>                  PDEOp;
///     Preconditioner: diagonal scaling with norm of operator
typedef H1NormPreconditioner1D<T, PrimalBasis>                      NormPrec;
///     Right Hand Side (RHS): basic 1D class for rhs integrals of the form $\int(f \cdot v)$,
///     possibly with additional peak contributions (not needed here)
typedef RHSWithPeaks1D<T, PrimalBasis>                              Rhs;


/// Reference solution
T
u_f(T x)
{
    return x*(x-1)*std::exp(x);
}

T
u_x_f(T x)
{
    return (x*x+x-1.)*std::exp(x);
}

T
u_xx_f(T x)
{
    return (x*x+3*x)*std::exp(x);
}

/// Non-constant reaction term
T
a_f(T x)
{
    return 1+exp(x);
}

/// Non-constant convection term
T
b_f(T x)
{
    if (x<=0.4) {
        return 1.;
    }
    else {
        return -1.;
    }
}

/// Forcing function of the form `T f(T x)` - here a constant function
T
rhs_f(T x)
{
    return -a_f(x)*u_xx_f(x) + b_f(x)*u_x_f(x) + a_f(x)*u_f(x);
}

/// Auxiliary function to print solution values, generates `.txt`-file with
/// columns: `x u(x)`
void
printU(const DenseVectorT u, const PrimalBasis& basis, const int J,
       const char* filename, const double deltaX=1./128.)
{
    ofstream file(filename);
    for(double x = 0; x <= 1.; x += deltaX){
        file << x << " " << u_f(x) << " " << evaluate(basis,J, u, x, 0) << endl;
    }
    file.close();
}

/// Auxiliary function to print solution values, generates `.txt`-file with
/// columns: `x u(x)`
void
H1errorU(const DenseVectorT u, const PrimalBasis& basis, const int j, T &L2error, T &H1error,
         const double deltaX=1./128.)
{
    L2error = 0.;
    H1error = 0.;
    T H1seminormerror = 0.;
    for(double x = 0; x <= 1.; x += deltaX) {
        T diff_u   = u_f(x)   - evaluate(basis,j, u, x, 0);
        T diff_u_x = u_x_f(x) - evaluate(basis,j, u, x, 1);
        L2error += diff_u*diff_u;
        H1seminormerror += diff_u_x*diff_u_x;
    }
    H1error = L2error + H1seminormerror;
    L2error *= deltaX;
    L2error = std::sqrt(L2error);
    H1error *= deltaX;
    H1error = std::sqrt(H1error);

}

int main()
{
    /// wavelet basis parameters:
    int d = 3;          // (d,d_)-wavelets
    int d_ = 3;
    int j0 = 3;         // minimal level
    int J = 8;          // maximal level

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, d_, j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    /// Operator initialization
    DenseVectorT a_singPts;
    DenseVectorT b_singPts(1); b_singPts = 0.5;
    Function<T> a(a_f, a_singPts);
    Function<T> b(b_f, b_singPts);
    PDEOp       op(basis, a, b, a);
    NormPrec    p(basis);

    /// Righthandside initialization
    DenseVectorT singPts;                      // singular points of the rhs forcing function: here none
    FullColMatrixT deltas;                     // peaks (and corresponding scaling coefficients): here none
    Function<T> F(rhs_f, singPts);             // Function object (wraps a function and its singular points)
    Rhs rhs(basis, F, deltas, 10, false, true); // RHS: specify integration order for Gauss quadrature (here: 4) and
                                               //      if there are singular parts (false) and/or smooth parts (true)
                                               //      in the integral

    ofstream file("pde1d_convergence.txt");
    for (int j=j0; j<=J; ++j) {
        /// Assembler: assemble the problem components
        Assembler1D<T, PrimalBasis> assembler(basis);
        SparseMatrixT   A = assembler.assembleStiffnessMatrix(op, j);
        DiagonalMatrixT P = assembler.assemblePreconditioner(p, j);
        DenseVectorT    f = assembler.assembleRHS(rhs, j);

        /// Initialize empty solution vector
        DenseVectorT u(basis.mra.rangeI(j));

        /// Solve problem using pcg
        cout << pgmres(P, A, u, f) << " pgmres iterations" << endl;

        /// Compute errors
        T L2error=0, H1error=0.;
        H1errorU(u, basis, j, L2error, H1error, pow2i<T>(-j-2));
        file << basis.mra.cardI(j) << " " << L2error << " " << H1error << endl;

        /// Print solution to file "u.txt"
        printU(u, basis, j, "u.txt");
    }
    file.close();
    return 0;
}
