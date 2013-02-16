#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/finance/operators/cgmyoperator1d.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;

typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;

typedef FinanceOperator1D<T, CGMY, PrimalBasis>                     CGMYOp;
typedef NoPreconditioner<T, Index1D>                                Preconditioner;

template <typename T>
struct RhsCGMY {
    int                  nr;
    const Kernel<T,CGMY> &kernel;

    RhsCGMY(int _nr, const Kernel<T,CGMY>& _kernel);

    T
    operator()(int j, int k, XType e);
};

int main()
{
    /// wavelet basis parameters:
    int d = 2;          // (d,d_)-wavelets
    int d_ = 2;
    int j0 = 2;         // minimal level
    int J = 2;          // maximal level

    Parameters<T,CGMY> parameters(0., 1., 2., 3., 1.5);

    PrimalBasis basis(d, d_, j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    Assembler1D<T, PrimalBasis> assembler(basis);




    //Implement deltas.
    //int j = 3, k = 4;
    //FullColMatrixT deltas = computeDeltas<T,PrimalBasis>(basis, j, k, XBSpline);
    //cout << deltas << endl;

    //Implement Kernel
    //Kernel<T,CGMY> kernel(parameters);

     //Implement CGMY-operator
    //CGMYOp          a(basis, parameters);
    //SparseMatrixT   A = assembler.assembleStiffnessMatrix(a, J);
    //flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
    //densify(cxxblas::NoTrans,A,A_dense);
    //cout << A_dense << endl;

    //Implement RHS for CGMY operator
    //RHSCGMY(1, kernel);

/*
    Preconditioner  p;

    DiagonalMatrixT P = assembler.assemblePreconditioner(p, J);
    DenseVectorT    f = assembler.assembleRHS(rhs, J);

    /// Initialize empty solution vector
    DenseVectorT u(basis.mra.rangeI(J));

    /// Solve problem using pcg
    cout << pcg(P, A, u, f) << " pcg iterations" << endl;


*/
    return 0;

}
