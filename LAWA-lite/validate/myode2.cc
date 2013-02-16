#include <cmath>
#include <fstream>
#include <iostream>

#include <lawa/lawa.h>

#include <validate/mybasis.h>
#include <validate/mybasis.tcc>

#include <validate/myintegral.h>
#include <validate/myintegral.tcc>

#include <validate/myrhsintegral.h>
#include <validate/myrhsintegral.tcc>


using namespace lawa;
using namespace std;

typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
typedef flens::DenseVector<Array<double> >               RealDenseVector;
typedef flens::DenseVector<Array<int> >                  IntegerDenseVector;



struct MyOperator
{
    MyOperator(int d, int d_, int J, int _jMax=8)
        : U(d, d_), V(d, d_),
          j0(std::min(U.j0(), V.j0())), j1(j0+J), jMax(_jMax),
          KU(U.getFirstAbsoluteIndex(j0-1), U.getLastAbsoluteIndex(j1)),
          KV(V.getFirstAbsoluteIndex(j0-1), V.getLastAbsoluteIndex(j1))
    {
        cerr << "U.j0() = " << U.j0() << endl;
        cerr << "V.j0() = " << V.j0() << endl;
        cerr << "j0 = " << j0 << endl;

        // V.basisRight.enforceBoundaryCondition<DirichletBC>();
    }

    double
    operator()(int row, int col) const
    {
        static MyIntegral<double> myIntegral(U,V);

        double value = 0;

        int k0 = V.getFirstAbsoluteIndex(j0-1);
        int k1 = V.getLastAbsoluteIndex(jMax);

        /*
        for (int k=k0; k<=k1; ++k) {
            value += (-myIntegral(row,0,k,1) + myIntegral(row,0,k,0))
                   * (-myIntegral(col,0,k,1) + myIntegral(col,0,k,0));
        }
        return value;
        */
        
        return -myIntegral(row,0,col,1) + myIntegral(row,0,col,0);
    }

    int
    firstRow() const
    {
        return U.getFirstAbsoluteIndex(j0-1);
    }

    int
    lastRow() const
    {
        return U.getLastAbsoluteIndex(j1);
    }

    int
    numRows() const
    {
        return lastRow() - firstRow() + 1;
    }

    int
    firstCol() const
    {
        return firstRow();
    }

    int
    lastCol() const
    {
        return lastRow();
    }

    int
    numCols() const
    {
        return numRows();
    }

    MyBasis<double>     U, V;
    int                 j0, j1, jMax;
    Range<int>          KU, KV;
};




struct MyRHS
{
    MyRHS(int d, int d_, int J)
    {
    }

};


const int    k  = 2;
const double u0 = 0;

double
g(double t)
{
//    return u0;
    return sin(2*M_PI*k*t) + 2*M_PI*k*cos(2*M_PI*k*t) + u0;
    //return -sin(t) + cos(t);
    //return cos(t) + sin(t);
}

int
main()
{

    MyOperator  B(2, 4, 2);

    for (int i=1; i<=14; ++i) {
        for (int j=1; j<=12; ++j) {
            cout.width(10);
            cout << B(i,j) << "  ";
        }
        cout << endl;
    }

    /*
    for (int i=B.firstRow(); i<=B.lastRow(); ++i) {
        for (int j=B.firstCol(); j<=B.lastCol(); ++j) {
            cout << B(i,j) << " ";
        }
        cout << endl;
    }
    */

}
