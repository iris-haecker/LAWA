#include <cmath>
#include <fstream>
#include <iostream>

#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;
typedef flens::DenseVector<Array<double> >               RealDenseVector;
typedef flens::DenseVector<Array<int> >                  IntegerDenseVector;


typedef Basis<double,Primal,Interval,Dijkema>            PrimalBasis;


struct MyPrimalBasis
{
    MyPrimalBasis(int d, int d_)
        : ULeft(d, d_), URight(d, d_)
    {
    }
    
    template <BoundaryCondition BCLeft, BoundaryCondition BCRight>
    void
    enforceBoundaryCondition()
    {
        if (BCLeft==DirichletBC) {
            ULeft.enforceBoundaryCondition<DirichletBC>();
        }
        if (BCRight==DirichletBC) {
            URight.enforceBoundaryCondition<DirichletBC>();
        }
    }
    
    int
    index(int j, int k, XType type)
    {        
    }
    
    
    BoundaryCondition BCLeft, BoundaryCondition BCRight>
    PrimalBasis ULeft, URight;
};

template <QuadratureType Quad>
struct MyIntegral
{
    double
    operator()(int j1, int k1, XType type1, int deriv1,
               int j2, int k2, XType type2, int deriv2)
    {        
    }

    const First  &first;
    const Second &second;
};



int
main()
{
    
}