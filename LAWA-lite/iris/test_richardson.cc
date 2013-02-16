#include <iostream>

#include <iris/iris.cxx>
#include <iris/richardson.h>
#include <iris/richardson.tcc>


using namespace lawa;
using namespace std;

int
main()
{
    typedef flens::DenseVector<Array<double> >      		 RealDenseVector;
    typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;

	RealGeMatrix  A(4,4);
	
	A.diag( 0) =  2;
	A.diag( 1) = -1;
	A.diag(-1) = -1;

	RealDenseVector  b(4);
	b = 1, 2, 3, 4;

	RealDenseVector  x(4);
	x = 0;

	double w   = 0.4;
	double tol = 0.000000001;
	int maxIt  = 1000;

	int it = richardson(A, b, w, x, tol, maxIt);
	
	cout << "it = " << it << endl;


}