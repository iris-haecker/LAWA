#include <iostream>
#define SOLVER_DEBUG

#include <iris/iris.cxx>


using namespace lawa;
using namespace std;

int
main()
{
    typedef flens::DenseVector<Array<double> >               RealDenseVector;
    typedef flens::DenseVector<Array<int> >                  IntDenseVector;


    RealDenseVector   x(5);
    IntDenseVector    xi;

    x = 1, -2, 3, 1.5, -5;

    myAbsSort(x, xi);

    cout << "x  = " << x << endl;
    cout << "xi = " << xi << endl;

    for (int k=1; k<=xi.length(); ++k) {
        cout << x(xi(k)) << "  ";
    }
    cout << endl;

}