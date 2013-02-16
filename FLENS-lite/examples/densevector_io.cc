#include <flens/flens.h>

using namespace flens;
using namespace std;

int
main()
{
    DenseVector<Array<double> > x(4), y(4), z;
    x = 1, 2, 3, 4;
    y = 4, 5, 6, 7;
    
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    z = x+y;
    cout << "z = " << z << endl;

    GeMatrix<FullStorage<double> > A(4,4), B(4,4), C;
    A = 1, 2, 3, 4,
        5, 6, 7, 9,
        9, 5, 3, 1,
        1, 5, 2, 9;

    B = 5, 2, 3, 4,
        5, 6, 7, 9,
        9, 5, 3, 1,
        1, 5, 2, 9;

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    C = A+B;
    cout << "C = " << C << endl;

    C = transpose(B);
    cout << "C = " << C << endl;
    C = transpose(transpose(B));
    cout << "C = " << C << endl;

    C += 2*transpose(B);
    cout << "C = " << C << endl;

    C += transpose(2*B);
    cout << "C = " << C << endl;

    C = A+transpose(B);
    cout << "C = " << C << endl;

    z = A*x;
    cout << "z = " << z << endl;

    z = A*x+y;
    cout << "z = " << z << endl;

    z = y+A*x;
    cout << "z = " << z << endl;

    z = transpose(2*A)*x+y;
    cout << "z = " << z << endl;

    z = y+2*transpose(A)*x;
    cout << "z = " << z << endl;
}