#include <cmath>
#include <fstream>
#include <iostream>

#include <lawa/lawa.h>

#include <validate/mybasis.h>
#include <validate/mybasis.tcc>

#include <validate/myintegral.h>
#include <validate/myintegral.tcc>


using namespace lawa;
using namespace std;

typedef flens::GeMatrix<FullStorage<double, ColMajor> >  RealGeMatrix;

int
main()
{
    MyBasis<double> U(2, 4);
    
    const int j0 = U.j0();
    
    cout << U.length(j0-1) << endl;
    cout << U.length(j0) << endl;
    cout << U.length(j0+1) << endl;
    
    for (int j=j0-1; j<=j0+2; ++j) {
        cerr << "j = " << U.getActualLevel(j)
             << ", type = " << U.getXType(j)
             << " :  " << U.firstIndex(j)
             << "... " << U.lastIndex(j)
             << ", length = " << U.length(j)
             << endl;
    }

    cerr << endl;


    U.basisLeft.enforceBoundaryCondition<DirichletBC>();
    
    
    //U.basisRight.enforceBoundaryCondition<DirichletBC>();

    for (int j=j0-1; j<=j0+2; ++j) {
        cerr << "j = " << U.getActualLevel(j)
             << ", type = " << U.getXType(j)
             << " :  " << U.firstIndex(j)
             << "... " << U.lastIndex(j)
             << ", length = " << U.length(j)
             << endl;
    }

//
//  Plot functions on level J
//
    const int N = 500;
    
    for (int j=j0-1; j<=j0+2; ++j) {
        for (int k=U.firstIndex(j); k<=U.lastIndex(j); ++k) {
            for (int i=0; i<=N; ++i) {
                const double x = double(i)/N;
                cout << x << " " << U(x, j, k) << endl;
            }
            cout << endl << endl;
        }
    }

    int count = 0;
    for (int j=j0-1; j<=j0+2; ++j) {
        for (int k=U.firstIndex(j); k<=U.lastIndex(j); ++k) {
            ++count;
            cerr << "U.getLevel(" << count << ") = " << U.getLevel(count)
                 << " should be " << j << endl
                 << "U.getIndex(" << count << ") = " << U.getIndex(count)
                 << " should be " << k << endl
                 << endl;
        }
    }

    MyIntegral<double>  myIntegral(U, U);
    
    const int J = j0+1;
    
    Range<int> K1(U.getFirstAbsoluteIndex(j0-1), U.getLastAbsoluteIndex(J));
    Range<int> K2(U.getFirstAbsoluteIndex(j0-1), U.getLastAbsoluteIndex(J));
    
    RealGeMatrix  A(K1, K2);
    
    for (int k1=K1.firstIndex(); k1<=K1.lastIndex(); ++k1) {
        for (int k2=K2.firstIndex(); k2<=K2.lastIndex(); ++k2) {
            A(k1, k2) = myIntegral(k1, 0, k2, 0);
        }
    }

    cerr << "A = " << A << endl;

}
