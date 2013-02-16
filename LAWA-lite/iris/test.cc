#include <iostream>
#include <iris/iris.cxx>
#include <vector>

using namespace lawa;
using namespace std;


int
main()
{
    typedef flens::DenseVector<Array<double> >      RealDenseVector;
    
    vector<RealDenseVector>  vec(5);

    vec[2].engine().resize(5);

    for (size_t i=0; i<vec.size(); ++i) {
        cout << "vec(" << i << ") =" << vec[i] << std::endl;
    }

}