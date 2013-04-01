#ifndef IRIS_MYEVAL_TCC
#define IRIS_MYEVAL_TCC 1

#include <lawa/flensforlawa.h>
#include <iris/mybasis.h>
#include <fstream>

namespace lawa {

template <typename T>
MyEval<T>::MyEval(const MyBasis<T> &_U, const CoeffVector &_u)
    : U(_U), u(_u)
{
    assert(u.firstIndex()==1);
}

template <typename T>
T
MyEval<T>::operator()(T x) const
{
    T y = 0;
    for (int p=1; p<=u.length(); ++p) {
        y += u(p)*U(x, p, 0);
    }
    return y;
}

template <typename T>
void
MyEval<T>::dump(int N, const char *file) const
{
    std::fstream   out(file, std::fstream::out);

    for (int i=0; i<=N; ++i) {
        const double x = double(i)/N;
        out << x << " " << operator()(x) << std::endl;
    }
}

template <typename T>
void
MyEval<T>::dump(int N, const Function<T> &refSol, const char *file) const
{
    std::fstream   out(file, std::fstream::out);

    for (int i=0; i<=N; ++i) {
        const double x = double(i)/N;
        out << x << " " << operator()(x) << " " << refSol(x) << std::endl;
    }
}

template <typename T>
T
MyEval<T>::diff_L1(int N, const Function<T> &f) const
{
    using std::abs;

    double diff1;
    double diff2;

    for (int k=1; k<=10; ++k) {

        diff1 = 0;
        diff2 = 0;

        for (int i=0; i<=2*N; ++i) {
            const double x     = double(i)/(2*N);
            const double value = abs(f(x) - operator()(x));
            
            if (i==0 || i==2*N) {
                diff1 += 0.5*value;
                diff2 += 0.5*value;
                continue;
            }

            if (i%2 == 0) {
                diff1 += value;
            }
            diff2 += value;
        }

        diff1 /= N;
        diff2 /= 2*N;
        
        if (abs(diff1-diff2)<0.00001) {
            break;
        }
    }
    return diff2;
}

template <typename T>
T
MyEval<T>::diff_L2(int N, const Function<T> &f) const
{
    using std::pow;

    double diff1;
    double diff2;

    for (int k=1; k<=10; ++k) {

        diff1 = 0;
        diff2 = 0;

        for (int i=0; i<=2*N; ++i) {
            const double x     = double(i)/(2*N);
            const double value = pow(f(x) - operator()(x), 2);
            
            if (i==0 || i==2*N) {
                diff1 += 0.5*value;
                diff2 += 0.5*value;
                continue;
            }

            if (i%2 == 0) {
                diff1 += value;
            }
            diff2 += value;
        }
        diff1 /= N;
        diff2 /= 2*N;

        if (abs(diff1-diff2)<0.00001) {
            break;
        }
    }
    return diff2;
}

template <typename T>
T
MyEval<T>::diff_LInf(int N, const Function<T> &f) const
{
    using std::abs;
    using std::max;

    double diff1;
    double diff2;

    for (int k=1; k<=10; ++k) {

        diff1 = 0;
        diff2 = 0;

        for (int i=0; i<=2*N; ++i) {
            const double x     = double(i)/(2*N);
            const double value = abs(f(x) - operator()(x));
            
            if (i%2 == 0) {
                diff1 = max(diff1, value);
            }
            diff2 = max(diff2, value);
        }

        if (abs(diff1-diff2)<0.00001) {
            break;
        }
    }
    return diff2;
}

} // namespace lawa

#endif // IRIS_MYEVAL_TCC
