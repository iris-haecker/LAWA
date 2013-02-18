#ifndef IRIS_MYSOLVE_TCC
#define IRIS_MYSOLVE_TCC 1

#include <map>

namespace lawa {

template <typename VX, typename TI>
void
support(const DenseVector<VX> &x, IndexSet<TI> &Lambda)
{
    using std::abs;

    typedef typename DenseVector<VX>::ElementType ElementType;
    typedef typename DenseVector<VX>::IndexType   IndexType;

    assert(x.firstIndex()==IndexType(1));

    const ElementType Zero(0);

    for (IndexType k=1; k<=x.length(); ++k) {
        if (abs(x(k))>Zero) {
            Lambda.insert(k);
        }
    }
}

template <typename VW, typename VLAMBDA, typename T>
void
MYGROW(const DenseVector<VW> &w, T nu_bar, T &nu, IndexSet<int> &Lambda)
{
    zeta = 2*omega*nu_bar/(1-omega);
    do {
        zeta = zeta/2;
        r = MYRHS(zeta/2);
        r -= MYAPPLY(w,zeta/2);

        rNorm = sqrt(r*r) + zeta;
        if (rNorm + zeta <= epsilon) {
            break;
        }
        if (zeta <= omega*rNorm) {
            break;
        }
        
    } while (true);

    if (nu>epsilon) {
        support(w, Lambda);

        T normRLambdaSquare = 0;

        IndexSet<int>::const_iterator it;
        for (it=Lambda.begin(); it!=Lambda.end(); ++it) {
            normRLambdaSquare += pow(r(k), 2);
        }

        if (normRLambdaSquare<pow(alpha*rNorm, 2)) {

            myAbsSort(r, rAbsSorted);

            for (int k=1; k<=rAbsSorted.length(); ++k) {
                if (Lambda.count(rAbsSorted(k))>0) {
                    continue;
                }
                normRLambdaSquare += pow(r(rAbsSorted(k)), 2);
                if (normRLambdaSquare>=pow(alpha*rNorm,2)) {
                    break;
                }
            }

        }

    }
}

} // namespace lawa

#endif // IRIS_MYSOLVE_H