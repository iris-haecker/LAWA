#ifndef IRIS_MYSOLVE_H
#define IRIS_MYSOLVE_H 1

#include <map>

namespace lawa {
	
	
	

template <typename VW, typename VLAMBDA, typename T>
void
MYGROW(const DenseVector<VW> &w, T nu_bar, T &nu, std::map<int, int> Lambda)
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
		for (std::map<int, T>::const_iterator it=Lambda.begin(); it!=Lambda.end(); ++it) {
			normRLambdaSquare += pow(r(it->first),2);
		}
		while (sqrt(normRLambdaSquare)<alpha*rNorm) {
			
		}
		
		
	}
}

} // namespace lawa

#endif // IRIS_MYSOLVE_H