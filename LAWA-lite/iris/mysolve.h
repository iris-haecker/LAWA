#ifndef IRIS_MYSOLVE_H
#define IRIS_MYSOLVE_H 1

#include <iris/mybasis.h>

namespace lawa {

template <typename VW, typename VLAMBDA, typename T>
void
MYGROW(const DenseVector<VW> &w, T nu_bar, T &nu, DenseVector<VLAMBDA> &Lambda, int &N)
{
	zeta = 2*omega*nu_bar/(1-omega);
	do {
		zeta = zeta/2;
		r = MYRHS(zeta/2);
		r -= MYAPPLY(w,zeta/2);
		
		if
		
	} while (true);
}

} // namespace lawa

#endif // IRIS_MYSOLVE_H