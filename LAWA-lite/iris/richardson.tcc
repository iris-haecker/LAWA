#ifndef IRIS_RICHARDSON_TCC
#define IRIS_RICHARDSON_TCC 1

template <typename MatrixA, typename VectorB, typename Omega, typename VectorX,
   		  typename Tol> 
int
richardson(const MatrixA &A, const VectorB &b, Omega w, VectorX &x, Tol tol,
  		   int maxIt)
{
	VectorX r = b - A*x;
	
	for (int k=1; k<=maxIt; ++k) {
		x += w*r;
		
		typename VectorX::ElementType rNorm = sqrt(r*r);
		
		std::cerr << "rNorm = " << rNorm << std::endl;
		
		if (rNorm<tol) {
			return k;
		}
		r = b - A*x;
	}
	return maxIt;
}

#endif // IRIS_DAMPEDRICHARDSON_TCC