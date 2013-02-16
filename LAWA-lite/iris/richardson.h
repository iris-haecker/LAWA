#ifndef IRIS_RICHARDSON_H
#define IRIS_RICHARDSON_H 1

template <typename MatrixA, typename VectorB, typename Omega,
 		  typename VectorX, typename Tol> 
	int
	richardson(const MatrixA &A, const VectorB &b, Omega w, VectorX &x,
	 		   Tol tol, int maxIt);

#endif // IRIS_DAMPEDRICHARDSON_H