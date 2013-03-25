/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <cassert>
#include <limits>
#include <iomanip>

namespace lawa {

template <typename MA, typename VX, typename VB>
int
cg(const MA &A, VX &x, const VB &b, typename _cg<VB>::T tol,
   int maxIterations)
{
    typename _cg<VB>::T alpha, beta, rNormSquare, rNormSquarePrev;
    typename _cg<VB>::AuxVector Ap, r, p;

    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
        std::cerr << "x.length() = " << x.length() << std::endl;
    }

    r = A*x - b;
    p = -1*r;
    rNormSquare = r*r;

    for (int k=1; k<=maxIterations; k++) {
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << sqrt(rNormSquare)
                << std::endl;
        #endif
        if (sqrt(rNormSquare)<=tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rNormSquare/(p * Ap);
        x += alpha*p;
        r += alpha*Ap;

        rNormSquarePrev = rNormSquare;
        rNormSquare = r*r;
        beta = rNormSquare/rNormSquarePrev;
        p = beta*p - r;

        /*
        std::cout << "k = " << k
                  << ", rNormSquare = " << rNormSquare
                  << ", rho = " << sqrt(rNormSquare)
                  << std::endl;
        */
    }
    return maxIterations;
}

template <typename MA, typename VX, typename VB>
int
cgls(const MA &A, VX &x, const VB &b, typename _cg<VB>::T tol,
     int maxIterations)
{
    typename _cg<VB>::T alpha, beta, gammaPrev, gamma, b_norm;
    typename _cg<VB>::AuxVector r, q, s, p;

    assert(b.length()==A.numRows());
/*
    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
    }
    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        x(i) = 0;
    }
*/
    b_norm = b*b;
    if (std::sqrt(b_norm) < 1e-15) {
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            x(i) = 0;
        }
        return 0;
    }

    //r = b;
    r = b-A*x;
    flens::blas::mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);
    p = s;
    gammaPrev = s*s;
    for (int k=1; k<=maxIterations; k++) {
        q = A*p;
        alpha = gammaPrev/(q*q);
        x +=   alpha *p;
        r += (-alpha)*q;
        flens::blas::mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);
        gamma = s*s;
        if (sqrt(gamma)<=tol) {
            std::cerr << "    cgls: gamma = " << gamma << std::endl;
            return k-1;
        }
        beta  = gamma/gammaPrev;
        p *= beta;
        p += s;
        gammaPrev = gamma;
    }
    return maxIterations;
}


// Algorithm 9.2, Y. Saad: Iterative Methods for Sparse Linear Systems
// for solving Ax=b with P^T A P u = P^T b, u=P^{-1} x
// Note the role of P and P^T is switched.
// This algorithm uses split preconditioning defined as in
// (6.2), K. Urban: Wavelet Methods for Elliptic PDEs
// Also non-symmetric preconditioning is possible as mentioned in remark 6.2
template <typename Prec, typename MA, typename VX, typename VB>
int
pcg(const Prec &P, const MA &A, VX &x, const VB &b,
    typename _cg<VB>::T tol, int maxIterations)
{
    typename _cg<VB>::T pNormSquare, alpha, beta, rHatq, rHatqPrev;
    typename _cg<VB>::AuxVector r, rHat, p, Ap;

    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
    }
    
    r.engine().resize(b.length(), 1);
    rHat.engine().resize(r.length(), 1);
    p.engine().resize(r.length(), 1);
    Ap.engine().resize(A.numRows(), 1);
    
    r = b - A*x;
    if(r*r == 0){
        return 0;
    }    
    rHat = transpose(P)*r;
    p = P*rHat;
    rHatq = rHat*rHat;
    for (int k=1; k<=maxIterations; k++) {
        Ap = A*p;
        alpha = rHatq/(Ap*p);
        x += alpha*p;
        rHat = rHat - alpha * transpose(P)*Ap;
        rHatqPrev = rHatq;
        rHatq =rHat*rHat;
        beta = rHatq/rHatqPrev;
        p *= beta;
        p += P*rHat;
        pNormSquare = p*p;
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " <<  sqrt(pNormSquare) 
                << std::endl;
        #endif

        /*
        std::cout << "k = " << k
                  << ", pNormSquare = " << pNormSquare
                  << ", rho = " << sqrt(pNormSquare)
                  << std::endl;
        */

        if (sqrt(pNormSquare)<=tol) {
            return k-1;
        }
    }
    return maxIterations;
}

/*
template <typename Prec, typename MA, typename VX, typename VB>
int
pcg(const Prec &B, const MA &A, VX &x, const VB &b,
    typename _cg<VB>::T tol, int maxIterations)
{
    typename _cg<VB>::T pNormSquare, alpha, beta, rq, rqPrev;
    typename _cg<VB>::AuxVector r, q, p, Ap;

    r = A*x - b;
    q = B*r;

    p = q;
    // TODO: next line results in an error with T = long double. WHY???
    // p = -q;
    p *= -1;
    rq = r*q;

    for (int k=1; k<=maxIterations; k++) {
        pNormSquare = p*p;
        if (sqrt(pNormSquare)<tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rq/(p*Ap);
        x += alpha*p;

        r += alpha*Ap;
        q = B*r;

        rqPrev = rq;
        rq = r*q;
        beta = rq/rqPrev;
        p = beta*p - q;
    }
    return maxIterations;
}
*/
} // namespace lawa

