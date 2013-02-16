/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#include <iostream>
#include <vector>
#include <lawa/lawa.h>

typedef double T;
using namespace lawa;
using namespace std;

typedef Basis<T,Primal,R,CDF> Basis1D;
typedef Wavelet<T,Primal,R,CDF> WaveletR;

//Operator definitions
typedef HelmholtzOperator1D<T, Basis1D>          HelmholtzBilinearForm1D;
typedef CompressionPDE1D<T, Basis1D>             Compression_W_XBSpline1D;
typedef CompressionPDE1D_WO_XBSpline<T, Basis1D> Compression_WO_XBSpline1D;
typedef H1NormPreconditioner1D<T,Basis1D>        Preconditioner1D;

//MapMatrix definition
typedef MapMatrix<T,Index1D,HelmholtzBilinearForm1D,
                  Compression_W_XBSpline1D,Preconditioner1D>   MA_W_XBSpline;
typedef MapMatrix<T,Index1D,HelmholtzBilinearForm1D,
                  Compression_WO_XBSpline1D,Preconditioner1D>  MA_WO_XBSpline;


template <typename T>
void
computeEV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB);

template <typename T>
void
computeSV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB);

template<typename T>
IndexSet<Index1D>
LambdaForEigenvalues(int jmin, int jmax, const Basis1D &basis, T radius, bool with_BSplines);

typedef IndexSet<Index1D>::const_iterator const_set_it;

int main (int argc, char *argv[]) {
    if (argc != 8) {
        cout << "usage " << argv[0] << " d d_ c radius jmin max_level W_XBSpline" << endl;
        exit(1);
    }

    int d=atoi(argv[1]), d_=atoi(argv[2]);
    T   c=atof(argv[3]);
    T   radius    =atof(argv[4]);
    int jmin =atoi(argv[5]);
    int max_level =atoi(argv[6]);
    int W_XBSpline=atoi(argv[7]);
    cout.precision(8);

    Basis1D basis(d,d_,jmin);
    HelmholtzBilinearForm1D Bil(basis,c);
    Preconditioner1D P(basis);
    Compression_W_XBSpline1D  compression_w_XBSpline(basis);
    Compression_WO_XBSpline1D compression_wo_XBSpline(basis);
    MA_W_XBSpline  A_w_XBSpline(Bil,P,compression_w_XBSpline);
    MA_WO_XBSpline A_wo_XBSpline(Bil,P,compression_wo_XBSpline);

    cout << "Chosen parameters: d=" << d << ", d_=" << d_ << ", c=" << c << ", radius="<< radius
         << ", jmin=" << jmin << ", max_level=" << max_level << endl;


    if (W_XBSpline==1) {
        std::stringstream filename;
        filename << "eigenvalues_helmholtz_W_XBSpline_" << jmin << "_" << d << "_" << d_  << ".dat";
        std::ofstream file_eigenvalues(filename.str().c_str());
        for (int jmax=0; jmax<=max_level; jmax+=1) {
            for (T r=1.; r<=radius; r+=1.) {
                IndexSet<Index1D> Lambda = LambdaForEigenvalues(jmin, jmax, basis, r, true);
                int N = Lambda.size();
                cout << "Size of Lambda: " << N << endl;
                flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
                toFlensSparseMatrix(A_w_XBSpline, Lambda, Lambda, A_flens);
                flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
                densify(cxxblas::NoTrans,A_flens,A_dense);

                T cB, CB;
                computeEV(A_dense, cB, CB);
                file_eigenvalues << " " << jmax << " " << r
                                 << " " << cB << " " << " " << CB << endl;
                cout             << " " << jmax << " " << r
                                 << " " << cB << " " << " " << CB << endl;
            }
        }
    }
    else {
        std::stringstream filename;
        filename << "eigenvalues_helmholtz_WO_XBSpline_" << jmin << "_" << d << "_" << d_  << ".dat";
        std::ofstream file_eigenvalues(filename.str().c_str());
        for (int jmax=0; jmax<=max_level; jmax+=1) {
            for (T r=1.; r<=radius; r+=1.) {
                IndexSet<Index1D> Lambda = LambdaForEigenvalues(jmin, jmax, basis, r, false);
                int N = Lambda.size();
                cout << "Size of Lambda: " << N << endl;
                flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
                toFlensSparseMatrix(A_wo_XBSpline, Lambda, Lambda, A_flens);
                flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
                densify(cxxblas::NoTrans,A_flens,A_dense);

                T cB, CB;
                computeEV(A_dense, cB, CB);
                file_eigenvalues << " " << jmax << " " << r
                                 << " " << cB << " " << " " << CB << endl;
                cout             << " " << jmax << " " << r
                                 << " " << cB << " " << " " << CB << endl;

            }
        }
    }


    return 0;
}

template <typename T>
void
computeEV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB) {
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U(A.numRows(),A.numRows()),
                                                               V(A.numCols(),A.numCols());
    int N = A.numRows();
    DenseVector<Array<T> > wr(N), wi(N);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > vl,vr;
    ev(false, false, A, wr, wi, vl, vr);
    cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
}

template <typename T>
void
computeSV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB) {
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U(A.numRows(),A.numRows()),
                                                               V(A.numCols(),A.numCols());
    DenseVector<Array<T> > s(A.numCols());
    int iterations = svd(A,s,U,V);
    CB = s(s.firstIndex());
    cB = s(s.lastIndex());
}

template<typename T>
IndexSet<Index1D>
LambdaForEigenvalues(int jmin, int jmax, const Basis1D &basis, T radius, bool with_BSplines)
{
    const BSpline<T,Primal,R,CDF> phi = basis.mra.phi;
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;
    IndexSet<Index1D> Lambda;
    int k_left, k_right;
    if (with_BSplines) {
        for (int j=jmin; j<=jmax; ++j) {
            k_left = std::floor(-pow2i<T>(j)*radius-psi.support(0,0).l2);
            k_right = std::ceil(pow2i<T>(j)*radius-psi.support(0,0).l1);
            for (int k=k_left; k<=k_right; ++k) {
                Lambda.insert(Index1D(j,k,XWavelet));
            }
        }

        k_left  = int(std::floor(-pow2i<T>(jmin)*radius-phi.support(0,0).l2));
        k_right = int(std::ceil(  pow2i<T>(jmin)*radius-phi.support(0,0).l1));
        for (int k=k_left; k<=k_right; ++k) {
            Lambda.insert(Index1D(jmin,k,XBSpline));
        }
    }
    else {
        for (int j=jmin; j<=jmax; ++j) {
                int k_left, k_right;
                if (j>=-6) {
                    k_left = std::floor(-pow2i<T>(j)*radius-psi.support(0,0).l2);
                    k_right = std::ceil(pow2i<T>(j)*radius-psi.support(0,0).l1);
                    if (k_left>=0) cout << "j=" << j << ", k_left=" << k_left << endl;
                    if (k_right<=0) cout << "j=" << j << ", k_right=" << k_right << endl;
                }
                else {
                    k_left = -psi.d-psi.d_;
                    k_right = psi.d+psi.d_;
                }
                for (int k=k_left; k<=k_right; ++k) {
                    Lambda.insert(Index1D(j,k,XWavelet));
                }
            }
    }
    return Lambda;
}
