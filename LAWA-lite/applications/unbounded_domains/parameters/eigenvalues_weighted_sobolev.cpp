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
typedef flens::DenseVector<flens::Array<T> >                 DenseVectorT;

//Operator definitions
typedef PDENonConstCoeffOperator1D<T, Basis1D>               PDEOperator1D;
typedef NoCompression<T, Index1D, Basis1D>                   Compression;
typedef WeightedSobolevNormPreconditioner1D<T,Basis1D>       NormPreconditioner1D;
typedef WeightedSobolevMidPointPreconditioner1D<T,Basis1D>   MidPointPreconditioner1D;

//MapMatrix definition
typedef MapMatrix<T,Index1D,PDEOperator1D,
                  Compression,NormPreconditioner1D>          MA_norm_prec;
typedef MapMatrix<T,Index1D,PDEOperator1D,
                  Compression,MidPointPreconditioner1D>      MA_midpoint_prec;


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

/// Non-constant weight function
T
weight_f(T x)
{
    return std::exp(-2.*std::fabs(x));
}

/// Zero convection term
T
b_f(T /*x*/)
{
    return 0.;
}



int main (int argc, char *argv[]) {
    if (argc != 6) {
        cout << "usage " << argv[0] << " d d_ radius jmin max_level" << endl;
        exit(1);
    }

    int d=atoi(argv[1]), d_=atoi(argv[2]);
    T   radius    =atof(argv[3]);
    int jmin =atoi(argv[4]);
    int max_level =atoi(argv[5]);
    cout.precision(8);

    Basis1D                  basis(d,d_,jmin);
    DenseVectorT             weight_singPts(1); weight_singPts = 0.;
    DenseVectorT             b_singPts;
    Function<T>              weight(weight_f, weight_singPts);
    Function<T>              b(b_f, b_singPts);
    PDEOperator1D            Bil(basis, weight, b, weight, 2*basis.d);
    NormPreconditioner1D     NormP(basis,weight,1);
    MidPointPreconditioner1D MidPointP(basis,weight,1);
    Compression              compression(basis);
    MA_norm_prec             A_norm_prec(Bil,NormP,compression);
    MA_midpoint_prec         A_midpoint_prec(Bil,MidPointP,compression);

    cout << "Chosen parameters: d=" << d << ", d_=" << d_  << ", radius="<< radius
         << ", jmin=" << jmin << ", max_level=" << max_level << endl;



    T cB_norm, CB_norm, cB_midpoint, CB_midpoint;
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
    std::stringstream filename;
    filename << "eigenvalues_weighted_sobolev_" << jmin << "_" << d << "_" << d_  << ".dat";
    std::ofstream file(filename.str().c_str());
    for (int jmax=0; jmax<=max_level; jmax+=1) {
        for (T r=1.; r<=radius; r+=1.) {
            IndexSet<Index1D> Lambda = LambdaForEigenvalues(jmin, jmax, basis, r, true);
            int N = Lambda.size();
            cout << "Size of Lambda: " << N << endl;

            cout << "Norm preconditioning ..." << endl;
            flens::SparseGeMatrix<CRS<T,CRS_General> > A_norm_prec_sparse(N,N);
            toFlensSparseMatrix(A_norm_prec, Lambda, Lambda, A_norm_prec_sparse);
            densify(cxxblas::NoTrans,A_norm_prec_sparse,A_dense);
            cB_norm=0., CB_norm=0.;
            computeEV(A_dense, cB_norm, CB_norm);
            cout             << " " << jmax << " " << r
                             << " " << cB_norm << " " << " " << CB_norm << endl;


            cout << "Midpoint preconditioning..." << endl;
            flens::SparseGeMatrix<CRS<T,CRS_General> > A_midpoint_prec_sparse(N,N);
            toFlensSparseMatrix(A_midpoint_prec, Lambda, Lambda, A_midpoint_prec_sparse);
            densify(cxxblas::NoTrans,A_midpoint_prec_sparse,A_dense);
            cB_midpoint=0., CB_midpoint=0.;
            computeEV(A_dense, cB_midpoint, CB_midpoint);
            cout             << " " << jmax << " " << r
                             << " " << cB_midpoint << " " << " " << CB_midpoint << endl << endl;

            file << r << " " << jmax-jmin << " " << cB_norm << " " << CB_norm << " "
                 << cB_midpoint << " " << CB_midpoint << endl;
        }
        file << endl;
    }
    file.close();
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
