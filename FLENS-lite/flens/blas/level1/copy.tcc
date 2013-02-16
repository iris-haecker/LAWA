/*
 *   Copyright (c) 2009, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FLENS_BLAS_LEVEL1_COPY_TCC
#define FLENS_BLAS_LEVEL1_COPY_TCC 1

#include <flens/aux/macros.h>
#include <flens/blas/debugmacro.h>
#include <flens/storage/storageinfo.h>

namespace flens { namespace blas {

//-- common interface for vectors ----------------------------------------------
template <typename VX, typename VY>
void
copy(const Vector<VX> &x, Vector<VY> &y)
{
    copy(x.impl(), y.impl());
}

//-- common interface for matrices ---------------------------------------------
template <typename MA, typename MB>
void
copy(cxxblas::Transpose trans, const Matrix<MA> &A, Matrix<MB> &B)
{
    copy(trans, A.impl(), B.impl());
}

//-- copy
template <typename VX, typename VY>
void
copy(const DenseVector<VX> &x, DenseVector<VY> &y)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(x, y);
    if (y.length()!=x.length()) {
        y.engine().resize(x.engine());
    }
    y.engine().changeIndexBase(x.firstIndex());

#   ifdef HAVE_CXXBLAS_COPY
    cxxblas::copy(x.length(),
                  x.engine().data(), x.engine().stride(),
                  y.engine().data(), y.engine().stride());
#   else
    ASSERT(0);
#   endif

    FLENS_CLOSURELOG_END_ENTRY
}

//-- gecopy
template <typename MA, typename MB>
void
copy(cxxblas::Transpose trans,
     const GeMatrix<MA> &A, GeMatrix<MB> &B)
{
    FLENS_CLOSURELOG_ADD_ENTRY_COPY(A, B);

    if (trans==cxxblas::NoTrans) {
        if ((A.numRows()!=B.numRows())
         || (A.numCols()!=B.numCols())) {
            B.engine().resize(A.engine());
        }
    } else {
        if ((A.numRows()!=B.numCols())
         || (A.numCols()!=B.numRows())) {
            B.engine().resize(A.numCols(),A.numRows(), A.firstCol(), A.firstRow());
        } 
    }
    trans = (StorageInfo<MA>::Order==StorageInfo<MB>::Order)
          ? cxxblas::Transpose(trans ^ cxxblas::NoTrans)
          : cxxblas::Transpose(trans ^ cxxblas::Trans);

#   ifdef HAVE_CXXBLAS_GECOPY
    cxxblas::gecopy(StorageInfo<MB>::Order, trans,
                    B.numRows(), B.numCols(),
                    A.engine().data(), A.engine().leadingDimension(),
                    B.engine().data(), B.engine().leadingDimension());
#   else
    ASSERT(0);
#   endif
    FLENS_CLOSURELOG_END_ENTRY
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_COPY_TCC
