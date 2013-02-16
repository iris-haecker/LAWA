/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef CXXBLAS_LEVEL2_TRMV_TCC
#define CXXBLAS_LEVEL2_TRMV_TCC 1

namespace cxxblas {

#ifdef HAVE_CBLAS

// strmv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trmv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const float *A, IndexType ldA,
     float *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_strmv");

    cblas_strmv(CBLAS::value(order), CBLAS::value(upLo),
                CBLAS::value(transA), CBLAS::value(diag),
                n,
                A, ldA,
                x, incX);
}

// dtrmv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trmv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const double *A, IndexType ldA,
     double *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_dtrmv");

    cblas_dtrmv(CBLAS::value(order), CBLAS::value(upLo),
                CBLAS::value(transA), CBLAS::value(diag),
                n,
                A, ldA,
                x, incX);
}

// ctrmv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trmv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const ComplexFloat *A, IndexType ldA,
     ComplexFloat *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_ctrmv");

    cblas_ctrmv(CBLAS::value(order), CBLAS::value(upLo),
                CBLAS::value(transA), CBLAS::value(diag),
                n,
                reinterpret_cast<const float *>(A), ldA,
                reinterpret_cast<const float *>(x), incX);
}

// ztrmv
template <typename IndexType>
typename If<IndexType>::isBlasCompatibleInteger
trmv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const ComplexDouble *A, IndexType ldA,
     ComplexDouble *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("[" BLAS_IMPL "] cblas_ztrmv");

    cblas_ztrmv(CBLAS::value(order), CBLAS::value(upLo),
                CBLAS::value(transA), CBLAS::value(diag),
                n,
                reinterpret_cast<const double *>(A), ldA,
                reinterpret_cast<const double *>(x), incX);
}

#endif // HAVE_CBLAS

} // namespace flens

#endif // CXXBLAS_LEVEL2_TRMV_TCC

