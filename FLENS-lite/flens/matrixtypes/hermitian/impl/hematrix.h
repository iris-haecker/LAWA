/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_MATRIXTYPES_HERMITIAN_IMPL_HEMATRIX_H
#define FLENS_MATRIXTYPES_HERMITIAN_IMPL_HEMATRIX_H 1

#include <cxxblas/typedefs.h>
#include <flens/matrixtypes/hermitian/hermitianmatrix.h>

namespace flens {

// forward declarations
template <typename FS>
    class GeMatrix;

template <typename FS>
    class SyMatrix;

template <typename FS>
    class TrMatrix;

template <typename FS>
class HeMatrix
    : public HermitianMatrix<HeMatrix<FS> >
{
    public:
        typedef FS                                  Engine;
        typedef typename Engine::ElementType        ElementType;
        typedef typename Engine::IndexType          IndexType;

        // view types from Engine
        typedef typename Engine::ConstView          EngineConstView;
        typedef typename Engine::View               EngineView;
        typedef typename Engine::NoView             EngineNoView;

        typedef typename Engine::ConstArrayView     ConstArrayView;
        typedef typename Engine::ArrayView          ArrayView;
        typedef typename Engine::Array              Array;

        // view types
        typedef DenseVector<ConstArrayView>         ConstVectorView;
        typedef DenseVector<ArrayView>              VectorView;
        typedef DenseVector<Array>                  Vector;

        typedef GeMatrix<EngineConstView>           ConstGeneralView;
        typedef GeMatrix<EngineView>                GeneralView;
        typedef GeMatrix<EngineNoView>              GeneralNoView;

        typedef HeMatrix<EngineConstView>           ConstView;
        typedef HeMatrix<EngineView>                View;
        typedef HeMatrix<EngineNoView>              NoView;

        typedef SyMatrix<EngineConstView>           ConstSymmetricView;
        typedef SyMatrix<EngineView>                SymmetricView;
        typedef SyMatrix<EngineNoView>              SymmetricNoView;

        typedef TrMatrix<EngineConstView>           ConstTriangularView;
        typedef TrMatrix<EngineView>                TriangularView;
        typedef TrMatrix<EngineNoView>              TriangularNoView;

        HeMatrix();

        explicit
        HeMatrix(IndexType dim);

        HeMatrix(IndexType dim,
                 IndexType firstRow, IndexType firstCol);

        HeMatrix(const Engine &engine, cxxblas::StorageUpLo upLo);

        HeMatrix(const HeMatrix &rhs);

        template <typename RHS>
            HeMatrix(const HeMatrix<RHS> &rhs);

        template <typename RHS>
            HeMatrix(HeMatrix<RHS> &rhs);

        template <typename RHS>
            HeMatrix(const Matrix<RHS> &rhs);

        // -- operators --------------------------------------------------------

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        // -- views ------------------------------------------------------------

        // general views
        ConstGeneralView
        general() const;

        GeneralView
        general();

        // -- methods ----------------------------------------------------------

        IndexType
        dim() const;

        IndexType
        firstRow() const;

        IndexType
        lastRow() const;

        IndexType
        firstCol() const;

        IndexType
        lastCol() const;

        // -- implementation ---------------------------------------------------

        const Engine &
        engine() const;

        Engine &
        engine();

        cxxblas::StorageUpLo
        upLo() const;

        cxxblas::StorageUpLo &
        upLo();

    private:
        Engine               _engine;
        cxxblas::StorageUpLo _upLo;
};

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HEMATRIX_H
