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

#ifndef FLENS_STORAGE_FULLSTORAGE_FULLSTORAGE_TCC
#define FLENS_STORAGE_FULLSTORAGE_FULLSTORAGE_TCC 1

#include <cxxblas/level1extensions/gecopy.h>

namespace flens {

//= Constructors

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
FullStorage<T, Order, I, A>::FullStorage()
    :  _data(0),
       _numRows(0), _numCols(0),
       _firstRow(I::defaultIndexBase), _firstCol(I::defaultIndexBase)
{
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
FullStorage<T, Order, I, A>::FullStorage(IndexType numRows, IndexType numCols,
                                         IndexType firstRow, IndexType firstCol,
                                         const ElementType &value,
                                         const Allocator &allocator)
    : _data(0), _allocator(allocator),
      _numRows(numRows), _numCols(numCols),
      _firstRow(firstRow), _firstCol(firstCol)
{
    _allocate(value);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
FullStorage<T, Order, I, A>::FullStorage(const FullStorage &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _firstRow(rhs.firstRow()), _firstCol(rhs.firstCol())
{
    _allocate(ElementType());
    cxxblas::gecopy(Order,
                    cxxblas::NoTrans, _numRows, _numCols,
                    rhs.data(), rhs.leadingDimension(),
                    data(), leadingDimension());
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
template <typename RHS>
FullStorage<T, Order, I, A>::FullStorage(const RHS &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _firstRow(rhs.firstRow()), _firstCol(rhs.firstCol())
{
    _allocate(ElementType());
    cxxblas::gecopy(Order,
                    cxxblas::NoTrans, _numRows, _numCols,
                    rhs.data(), rhs.leadingDimension(),
                    data(), leadingDimension());
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
FullStorage<T, Order, I, A>::~FullStorage()
{
    _release();
}

//-- operators -----------------------------------------------------------------

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ElementType &
FullStorage<T, Order, I, A>::operator()(IndexType row, IndexType col) const
{
    ASSERT(row>=_firstRow);
    ASSERT(row<_firstRow+_numRows);
    ASSERT(col>=_firstCol);
    ASSERT(col<_firstCol+_numCols);

    if (Order==cxxblas::ColMajor) {
        return _data[col*_numRows+row];
    }
    return _data[row*_numCols+col];
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::ElementType &
FullStorage<T, Order, I, A>::operator()(IndexType row, IndexType col)
{
    ASSERT(row>=_firstRow);
    ASSERT(row<_firstRow+_numRows);
    ASSERT(col>=_firstCol);
    ASSERT(col<_firstCol+_numCols);

    if (Order==cxxblas::ColMajor) {
        return _data[col*_numRows+row];
    }
    return _data[row*_numCols+col];
}

//== Methods
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::firstRow() const
{
    return _firstRow;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::lastRow() const
{
    return _firstRow+_numRows-1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::firstCol() const
{
    return _firstCol;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::lastCol() const
{
    return _firstCol+_numCols-1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::numRows() const
{
    return _numRows;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::numCols() const
{
    return _numCols;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::leadingDimension() const
{
    return (Order==cxxblas::ColMajor) ? _numRows
                                      : _numCols;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::strideRow() const
{
    return (Order==cxxblas::ColMajor) ? 1
                                      : leadingDimension();
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::IndexType
FullStorage<T, Order, I, A>::strideCol() const
{
    return (Order==cxxblas::ColMajor) ? leadingDimension()
                                      : 1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ElementType *
FullStorage<T, Order, I, A>::data() const
{
    ASSERT(_data);
    return &(this->operator()(_firstRow, _firstCol));
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::ElementType *
FullStorage<T, Order, I, A>::data()
{
    ASSERT(_data);
    return &(this->operator()(_firstRow, _firstCol));
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::Allocator &
FullStorage<T, Order, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
bool
FullStorage<T, Order, I, A>::resize(IndexType numRows, IndexType numCols,
                                    IndexType firstRow, IndexType firstCol,
                                    const ElementType &value)
{
    if ((_numRows!=numRows) || (_numCols!=numCols)) {
        _release();
        _numRows = numRows;
        _numCols = numCols;
        _firstRow = firstRow;
        _firstCol = firstCol;
        _allocate(value);
        return true;
    }
    changeIndexBase(firstRow, firstCol);
    return false;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
bool
FullStorage<T, Order, I, A>::resize(const Range<IndexType> &rows,
                                    const Range<IndexType> &cols,
                                    const ElementType &value)
{
    if ((_numRows!=rows.length()) || (_numCols!=cols.length())) {
        _release();
        _numRows = rows.length();
        _numCols = cols.length();
        _firstRow = rows.firstIndex();
        _firstCol = cols.firstIndex();
        _allocate(value);
        return true;
    }
    changeIndexBase(rows.firstIndex(), cols.firstIndex());
    return false;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
template <typename FS>
bool
FullStorage<T, Order, I, A>::resize(const FS &rhs, const ElementType &value)
{
    return resize(rhs.numRows(), rhs.numCols(),
                  rhs.firstRow(), rhs.firstCol(),
                  value);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
bool
FullStorage<T, Order, I, A>::fill(const ElementType &value)
{
    ASSERT(_data);
    std::fill_n(data(), numRows()*numCols(), value);
    return true;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
void
FullStorage<T, Order, I, A>::changeIndexBase(IndexType firstRow,
                                             IndexType firstCol)
{
    if (Order==cxxblas::RowMajor) {
        _data = data() - (firstRow*leadingDimension() + firstCol);
    }
    if (Order==cxxblas::ColMajor) {
        _data = data() - (firstCol*leadingDimension() + firstRow);
    }
    _firstRow = firstRow;
    _firstCol = firstCol;
}

// view of rectangular part
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ConstView
FullStorage<T, Order, I, A>::view(IndexType fromRow, IndexType fromCol,
                                  IndexType toRow, IndexType toCol,
                                  IndexType firstViewRow,
                                  IndexType firstViewCol) const
{
    ASSERT(fromRow>=firstRow());
    ASSERT(fromRow<=toRow);
    ASSERT(toRow<=lastRow());

    ASSERT(fromCol>=firstCol());
    ASSERT(fromCol<=toCol);
    ASSERT(toCol<=lastCol());

    return ConstView(&(this->operator()(fromRow, fromCol)),// data
                     allocator(),                          // allocator
                     toRow-fromRow+1,                      // # rows
                     toCol-fromCol+1,                      // # cols
                     leadingDimension(),                   // leading dimension
                     firstViewRow,                         // firstRow
                     firstViewCol);                        // firstCol
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::View
FullStorage<T, Order, I, A>::view(IndexType fromRow, IndexType fromCol,
                                  IndexType toRow, IndexType toCol,
                                  IndexType firstViewRow,
                                  IndexType firstViewCol)
{
    ASSERT(fromRow>=firstRow());
    ASSERT(fromRow<=toRow);
    ASSERT(toRow<=lastRow());

    ASSERT(fromCol>=firstCol());
    ASSERT(fromCol<=toCol);
    ASSERT(toCol<=lastCol());

    return View(&(this->operator()(fromRow, fromCol)),// data
                allocator(),                          // allocator
                toRow-fromRow+1,                      // # rows
                toCol-fromCol+1,                      // # cols
                leadingDimension(),                   // leading dimension
                firstViewRow,                         // firstRow
                firstViewCol);                        // firstCol
}

// view of single row
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ConstArrayView
FullStorage<T, Order, I, A>::viewRow(IndexType row,
                                     IndexType firstViewIndex) const
{
    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ConstArrayView(&(this->operator()(row, _firstCol)),
                          allocator(),
                          numCols(),
                          strideCol(),
                          firstViewIndex);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::ArrayView
FullStorage<T, Order, I, A>::viewRow(IndexType row,
                                     IndexType firstViewIndex)
{
    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ArrayView(&(this->operator()(row, _firstCol)),
                     allocator(),
                     numCols(),
                     strideCol(),
                     firstViewIndex);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ConstArrayView
FullStorage<T, Order, I, A>::viewRow(IndexType row,
                                     IndexType firstCol, IndexType lastCol,
                                     IndexType firstViewIndex) const
{
    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ConstArrayView(&(this->operator()(row, firstCol)),
                          allocator(),
                          lastCol-firstCol+1,
                          strideCol(),
                          firstViewIndex);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::ArrayView
FullStorage<T, Order, I, A>::viewRow(IndexType row,
                                     IndexType firstCol, IndexType lastCol,
                                     IndexType firstViewIndex)
{
    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ArrayView(&(this->operator()(row, firstCol)),
                     allocator(),
                     lastCol-firstCol+1,
                     strideCol(),
                     firstViewIndex);
}

// view of single column
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ConstArrayView
FullStorage<T, Order, I, A>::viewCol(IndexType col,
                                     IndexType firstViewIndex) const
{
    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ConstArrayView(&(this->operator()(_firstRow, col)),
                          allocator(),
                          numRows(),
                          strideRow(),
                          firstViewIndex);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::ArrayView
FullStorage<T, Order, I, A>::viewCol(IndexType col,
                                     IndexType firstViewIndex)
{
    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ArrayView(&(this->operator()(_firstRow, col)),
                     allocator(),
                     numRows(),
                     strideRow(),
                     firstViewIndex);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ConstArrayView
FullStorage<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
                                     IndexType col,
                                     IndexType firstViewIndex) const
{
    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ConstArrayView(&(this->operator()(firstRow, col)),
                          allocator(),
                          lastRow-firstRow+1,
                          strideRow(),
                          firstViewIndex);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::ArrayView
FullStorage<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
                                     IndexType col,
                                     IndexType firstViewIndex)
{
    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ArrayView(&(this->operator()(firstRow, col)),
                     allocator(),
                     lastRow-firstRow+1,
                     strideRow(),
                     firstViewIndex);
}

// view of d-th diagonal
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename FullStorage<T, Order, I, A>::ConstArrayView
FullStorage<T, Order, I, A>::viewDiag(IndexType d,
                                      IndexType firstViewIndex) const
{
    IndexType col = firstCol() + ( (d>0) ? d : 0 );
    IndexType row = firstRow() + ( (d>0) ? 0 : -d );
    return ConstArrayView(&(this->operator()(row,col)),
                          allocator(),
                          std::min(numRows(),numCols()) - std::abs(d),
                          leadingDimension()+1,
                          firstViewIndex);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename FullStorage<T, Order, I, A>::ArrayView
FullStorage<T, Order, I, A>::viewDiag(IndexType d,
                                      IndexType firstViewIndex)
{
    IndexType col = firstCol() + ( (d>0) ? d : 0 );
    IndexType row = firstRow() + ( (d>0) ? 0 : -d );
    return ArrayView(&(this->operator()(row,col)),
                     allocator(),
                     std::min(numRows(),numCols()) - std::abs(d),
                     leadingDimension()+1,
                     firstViewIndex);
}

//-- Private Methods -----------------------------------------------------------

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
void
FullStorage<T, Order, I, A>::_setIndexBase(IndexType firstRow,
                                           IndexType firstCol)
{
    // assume: _data points to allocated memory
    
    if (Order==cxxblas::RowMajor) {
        _data -= firstRow*_numCols + firstCol;
    }
    if (Order==cxxblas::ColMajor) {
        _data -= firstCol*_numRows + firstRow;
    }
    _firstRow = firstRow;
    _firstCol = firstCol;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
void
FullStorage<T, Order, I, A>::_raw_allocate()
{
    ASSERT(!_data);
    ASSERT(_numRows>0);
    ASSERT(_numCols>0);

    _data = _allocator.allocate(_numRows*_numCols);
#ifndef NDEBUG
    ElementType *p = _data;
#endif
    _setIndexBase(_firstRow, _firstCol);
    ASSERT(_data);
#ifndef NDEBUG
    ASSERT(p==data());
#endif
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
void
FullStorage<T, Order, I, A>::_allocate(const ElementType &value)
{
    _raw_allocate();
    T *p = data();
    for (IndexType i=0; i<numRows()*numCols(); ++i) {
        _allocator.construct(p++, value);
    }
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
void
FullStorage<T, Order, I, A>::_release()
{
    if (_data) {
        T *p = data();
        for (IndexType i=0; i<numRows()*numCols(); ++i) {
            _allocator.destroy(p++);
        }
        _allocator.deallocate(data(), numRows()*numCols());
        _data = 0;
    }
    ASSERT(_data==0);
}

} // namespace flens

#endif // FLENS_STORAGE_FULLSTORAGE_FULLSTORAGE_TCC
