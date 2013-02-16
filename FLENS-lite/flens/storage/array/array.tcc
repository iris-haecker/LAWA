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

#ifndef FLENS_STORAGE_ARRAY_ARRAY_TCC
#define FLENS_STORAGE_ARRAY_ARRAY_TCC 1

#include <cassert>
#include <cxxblas/level1/copy.h>
#include <flens/storage/array/arrayview.h>
#include <flens/storage/array/constarrayview.h>

namespace flens {

template <typename T, typename I, typename A>
Array<T, I, A>::Array()
    : _data(0), _length(0), _firstIndex(1)
{
}

template <typename T, typename I, typename A>
Array<T, I, A>::Array(IndexType length, IndexType firstIndex,
                      const ElementType &value, const Allocator &allocator)
    : _data(0), _allocator(allocator), _length(length), _firstIndex(firstIndex)
{
    _allocate(value);
}

template <typename T, typename I, typename A>
Array<T, I, A>::Array(const Array &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _length(rhs.length()), _firstIndex(rhs.firstIndex())
{
    if (length()>0) {
        _allocate();
        cxxblas::copy(length(), rhs.data(), rhs.stride(), data(), stride());
    }
}

template <typename T, typename I, typename A>
template <typename RHS>
Array<T, I, A>::Array(const RHS &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _length(rhs.length()), _firstIndex(rhs.firstIndex())
{
    ASSERT(length()>0);
    _allocate();
    cxxblas::copy(length(), rhs.data(), rhs.stride(), data(), stride());
}

template <typename T, typename I, typename A>
Array<T, I, A>::~Array()
{
    _release();
}

//-- operators -----------------------------------------------------------------

template <typename T, typename I, typename A>
const typename Array<T, I, A>::ElementType &
Array<T, I, A>::operator()(IndexType index) const
{
    ASSERT(index>=firstIndex());
    ASSERT(index<=lastIndex());
    return _data[index];
}

template <typename T, typename I, typename A>
typename Array<T, I, A>::ElementType &
Array<T, I, A>::operator()(IndexType index)
{
    ASSERT(index>=firstIndex());
    ASSERT(index<=lastIndex());
    return _data[index];
}

template <typename T, typename I, typename A>
typename Array<T, I, A>::IndexType
Array<T, I, A>::firstIndex() const
{
    return _firstIndex;
}

template <typename T, typename I, typename A>
typename Array<T, I, A>::IndexType
Array<T, I, A>::lastIndex() const
{
    return _firstIndex+_length-IndexType(1);
}

template <typename T, typename I, typename A>
typename Array<T, I, A>::IndexType
Array<T, I, A>::length() const
{
    return _length;
}

template <typename T, typename I, typename A>
typename Array<T, I, A>::IndexType
Array<T, I, A>::stride() const
{
    return IndexType(1);
}

template <typename T, typename I, typename A>
const typename Array<T, I, A>::ElementType *
Array<T, I, A>::data() const
{
    return &_data[_firstIndex];
}

template <typename T, typename I, typename A>
typename Array<T, I, A>::ElementType *
Array<T, I, A>::data()
{
    return &_data[_firstIndex];
}

template <typename T, typename I, typename A>
const typename Array<T, I, A>::Allocator &
Array<T, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, typename I, typename A>
bool
Array<T, I, A>::resize(IndexType length, IndexType firstIndex,
                       const ElementType &value)
{
    if (length!=_length) {
        _release();
        _length = length;
        _firstIndex = firstIndex;
        _allocate(value);
        return true;
    }
    changeIndexBase(firstIndex);
    return false;
}

template <typename T, typename I, typename A>
template <typename ARRAY>
bool
Array<T, I, A>::resize(const ARRAY &rhs, const ElementType &value)
{
    return resize(rhs.length(), rhs.firstIndex(), value);
}

template <typename T, typename I, typename A>
bool
Array<T, I, A>::fill(const ElementType &value)
{
    std::fill_n(data(), length(), value);
    return true;
}

template <typename T, typename I, typename A>
void
Array<T, I, A>::changeIndexBase(IndexType firstIndex)
{
    _data += _firstIndex - firstIndex;
    _firstIndex = firstIndex;
}

template <typename T, typename I, typename A>
const typename Array<T, I, A>::ConstView
Array<T, I, A>::view(IndexType from, IndexType to,
                     IndexType stride, IndexType firstViewIndex) const
{
    ASSERT(firstIndex()<=from);
    ASSERT(lastIndex()>=to);
    ASSERT(from<=to);
    ASSERT(stride>=1);
    return ConstView(&operator()(from),     // data
                     allocator(),           // allocator
                     (to-from)/stride+1,    // length
                     stride,                // stride
                     firstViewIndex);       // firstIndex
}

template <typename T, typename I, typename A>
typename Array<T, I, A>::View
Array<T, I, A>::view(IndexType from, IndexType to,
                     IndexType stride, IndexType firstViewIndex)
{
    ASSERT(firstIndex()<=from);
    ASSERT(lastIndex()>=to);
    ASSERT(from<=to);
    ASSERT(stride>=1);
    return View(&operator()(from),          // data
                allocator(),                // allocator
                (to-from)/stride+1,         // length
                stride,                     // stride
                firstViewIndex);            // firstIndex in view
}

//-- private methods -----------------------------------------------------------

template <typename T, typename I, typename A>
void
Array<T, I, A>::_raw_allocate()
{
    ASSERT(!_data);
    ASSERT(length()>=0);

    if (length()>0) {
        _data = _allocator.allocate(_length) - _firstIndex;
        ASSERT(_data+_firstIndex);
        ASSERT(_data);
    }
}

template <typename T, typename I, typename A>
void
Array<T, I, A>::_allocate(const ElementType &value)
{
    _raw_allocate();
    for (IndexType i=firstIndex(); i<=lastIndex(); ++i) {
        _allocator.construct(_data+i, value);
    }
}

template <typename T, typename I, typename A>
void
Array<T, I, A>::_release()
{
    if (_data) {
        for (IndexType i=firstIndex(); i<=lastIndex(); ++i) {
            _allocator.destroy(_data+i);
        }
        _allocator.deallocate(_data + _firstIndex, _length);   
        _data = 0;     
    }
    ASSERT(_data==0);
}

} // namespace flens

#endif // FLENS_STORAGE_ARRAY_ARRAY_TCC
