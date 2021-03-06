/*
 *   Copyright (c) 2004, Alexander Stippler
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

#ifndef FLENS_AUX_RANGE_TCC
#define FLENS_AUX_RANGE_TCC 1

namespace flens {

template <typename IndexType>
Range<IndexType>::Range(IndexType firstIndex, IndexType lastIndex)
    : _firstIndex(firstIndex), _stride(1), _lastIndex(lastIndex)
{
}

template <typename IndexType>
Range<IndexType>::Range(IndexType firstIndex,
                        IndexType stride,
                        IndexType lastIndex)
    : _firstIndex(firstIndex), _stride(stride), _lastIndex(lastIndex)
{
}

template <typename IndexType>
IndexType
Range<IndexType>::firstIndex() const
{
    return _firstIndex;
}

template <typename IndexType>
IndexType
Range<IndexType>::stride() const
{
    return _stride;
}

template <typename IndexType>
IndexType
Range<IndexType>::lastIndex() const
{
    return _lastIndex;
}

template <typename IndexType>
IndexType
Range<IndexType>::numTicks() const
{
    return (_lastIndex-_firstIndex)/_stride + IndexType(1);
}

template <typename IndexType>
IndexType
Range<IndexType>::length() const
{
    return (_lastIndex-_firstIndex)/_stride + IndexType(1);
}

//------------------------------------------------------------------------------

template <typename IndexType>
std::ostream &
operator<<(std::ostream &out, const Range<IndexType> &range)
{
    out << "[" << range.firstIndex() << ","; 
    if (range.stride()!=1) {
        out << range.stride();
    }
    out << range.lastIndex() << "]";
    return out;
}

template <typename IndexType>
bool
operator==(const Range<IndexType> &lhs, const Range<IndexType> &rhs)
{
    return (lhs.firstIndex()==rhs.firstIndex())
        && (lhs.lastIndex()==rhs.lastIndex())
        && (lhs.stride()==rhs.stride());
}


} // namespace flens

#endif // FLENS_AUX_RANGE_TCC
