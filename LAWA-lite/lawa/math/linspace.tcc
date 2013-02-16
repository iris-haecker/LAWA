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

namespace lawa {

template <typename T>
flens::DenseVector<flens::Array<T> >
linspace(T from, T to, int numTicks)
{
    T step = (to-from) / (numTicks-1);
    flens::DenseVector<flens::Array<T> > x(numTicks,0);
    x(0) = from;
    for (int i=1; i<numTicks-1; ++i) {
        x(i) = from + i*step;
    }
    x(x.lastIndex()) = to;
    return x;
}

} // namespace lawa

