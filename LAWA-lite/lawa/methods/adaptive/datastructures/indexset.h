/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEXSET_H_
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEXSET_H_ 1

#include <set>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/bspline.h>
#include <lawa/settings/enum.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename Index>
struct IndexSet : std::set<Index, lt<Lexicographical, Index > >
{
    using std::set<Index, lt<Lexicographical, Index > >::erase;
    using std::set<Index, lt<Lexicographical, Index > >::insert;

    IndexSet(void);

    IndexSet<Index>&
    operator= (const IndexSet<Index> &_set);

    IndexSet<Index>
    operator+ (const IndexSet<Index> &_set) const;
};

template <typename Index>
std::ostream& operator<< (std::ostream &s, const IndexSet<Index> &i);

void
getMinAndMaxLevel(const IndexSet<Index1D> &Lambda, int &jmin, int &jmax);

void
split(const IndexSet<Index2D> &Lambda, IndexSet<Index1D> &Lambda_x, IndexSet<Index1D> &Lambda_y);

IndexSet<Index1D>
extractSpaceIndices(const IndexSet<Index2D> &Lambda);

}   // namespace lawa


#include <lawa/methods/adaptive/datastructures/indexset.tcc>


#endif //  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEXSET_H

