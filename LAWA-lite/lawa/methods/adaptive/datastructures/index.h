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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEX_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEX_H 1

#include <lawa/settings/enum.h>
#include <iostream>

namespace lawa {

struct Index1D
{
    short j;
    int k;
    XType xtype;
    long val;
    mutable unsigned int linearindex;

    Index1D(void);
    Index1D(int j, int k, XType _xtype);
    Index1D(const Index1D &index);
    ~Index1D();
};

std::ostream& operator<<(std::ostream &s, const Index1D &_Index);

struct Index2D
{
    mutable unsigned int linearindex;
    Index2D(const Index1D &index1, const Index1D &index2);
    ~Index2D();
    Index1D index1, index2;

};

std::ostream& operator<<(std::ostream &s, const Index2D &_Index);

struct Index3D
{
    mutable unsigned int linearindex;
    Index3D(const Index1D &index1, const Index1D &index2, const Index1D &index3);
    ~Index3D();
    Index1D index1, index2, index3;

};

std::ostream& operator<<(std::ostream &s, const Index3D &_Index);

template <typename Index>
struct Entry
{
    Entry(const Index &row_index, const Index &col_index);
    const Index row_index, col_index;    //todo: no copy, but only a reference possible ?!
};

template <typename Index>
std::ostream& operator<<(std::ostream &s, const Entry<Index> &entry);

template <SortingCriterion S, typename SortingType>
struct lt
{
};

//Bitmask implementation
template<> 
struct lt<Lexicographical, Index1D>
{
    bool operator()(const Index1D &left, const Index1D &right) const;

    bool operator()(const Entry<Index1D> &left, const Entry<Index1D> &right) const;
};

//-- Index2D: Lexicographical --------------------------------------------------

template <>
struct lt<Lexicographical, Index2D>
{
    bool operator()(const Index2D &left, const Index2D &right) const;

    bool operator()(const Entry<Index2D> &left, const Entry<Index2D> &right) const;
};

/*
//-- Index2D: LexicographicalByRows --------------------------------------------

template <>
struct lt<LexicographicalByRows, Index2D>
{
    bool operator()(const Index2D &left, const Index2D &right) const;

    bool operator()(const Entry<Index2D> &left, const Entry<Index2D> &right) const;
};

//-- Index2D: LexicographicalByCols --------------------------------------------

template <>
struct lt<LexicographicalByCols, Index2D>
{
    bool operator()(const Index2D &left, const Index2D &right) const;

    bool operator()(const Entry<Index2D> &left, const Entry<Index2D> &right) const;
};
*/

//-- Index3D: Lexicographical --------------------------------------------------

template <>
struct lt<Lexicographical, Index3D>
{
    bool operator()(const Index3D &left, const Index3D &right) const;

    bool operator()(const Entry<Index3D> &left, const Entry<Index3D> &right) const;
};

} //namespace lawa

#include <lawa/methods/adaptive/datastructures/index.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEX_H
