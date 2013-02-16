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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFICIENTS_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFICIENTS_H 1

#include <map>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa {

template <SortingCriterion S, typename T, typename Index>
struct Coefficients
{
};

template <typename T, typename Index>
struct Coefficients<Lexicographical,T,Index> : std::map<Index,T,lt<Lexicographical,Index> >
{
    using std::map<Index,T,lt<Lexicographical,Index> >::insert;
    using std::map<Index,T,lt<Lexicographical,Index> >::erase;
    
    Coefficients();        //required in rhs.h

    Coefficients<Lexicographical,T,Index>&
    operator=(const Coefficients<Lexicographical,T,Index> &_coeff);

    Coefficients<Lexicographical,T,Index>&
    operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff);

    Coefficients<Lexicographical,T,Index>
    operator-(const Coefficients<Lexicographical,T,Index> &_coeff) const;

    Coefficients<Lexicographical,T,Index>
    operator+(const Coefficients<Lexicographical,T,Index> &_coeff) const;

    T
    operator*(const Coefficients<Lexicographical,T,Index> &_coeff) const;

    T
    norm(T tau=2.0) const;
};

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<Lexicographical,T,Index> &c);

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
operator*(T alpha, const Coefficients<Lexicographical,T,Index> &_coeff);

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
P(const Coefficients<Lexicographical,T,Index> &v, const IndexSet<Index> &Lambda);

template <typename T, typename Index>
IndexSet<Index>
supp(const Coefficients<Lexicographical,T,Index> &v);

template <typename T, typename Index>
void
FillWithZeros(const IndexSet<Index> &Lambda, Coefficients<Lexicographical,T,Index> &_coeff);


template <typename T, typename Index>
struct Coefficients<AbsoluteValue,T,Index> : std::multimap<T,Index,lt<AbsoluteValue,T> >
{
    using std::multimap<T,Index,lt<AbsoluteValue,T> >::insert;
    using std::multimap<T,Index,lt<AbsoluteValue,T> >::erase;
    
    Coefficients();

    Coefficients<AbsoluteValue,T,Index>&
    operator=(const Coefficients<Lexicographical,T,Index> &_coeff);

    Coefficients<AbsoluteValue,T,Index>&
    operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff);

    T
    norm(T tau=2.0) const;

    T
    l2bestnterm(int n) const;

    T
    wtauNorm(T tau) const;

    DenseVector<Array<T> >
    norm_sections() const;
};

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<AbsoluteValue,T,Index> &c);

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/coefficients.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFICIENTS_H

