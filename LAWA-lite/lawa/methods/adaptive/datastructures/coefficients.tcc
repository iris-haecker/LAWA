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

namespace lawa {

/*
 * Coefficient vector ordered by index
 */

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>::Coefficients()
{

}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>&
Coefficients<Lexicographical,T,Index>::operator=(const Coefficients<Lexicographical,T,Index> &_coeff)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::value_type val_type;
    erase(Coefficients<Lexicographical,T,Index>::begin(), Coefficients<Lexicographical,T,Index>::end());

    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).first, (*lambda).second));
        }
    }

    return *this;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>&
Coefficients<Lexicographical,T,Index>::operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::value_type val_type;
    erase(Coefficients<Lexicographical,T,Index>::begin(), Coefficients<Lexicographical,T,Index>::end());

    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).second, (*lambda).first));
        }
    }

    return *this;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
Coefficients<Lexicographical,T,Index>::operator-(const Coefficients<Lexicographical,T,Index> &_coeff) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret = *this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            ret.operator[]((*lambda).first) -= (*lambda).second;
        }
    }
    return ret;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
Coefficients<Lexicographical,T,Index>::operator+(const Coefficients<Lexicographical,T,Index> &_coeff) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret = *this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            ret.operator[]((*lambda).first) += (*lambda).second;
        }
    }
    return ret;
}

template <typename T, typename Index>
T
Coefficients<Lexicographical,T,Index>::operator*(const Coefficients<Lexicographical,T,Index> &_coeff2) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> _coeff1 = *this;
    T ret = 0.0;
    if (_coeff2.size() > 0) {
        for (const_it lambda = _coeff2.begin(); lambda != _coeff2.end(); ++lambda) {
            ret += _coeff1.operator[]((*lambda).first) * (*lambda).second;
        }
    }
    return ret;
}

template <typename T, typename Index>
T
Coefficients<Lexicographical,T,Index>::norm(T tau) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    T result=0.0;
    if (Coefficients<Lexicographical,T,Index>::size() > 0) {
        for (const_it mu=Coefficients<Lexicographical,T,Index>::begin();
             mu!=Coefficients<Lexicographical,T,Index>::end(); ++mu) {
            result+=std::pow((T)fabs((*mu).second), (T)tau);
        }
    }

    return std::pow(result, 1.0/tau);
}

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<Lexicographical,T,Index> &c)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    s << std::endl << "Coefficient_Rn<Lexicographical,T,Index>:" << std::endl;
    if (c.size() > 0) {
        for (const_it lambda = c.begin(); lambda != c.end(); ++lambda) {
            s << "  [" << (*lambda).first << "\t | " << (*lambda).second
                << "]" << std::endl;
        }
    }
    return s << std::endl;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
operator*(T alpha, const Coefficients<Lexicographical,T,Index> &_coeff) {
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret(_coeff.d,_coeff.d_);
    for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
        ret.operator[]((*lambda).first) = alpha*(*lambda).second;
    }
    return ret;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
P(const Coefficients<Lexicographical,T,Index> &v, const IndexSet<Index> &Lambda)
{
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret;
    const_it v_end = v.end();
    for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        const_it it = v.find(*lambda);
        if (it != v_end) {
            ret[*lambda] = (*it).second;
        }
        else {
            ret[*lambda] = 0.;
        }
    }
    return ret;
}

template <typename T, typename Index>
IndexSet<Index>
supp(const Coefficients<Lexicographical,T,Index> &v)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    IndexSet<Index> ret;
    if (v.size() > 0) {
        for (const_it lambda = v.begin(); lambda != v.end(); ++lambda) {
            ret.insert((*lambda).first);
        }
    }

    return ret;
}

template <typename T, typename Index>
void
FillWithZeros(const IndexSet<Index> &Lambda, Coefficients<Lexicographical,T,Index> &u)
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    for (const_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        if (u.count((*lambda)) == 0) u[*lambda] = 0.0;
    }
}


/*
 * Coefficient vector ordered by absolute value
 */
template <typename T, typename Index>
Coefficients<AbsoluteValue,T,Index>::Coefficients()
{

}

template <typename T, typename Index>
Coefficients<AbsoluteValue,T,Index>&
Coefficients<AbsoluteValue,T,Index>::operator= (const Coefficients<Lexicographical,T,Index> &_coeff)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;

    erase(Coefficients<AbsoluteValue,T,Index>::begin(), Coefficients<AbsoluteValue,T,Index>::end());
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).second, (*lambda).first));
        }
    }
    return *this;
}

template <typename T, typename Index>
Coefficients<AbsoluteValue,T,Index>&
Coefficients<AbsoluteValue,T,Index>::operator= (const Coefficients<AbsoluteValue,T,Index> &_coeff)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;

    erase(Coefficients<AbsoluteValue,T,Index>::begin(), Coefficients<AbsoluteValue,T,Index>::end());
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).first, (*lambda).second));
        }
    }
    return *this;
}

template <typename T, typename Index>
T
Coefficients<AbsoluteValue,T,Index>::norm(T tau) const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    T result=0.0;
    if (Coefficients<AbsoluteValue,T,Index>::size() > 0) {
        for (const_it mu=Coefficients<AbsoluteValue,T,Index>::begin();
             mu!=Coefficients<AbsoluteValue,T,Index>::end(); ++mu) {
            result+=std::pow((T)fabs((*mu).first),(T) tau);
        }
    }

    return std::pow(result, 1.0/tau);
}

template <typename T, typename Index>
T
Coefficients<AbsoluteValue,T,Index>::l2bestnterm(int n) const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_reverse_iterator const_rev_it;
    T result=0.0;
    int count=(*this).size();
    for (const_rev_it lambda=Coefficients<AbsoluteValue,T,Index>::rbegin(); lambda!=Coefficients<AbsoluteValue,T,Index>::rend(); ++lambda) {
        if (count>=n) {
            result += fabs((*lambda).first) * fabs((*lambda).first);
        }
        --count;
    }
    return sqrt(result);
}

template <typename T, typename Index>
T
Coefficients<AbsoluteValue,T,Index>::wtauNorm(T tau) const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    T result=0.0, temp=0.0;
    int n=1;
    for (const_it lambda=Coefficients<AbsoluteValue,T,Index>::begin(); lambda!=Coefficients<AbsoluteValue,T,Index>::end(); ++lambda) {
        temp=std::pow(T(n),1.0/tau) * fabs((*lambda).first);
        if (temp>result) {
            result=temp;
        }
        ++n;
    }

    return result + norm();
}

template <typename T, typename Index>
DenseVector<Array<T> >
Coefficients<AbsoluteValue,T,Index>::norm_sections() const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;

    DenseVector<Array<T> > result;

    if (Coefficients<AbsoluteValue,T,Index>::size() > 0) {
        //result.engine().resizeOrClear(_(0, int(log(T(Coefficients<AbsoluteValue,T,Index>::size()))/log(T(2))+1)));
        result.engine().resize(int(log(T(Coefficients<AbsoluteValue,T,Index>::size()))/log(T(2))+1));
        T help=0.0;
        int count=0, s=0;

        for (const_it lambda = Coefficients<AbsoluteValue,T,Index>::begin(); lambda != Coefficients<AbsoluteValue,T,Index>::end(); ++lambda) {
            help += ((*lambda).first)*((*lambda).first);
            count++;

            if (count==(1<<s)) {
                result(s+1) = sqrt(help);
                s++;
                help=0.0;
            }
        }
        if (help != 0.0) {
            result(s)=sqrt(help);
        }
    } else {
        result.engine().resize(0);
    }
    return result;
}

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<AbsoluteValue,T,Index> &c)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    s << std::endl << "Coefficients<AbsoluteValue,T,Index>:" << std::endl;
    if (c.size() > 0) {
        for (const_it it=c.begin(); it!=c.end(); ++it) {
            s << "  [" << (*it).first << "\t | " << (*it).second << "]" << std::endl;
        }
    }
    return s << std::endl;
}


}   //namespace lawa

