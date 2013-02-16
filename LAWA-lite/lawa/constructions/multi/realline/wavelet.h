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

#ifndef LAWA_CONSTRUCTIONS_MULTI_REALLINE_WAVELET_H
#define LAWA_CONSTRUCTIONS_MULTI_REALLINE_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Wavelet<_T,Orthogonal,R,Multi>
    : public BasisFunction<_T,Orthogonal,R,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = Multi;
    
        Wavelet(int _d);
        
//TODO        Wavelet(const Basis<T,Orthogonal,R,Multi> &basis);
        
        virtual
        ~Wavelet();
        
        T
        operator()(T x, int j, long k, unsigned short deriv) const;
        
        Support<T>
        support(int j, int k) const;
        
        DenseVector<Array<T> >
        singularSupport(int j, int k) const;
        
//TODO        int polynomialOrder;
        const int d;
        const int vanishingMoments;

    private:
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        long
        _shift(long k) const;
        
        int
        _type(long k) const;
    
        unsigned int _numSplines;
        Evaluator *_evaluator;
        Support<T> *_support;
        DenseVector<Array<T> > *_singularSupport;

    
};
    
} // namespace lawa

#include <lawa/constructions/multi/realline/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_MULTI_REALLINE_WAVELET_H
