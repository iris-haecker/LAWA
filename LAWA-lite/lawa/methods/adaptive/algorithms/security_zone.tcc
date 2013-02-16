#include <algorithm>

namespace lawa {

template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,Domain,Cons> &basis) {
    IndexSet<Index1D> ret, tmp;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        tmp = C((*lambda),c,basis);
        for (const_it mu=tmp.begin(); mu!=tmp.end(); ++mu) {
            if (Lambda.count(*mu) == 0) ret.insert(*mu);
        }
    }
    return ret;
}

template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const Index1D &lambda, T c, const Basis<T,Primal,Domain,Cons> &basis) {
    IndexSet<Index1D> ret;
    C(lambda,c,basis.mra,basis,ret);
    return ret;
}

template <typename T>
IndexSet<Index1D>
C_WO_XBSpline(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,R,CDF> &basis, bool only_pos)
{
    IndexSet<Index1D> ret, tmp;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        if (only_pos && (*lambda).j<=0) continue;
        tmp = C_WO_XBSpline((*lambda),c,basis);
        for (const_it mu=tmp.begin(); mu!=tmp.end(); ++mu) {
            if (Lambda.count(*mu) == 0) ret.insert(*mu);
        }
    }
    return ret;
}

template <typename T>
IndexSet<Index1D>
C_WO_XBSpline(const Index1D &lambda, T c, const Basis<T,Primal,R,CDF> &basis) {
    const Wavelet<T,Primal,R,CDF> psi(basis,0);
    int j = lambda.j, k = lambda.k;
    IndexSet<Index1D> ret;
    Support<T> contractedSupp, supp = basis.psi.support(j,k);
    T center = 0.5*(supp.l1 + supp.l2);
    contractedSupp.l1 = c*supp.l1 + (1-c)*center;
    contractedSupp.l2 = c*supp.l2 + (1-c)*center;
    long int kMin, kMax;

    kMin = floor( pow2i<T>(j-1)*contractedSupp.l1 - basis.psi.support(0,0).l2);
    kMax = ceil(pow2i<T>(j-1)*contractedSupp.l2 - basis.psi.support(0,0).l1);
    for (long int k1=kMin; k1<=kMax; ++k1) {
        if (overlap(contractedSupp, basis.psi.support(j-1,k1))>0) ret.insert(Index1D(j-1,k1,XWavelet));
    }

    kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.support(0,0).l2);
    kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - basis.psi.support(0,0).l1);
    for (long int k1=kMin; k1<=kMax; ++k1) {
        if (overlap(contractedSupp, basis.psi.support(j,k1))>0) ret.insert(Index1D(j,k1,XWavelet));
    }

    kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.support(0,0).l2);
    kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.support(0,0).l1);
    for (long int k1=kMin; k1<=kMax; ++k1) {
        if (overlap(contractedSupp, basis.psi.support(j+1,k1))>0) ret.insert(Index1D(j+1,k1,XWavelet));
    }
    return ret;
}

// Security zone interval
template <typename T, Construction Cons>
void
C(const Index1D &lambda, T c, const MRA<T,Primal,Interval,Cons> &/*mra*/,
  const Basis<T,Primal,Interval,Cons> &basis, IndexSet<Index1D> &ret)
{
    using std::min;
    using std::max;
    //ret.insert(Index1D(lambda.j,lambda.k,lambda.xtype));

    int j=lambda.j, jP1, k=lambda.k;
    XType xtype = lambda.xtype;
    int kMin_mra = basis.mra.rangeI(j).firstIndex(), kMax_mra = basis.mra.rangeI(j).lastIndex();
    if (lambda.xtype==XBSpline) {
        jP1=j;
        ret.insert(Index1D(j,std::max(k-2,kMin_mra),xtype));
        ret.insert(Index1D(j,std::max(k-1,kMin_mra),xtype));
        ret.insert(Index1D(j,std::min(k+1,kMax_mra),xtype));
        ret.insert(Index1D(j,std::min(k+2,kMax_mra),xtype));
    } else {
        jP1=j+1;
    }
    Support<T> supp = basis.generator(xtype).support(j,k);

    T zLambda=0.5*(supp.l2+supp.l1);
    Support<T> contractedSupp(c*supp.l1 + (1-c)*zLambda, c*supp.l2 + (1-c)*zLambda);

    int kMin = basis.rangeJ(jP1).firstIndex(), kMax = basis.rangeJ(jP1).lastIndex();
    int kStart = std::min(std::max(iceil(contractedSupp.l1 * pow2i<T>(jP1)), kMin), kMax);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kStart))>0));
    while ((kStart-1 >= kMin) && (overlap(contractedSupp, basis.psi.support(jP1,std::max(kStart-1, kMin)))>0)) {
        --kStart;
    }
    int kEnd = std::max(std::min(ifloor(contractedSupp.l2 * pow2i<T>(jP1)), kMax), kMin);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kEnd))>0));
    while ((kEnd+1 <= kMax) && (overlap(contractedSupp, basis.psi.support(jP1,std::min(kEnd+1, kMax)))>0)) {
        ++kEnd;
    }

    for (int k=kStart; k<=kEnd; ++k) {
        ret.insert(Index1D(jP1,k,XWavelet));
    }
}

// Security zone periodic
template <typename T>
void
C(const Index1D &lambda, T c, const MRA<T,Primal,Periodic,CDF> &mra,
  const Basis<T,Primal,Periodic,CDF> &basis, IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;
    //ret.insert(Index1D(j,k,xtype));
    if (xtype==XBSpline) {
        ret.insert(Index1D(j,(k-1 >= mra.rangeI(j).firstIndex()) ? k-1 : mra.rangeI(j).lastIndex()
                                    + ((1 - (mra.rangeI(j).firstIndex() - k+1))%mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+1 <= mra.rangeI(j).lastIndex() ? k+1 : mra.rangeI(j).firstIndex()
                                    - ((1 - (k+1 - mra.rangeI(j).lastIndex()))%mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k-2 >= mra.rangeI(j).firstIndex() ? k-2 : mra.rangeI(j).lastIndex()
                                    + ((1 - (mra.rangeI(j).firstIndex() - k+2))%mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+2 <= mra.rangeI(j).lastIndex() ? k+2 : mra.rangeI(j).firstIndex()
                                    - ((1 - (k+2 - mra.rangeI(j).lastIndex()))%mra.cardI(j)),xtype));
        Support<T> contractedSupp, supp = mra.phi.phiR.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp = c*supp + (1-c)*center;

        int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

        for (int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
                int k = k1;
                if(k < basis.rangeJ(j).firstIndex()){
                    k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
                }
                if(k > basis.rangeJ(j).lastIndex()){
                    k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
                }
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    else {
        Support<T> contractedSupp, supp = basis.psi.psiR.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp = c*supp + (1-c)*center;
        /*    no wavelet indices on the same level?!
        long int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
                int k = k1;
                if(k < basis.rangeJ(j).firstIndex()){
                    k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
                }
                if(k > basis.rangeJ(j).lastIndex()){
                    k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
                }
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
        */

        long int kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.psiR.support(j+1,k1))>0)
            {
                int k = k1;
                if(k < basis.rangeJ(j+1).firstIndex()){
                    k = basis.rangeJ(j+1).lastIndex() + ((1 - (basis.rangeJ(j+1).firstIndex() - k))%basis.cardJ(j+1));
                }
                if(k > basis.rangeJ(j+1).lastIndex()){
                    k = basis.rangeJ(j+1).firstIndex() - ((1 - (k - basis.rangeJ(j+1).lastIndex()))%basis.cardJ(j+1));
                }
                ret.insert(Index1D(j+1,k,XWavelet));
            }
        }

    }
}

// Security zone realline
template <typename T>
void
C(const Index1D &lambda, T c, const MRA<T,Primal,R,CDF> &mra,
  const Basis<T,Primal,R,CDF> &basis, IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;
    if (xtype==XBSpline) {
        ret.insert(Index1D(j,k-1,xtype));
        ret.insert(Index1D(j,k+1,xtype));
        ret.insert(Index1D(j,k-2,xtype));
        ret.insert(Index1D(j,k+2,xtype));

        Support<T> contractedSupp, supp = mra.phi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;

        int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.support(0,0).l2);
        int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - basis.psi.support(0,0).l1);
        for (int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.support(j,k1))>0) ret.insert(Index1D(j,k1,XWavelet));
        }
    }
    else {
        Support<T> contractedSupp, supp = basis.psi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;
        long int kMin, kMax;

        kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.support(0,0).l2);
        kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.support(0,0).l1);
        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.support(j+1,k1))>0) ret.insert(Index1D(j+1,k1,XWavelet));
        }
    }
}


//Security zone 2D
template <typename T, typename Basis2D>
IndexSet<Index2D>
C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis)
{
    typedef typename IndexSet<Index2D>::const_iterator const_it_2d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    IndexSet<Index2D>  ret;

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_2d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1, C_index2;
        C_index1 = C((*lambda).index1, c, basis.first);
        C_index2 = C((*lambda).index2, c, basis.second);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index2D((*it_C_index1), (*lambda).index2))==0) {
                ret.insert(Index2D((*it_C_index1), (*lambda).index2));
            }
        }
        for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
            if (Lambda.count(Index2D((*lambda).index1, (*it_C_index2)))==0) {
                ret.insert(Index2D((*lambda).index1, (*it_C_index2)));
            }
        }
/*
        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                if (Lambda.count(Index2D((*it_C_index1), (*it_C_index2)))==0) {
                    ret.insert(Index2D((*it_C_index1), (*it_C_index2)));
                }
            }
        }
*/
    }
    return ret;
}

template <typename T, typename Basis2D>
IndexSet<Index2D>
C_t(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis)
{
    typedef typename IndexSet<Index2D>::const_iterator const_it_2d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    IndexSet<Index2D>  ret;

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_2d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1, C_index2;
        C_index1 = C((*lambda).index1, c, basis.first);
        //C_index2 = C((*lambda).index2, c, basis.second);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index2D((*it_C_index1), (*lambda).index2))==0) {
                ret.insert(Index2D((*it_C_index1), (*lambda).index2));
            }
        }
    }
    return ret;
}


//Security zone 3D
template <typename T, typename Basis3D>
IndexSet<Index3D>
C(const IndexSet<Index3D> &Lambda, T c, const Basis3D &basis)
{
    typedef typename IndexSet<Index3D>::const_iterator const_it_3d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    IndexSet<Index3D>  ret;

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_3d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1, C_index2, C_index3;
        C_index1 = C((*lambda).index1, c, basis.first);
        C_index2 = C((*lambda).index2, c, basis.second);
        C_index3 = C((*lambda).index3, c, basis.third);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index3D((*it_C_index1), (*lambda).index2, (*lambda).index3) )==0) {
                ret.insert(Index3D((*it_C_index1),   (*lambda).index2,  (*lambda).index3) );
            }
        }
        for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
            if (Lambda.count(Index3D((*lambda).index1, (*it_C_index2), (*lambda).index3) )==0) {
                ret.insert(Index3D((*lambda).index1,   (*it_C_index2), (*lambda).index3) );
            }
        }
        for (const_it it_C_index3=C_index3.begin(); it_C_index3!=C_index3.end(); ++it_C_index3) {
            if (Lambda.count(Index3D((*lambda).index1, (*lambda).index2, (*it_C_index3) ) )==0) {
                ret.insert(Index3D((*lambda).index1,   (*lambda).index2, (*it_C_index3) ) );
            }
        }
    }
    return ret;
}

} // namespace lawa

