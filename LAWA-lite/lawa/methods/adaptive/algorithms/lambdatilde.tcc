namespace lawa {

/*
 * Realizations of lambdaTilde1d
 */

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis, 
                  int s_tilde, int jmin, int jmax, bool update)
{
    const BSpline<T,Primal,R,CDF> phi = basis.mra.phi;
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;

    int j = lambda.j, k = lambda.k;
    int d = psi.d;
    IndexSet<Index1D> ret;
    Support<T> support_refbspline = phi.support(0,0);
    Support<T> support_refwavelet = psi.support(0,0);

    if (!update) {

        if (lambda.xtype == XBSpline) {
            Support<T> supp = phi.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << phi_col.singularSupport(j,k) << endl;

            // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
            BSpline<T,Primal,R,CDF> phi_row(d);
            int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2)-1;
            int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1)+1;
            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                if (overlap(supp, phi.support(j,k_row)) > 0) {
                    //std::cout << "lambdaTilde: BSpline (" << j << ", " << k_row << "): " << phi_row.support(j,k_row) << " " << supp  << std::endl;
                    ret.insert(Index1D(j,k_row,XBSpline));
                }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {        // realization of matrix compression via level threshold
                T Pow2i_Mjrow = pow2i<T>(-j_row);
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = phi.singularSupport(j,k);
                    //cout << "LambdaTilde: Singular support phi_col = " << singpts;
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                            if (((overlap(supp_row, supp) > 0)) && (!(distance(singsupp,supp_row) > 0 ))) {
                                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp  << std::endl;
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if (overlap(supp, supp_row) > 0)  {
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
        else {
            Support<T> supp = psi.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

            // Inserting all indices corresponding to Bsplines with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            if (fabs(j - jmin) <= s_tilde) {
                int kMin = floor( pow2i<T>(jmin)*supp.l1 - phi.support(0,0).l2)-1;
                int kMax =  ceil( pow2i<T>(jmin)*supp.l2 - phi.support(0,0).l1)+1;
                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    if (overlap(supp, phi.support(jmin,k_row)) > 0) {
                        //std::cout << "lambdaTilde: BSpline (" << jmin << ", " << k_row << "): " << phi.support(jmin,k_row) << " " << supp  << std::endl;
                        ret.insert(Index1D(jmin,k_row,XBSpline));
                    }
                }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
                T Pow2i_Mjrow = pow2i<T>(-j_row);
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = psi.optim_singularSupport(j,k);
                    //cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                            if ((overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) > 0 ))){
                                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if ((overlap(supp, supp_row) > 0) && (!(distance(psi.optim_singularSupport(j_row,k_row),supp) > 0 ))) {
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }
                }
            }
        }
    }
/*
    // We assume that indices corresponding to non-zero values are already calculated up to level jmax-1 > jmin.
    // Therefore, we only have to add wavelets on level jmax.
    else {
        int j = lambda.j, k = lambda.k;
        int d = lambda.d, d_= lambda.d_;

        if (lambda.xtype == XBSpline) {
            assert(j == jmin);
            BSpline<T,Primal,R> phi_col(d);
            Support<T> supp = phi_col.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " <<  phi_col.support(j,k) << " " << phi_col.singularSupport(j,k) << endl;

            if (jmax <= j+s_tilde) {
                // Inserting all indices corresponding to Wavelets with intersecting support using
                // a) local compactness  b) matrix compression  c) vanishing moments
                Wavelet<T,Primal,R> psi_row(d,d_);
                DenseVector<Array<T> > singpts = phi_col.singularSupport(j,k);
                for (int i=singpts.firstIndex(); i<=singpts.lastIndex(); ++i) {
                    int kMin =  ceil(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l2);
                    int kMax = floor(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l1);
                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        if ((overlap(psi_row.support(jmax,k_row), supp) == true)) {
                            //cout << "LambdaTilde: kMin = " << kMin << ", kMax = " << kMax << ", k_row = " << k_row << ": " << psi_row.support(j_row,k_row) << endl;
                            ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
                        }
                    }
                }
            }
        }

        else {
            assert(j >= jmin);
            Wavelet<T,Primal,R> psi_col(d,d_);
            Support<T> supp = psi_col.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            Wavelet<T,Primal,R> psi_row(d,d_);
            if (jmax <= j+s_tilde) {
                int level_diff = 3;
                if (d == 2) level_diff = 2;
                if (d == 3) level_diff = 4;
                if (jmax >= j+level_diff) {                                                // level difference has to be large enough for vanishing entries due to vanishing moments
                    DenseVector<Array<T> > singpts = psi_col.singularSupport(j,k);
                    //cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
                    for (int i=singpts.firstIndex(); i<=singpts.lastIndex(); ++i) {
                        int kMin =  ceil(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l2);
                        int kMax = floor(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l1);
                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            if ((overlap(psi_row.support(jmax,k_row), supp) == true)) {
                                //cout << "LambdaTilde update: jmax = " << jmax << ", kMin = " << kMin << ", kMax = " << kMax << ", k_row = " << k_row << ": " << psi_row.support(jmax,k_row) << endl;
                                ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = ceil( pow2i(jmax)*supp.l1 - psi_row.support(0,0).l2);
                    int kMax = floor(pow2i(jmax)*supp.l2 - psi_row.support(0,0).l1);
                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        if ((overlap(supp, psi_row.support(jmax,k_row)) == true) && !((distSingularSupport(psi_row,psi_col,jmax,k_row,j,k) > 0.0 ))){
                            //cout << "lambdaTilde_R: Wavelet (" << j_row << ", " << k_row << "): " << psi_row.support(j_row,k_row) << endl;
                            ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
                        }
                    }
                }
            }
        }


    }
    */
    return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Periodic,CDF> &basis, 
                  int s_tilde, int jmin, int jmax, bool update)
{
    BSpline<T,Primal,Periodic,CDF> phi_col(basis.mra), phi_row(basis.mra);
    Wavelet<T,Primal,Periodic,CDF> psi_col(basis), psi_row(basis);
    int j = lambda.j, k = lambda.k;

    IndexSet<Index1D> ret;
    Support<T> support_refbspline = phi_col.phiR.support(0,0);
    Support<T> support_refwavelet = psi_col.psiR.support(0,0);

    if (!update) {

        if (lambda.xtype == XBSpline) {
            Support<T> supp = phi_col.phiR.support(j,k);

            // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
            int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2)-1;
            int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1)+1;
            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                int k_row_per = k_row;
                if(k_row_per < basis.mra.rangeI(j).firstIndex()){
                    k_row_per = basis.mra.rangeI(j).lastIndex() + ((1 - (basis.mra.rangeI(j).firstIndex() - k_row_per))%basis.mra.cardI(j));
                }
                if(k_row_per > basis.mra.rangeI(j).lastIndex()){
                    k_row_per = basis.mra.rangeI(j).firstIndex() - ((1 - (k_row_per - basis.mra.rangeI(j).lastIndex()))%basis.mra.cardI(j));
                }
                if (overlap(phi_row.support(j,k_row_per), phi_col.support(j,k)) > 0) {
                    ret.insert(Index1D(jmin,k_row_per,XBSpline));
               }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {        // realization of matrix compression via level threshold
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = phi_col.phiR.singularSupport(j,k);
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            int k_row_per = k_row;
                            if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                                k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                               }
                               if(k_row_per> basis.rangeJ(j_row).lastIndex()){
                                   k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                               }
                               Support<T> supp_row = psi_row.support(j_row,k_row_per);
                            if (((overlap(supp_row, phi_col.support(j,k)) > 0)) &&
                                (!(distance(phi_col.singularSupport(j,k),supp_row) >= 0 ))) {
                                ret.insert(Index1D(j_row,k_row_per,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        int k_row_per = k_row;
                        if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                            k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                        }
                        if(k_row_per > basis.rangeJ(j_row).lastIndex()){
                            k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                        }
                        if (overlap(psi_row.support(j_row,k_row_per), phi_col.support(j,k)) > 0) {
                            ret.insert(Index1D(j_row,k_row_per,XWavelet));
                        }
                    }
                }
            }
        }

        else {
            Support<T> supp = psi_col.support(j,k);

            // Inserting all indices corresponding to Bsplines with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            if (fabs(j - jmin) <= s_tilde) {
                int kMin = floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2);
                int kMax =  ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1);
                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    int k_row_per = k_row;
                    if(k_row_per < basis.mra.rangeI(jmin).firstIndex()){
                        k_row_per = basis.mra.rangeI(jmin).lastIndex() + ((1 - (basis.mra.rangeI(jmin).firstIndex() - k_row_per))%basis.mra.cardI(jmin));
                    }
                    if(k_row_per > basis.mra.rangeI(jmin).lastIndex()){
                        k_row_per = basis.mra.rangeI(jmin).firstIndex() - ((1 - (k_row_per - basis.mra.rangeI(jmin).lastIndex()))%basis.mra.cardI(jmin));
                    }

                    if (overlap(phi_row.support(jmin,k_row_per), psi_col.support(j,k)) > 0) {
                        ret.insert(Index1D(jmin,k_row_per,XBSpline));
                    }
                }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = psi_col.singularSupport(j,k);
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            int k_row_per = k_row;
                            if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                                k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                            }
                            if(k_row_per> basis.rangeJ(j_row).lastIndex()){
                                k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                            }
                            Support<T> supp_row = psi_row.support(j_row,k_row_per);
                            if (((overlap(supp_row, psi_col.support(j,k)) > 0)) &&
                                (!(distance(psi_col.singularSupport(j,k),supp_row) >= 0 )) ) {
                                    ret.insert(Index1D(j_row,k_row_per,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        int k_row_per = k_row;
                        if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                            k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                        }
                        if(k_row_per > basis.rangeJ(j_row).lastIndex()){
                            k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                        }
                        if ((overlap(psi_row.support(j_row,k_row_per), psi_col.support(j,k)) > 0) &&
                            !(distance(psi_row.singularSupport(j_row,k_row_per),supp) >= 0 ) ) {
                           ret.insert(Index1D(j_row,k_row_per,XWavelet));
                        }
                    }
                }
            }
        }
    }
    return ret;
}

template <typename T, Construction Cons>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Interval,Cons> &basis,
                  int s_tilde, int jmin, int jmax, bool /*update*/) 
{
        using std::min;
        using std::max;

        BSpline<T,Primal,Interval,Cons> phi_col(basis.mra), phi_row(basis.mra);
        Wavelet<T,Primal,Interval,Cons> psi_col(basis), psi_row(basis);
        int j = lambda.j, k = lambda.k;
        IndexSet<Index1D> ret;

        if (lambda.xtype==XBSpline) {

            Support<T> supp_col = phi_col.support(j,k);

            //Adding B-Splines
            int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
            int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(jmin)),kMin), kMax);
            //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
            while ((kStart-1 >= kMin) && (overlap(supp_col, phi_row.support(jmin,max(kStart-1, kMin)))>0)) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)),kMax), kMin);
            assert((overlap(supp_col, phi_col.support(jmin,kEnd))>0));
            while ((kEnd+1 <= kMax) && (overlap(supp_col, phi_row.support(jmin,min(kEnd+1,kMax)))>0)) {
                ++kEnd;
            }

            for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                ret.insert(Index1D(jmin,k_row,XBSpline));
            }

            //Adding Wavelets
            for (int j_row=jmin; j_row<=min(jmin+s_tilde, jmax); ++j_row) {

                int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
                int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
                //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
                while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
                    --kStart;
                }
                int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
                //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
                    ++kEnd;
                }
                for (int k_row=kStart; k_row<=kEnd; k_row++) {
                    Range<int> rangeL = basis.rangeJL(j_row);
                    Range<int> rangeR = basis.rangeJR(j_row);
                    if ( (k_row <= rangeL.lastIndex()) || (k_row >= rangeR.firstIndex()) ) {
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                        continue;
                    }
                    if  (!(distance(phi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) >= 0 )) {    //singsupp!
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                }
            }
        }
        else {
            Support<T> supp_col = psi_col.support(j,k);

            //Adding B-Splines
            if (fabs(j - jmin) <= s_tilde) {
                int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
                int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(jmin)), kMin), kMax);
                //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
                while (kStart-1>=kMin && overlap(supp_col,phi_row.support(jmin,max(kStart-1,kMin)))>0) {
                    --kStart;
                }
                int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)), kMax), kMin);
                //assert((overlap(supp_col, phi_row.support(jmin,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp_col,phi_row.support(jmin,min(kEnd+1,kMax)))>0) {
                    ++kEnd;
                }

                for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                    if  (distance(psi_col.singularSupport(j,k),phi_row.support(jmin,k_row)) < 0 ) {        //singsupp!
                        ret.insert(Index1D(jmin,k_row,XBSpline));
                    }
                    else {
                        if ( (k <= basis.rangeJL(j).lastIndex() )  ||
                             (k >= basis.rangeJR(j).firstIndex() )     ) {
                                ret.insert(Index1D(jmin,k_row,XBSpline));
                        }
                    }
                }
            }

            //Adding Wavelets
            for (int j_row=max(j-s_tilde,jmin); j_row<=min(j+s_tilde,jmax); ++j_row) {

                int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
                int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
                //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
                while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
                    --kStart;
                }
                int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
                //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
                    ++kEnd;
                }
                for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                    if (distance(psi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) < 0 ) {
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                        continue;
                    }
                    else {
                        if ( (k_row <= basis.rangeJL(j_row).lastIndex() )  ||
                             (k_row >= basis.rangeJR(j_row).firstIndex() )     ) {
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                        continue;
                    }

                    if (distance(psi_row.singularSupport(j_row,k_row),supp_col) < 0 ) {
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                    else {
                        if ( (k <= basis.rangeJL(j).lastIndex() )  ||
                             (k >= basis.rangeJR(j).firstIndex() )     ) {
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                        }
                    }

                    /*
                    if  ( (!(distance(psi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) >= 0 )) &&    //singsupp!
                          (!(distance(psi_row.singularSupport(j_row,k_row), supp_col) >= 0 )) ) {            //singsupp!
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                    */
                }
            }
        }
        return ret;
}

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE_WO_XBSpline(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis, 
                              int s_tilde, int jmin, int jmax)
{
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;
    int j = lambda.j, k = lambda.k;
    IndexSet<Index1D> ret;
    Support<T> support_refwavelet = psi.support(0,0);
    Support<T> supp = psi.support(j,k);

    // Inserting all indices corresponding to Wavelets with intersecting support using
    // a) local compactness  b) matrix compression  c) vanishing moments
    for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
        T Pow2i_Mjrow = pow2i<T>(-j_row);
        if (j_row>=j+2) {
            DenseVector<Array<T> > singsupp = psi.optim_singularSupport(j,k);
            for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                    Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));
                    if ((overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) >= 0 ))){
                        ret.insert(Index1D(j_row,k_row,XWavelet));
                    }
                }
            }
        }
        else {
            int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
            int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                if (T(Pow2i_Mjrow*(support_refwavelet.l1+k_row)) >= T(Pow2i_Mjrow*(support_refwavelet.l2+k_row))) {
                    std::cout << "Warning3: Translation indices too large!!" << std::endl;
                    continue;
                }
                Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));
                if ((overlap(supp, supp_row) > 0) && (!(distance(psi.optim_singularSupport(j_row,k_row),supp) > 0 ))) {
                    ret.insert(Index1D(j_row,k_row,XWavelet));
                }
            }
        }
    }
    return ret;
}

} // namespace lawa

