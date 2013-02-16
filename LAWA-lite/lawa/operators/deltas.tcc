namespace lawa {
/*
template <typename T, DomainType Domain, Construction Cons>
GeMatrix<FullStorage<T,ColMajor> >
computeDeltas(const BSpline<T,Primal,Domain,Cons> &phi, int j, int k)
{
    GeMatrix<FullStorage<T,ColMajor> > ret;

    ret.engine().resize(phi.singularSupport(j,k).length(),2);
    ret(_,1) = phi.singularSupport(j,k);
    Support<T> supp = phi.support(j,k);
    T step = pow2i<T>(-(j+5)); // 1.0/(1<<(j+1));
    for (int i = 1; i<=ret.numRows(); ++i) {
        ret(i,2) = ((ret(i,1)==supp.l2) ? 0.0 : phi(std::min(ret(i,1)+step, supp.l2),j,k,0))
                 - ((ret(i,1)==supp.l1) ? 0.0 : phi(std::max(ret(i,1)-step, supp.l1),j,k,0));
    }
    return ret;
}

template <typename T, DomainType Domain, Construction Cons>
GeMatrix<FullStorage<T,ColMajor> >
computeDeltas(const Wavelet<T,Primal,Domain,Cons> &psi, int j, int k)
{
    GeMatrix<FullStorage<T,ColMajor> > ret;

    ret.engine().resize(psi.singularSupport(j,k).length(),2);
    ret(_,1) = psi.singularSupport(j,k);
    Support<T> supp = psi.support(j,k);
    T step = pow2i<T>(-(j+5)); // 1.0/(1<<(j+1));
    for (int i = 1; i<=ret.numRows(); ++i) {
        ret(i,2) = ((ret(i,1)==supp.l2) ? 0.0 : psi(std::min(ret(i,1)+step, supp.l2),j,k,0))
                 - ((ret(i,1)==supp.l1) ? 0.0 : psi(std::max(ret(i,1)-step, supp.l1),j,k,0));
    }
    return ret;
}
*/

}    //namespace lawa
