#ifndef IRIS2_MY_STATICRHS_TCC
#define IRIS2_MY_STATICRHS_TCC 1

namespace lawa {

template <typename T, typename Index, typename Preconditioner>
StaticRHS<T,Index,Preconditioner>::StaticRHS(const Preconditioner &_P,
                                             const Coefficients<Lexicographical,T,Index> &_rhs_data)
    :    P(_P), rhs_data(), rhs_abs_data()
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;

    for (const_coeff_it it=_rhs_data.begin(); it!=_rhs_data.end(); ++it) {
        rhs_data[(*it).first] = (*it).second;
        rhs_abs_data.insert(val_type((*it).second, (*it).first));
    }
}

template <typename T, typename Index, typename Preconditioner>
T
StaticRHS<T,Index,Preconditioner>::operator()(const Index &lambda)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;
    const_coeff_it it_end       = rhs_data.end();
    const_coeff_it it_index     = rhs_data.find(lambda);

    if (it_index != it_end) {
        return (*it_index).second;
    }
    
    //std::cerr << "Index not found: lambda = " << lambda << std::endl;
    //assert(0);

    return 0;
}

template <typename T, typename Index, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
StaticRHS<T,Index,Preconditioner>::operator()(const IndexSet<Index> &Lambda)
{
    typedef typename IndexSet<Index>::iterator const_it;
    Coefficients<Lexicographical,T,Index> ret;
    for (const_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        T tmp = operator()(*lambda);
        ret[*lambda] = tmp;
    }
    return ret;
}

template <typename T, typename Index, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
StaticRHS<T,Index,Preconditioner>::operator()(T tol)
{
    Coefficients<Lexicographical,T,Index> ret;
    ret = THRESH(rhs_abs_data,tol);
    return ret;
}

}    //namespace lawa

#endif // IRIS2_MY_STATICRHS_TCC