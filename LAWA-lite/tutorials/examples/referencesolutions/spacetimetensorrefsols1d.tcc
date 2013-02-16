namespace lawa {

template < typename T, typename Basis2D>
int
SpaceTimeTensorRefSols1D<T, Basis2D>::nr;

template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::c;

template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::k;

template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::r;

template < typename T, typename Basis2D>
flens::DenseVector<flens::Array<T> >
SpaceTimeTensorRefSols1D<T, Basis2D>::sing_pts_t;

template < typename T, typename Basis2D>
flens::DenseVector<flens::Array<T> >
SpaceTimeTensorRefSols1D<T, Basis2D>::sing_pts_x;

template < typename T, typename Basis2D>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
SpaceTimeTensorRefSols1D<T, Basis2D>::deltas_t;

template < typename T, typename Basis2D>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
SpaceTimeTensorRefSols1D<T, Basis2D>::deltas_x;



    
template < typename T, typename Basis2D>
void
SpaceTimeTensorRefSols1D<T, Basis2D>::setExample(int ex, T _c, T _k, T _r)
{
    nr = ex;
    c = _c;
    k = _k;
    r = _r;
    
    if (nr <= 3) {
        if(!(Basis2D::FirstBasisType::Domain == Periodic)) {
            std::cerr << "First Basis Type must be Periodic" << std::endl;
        }
        if(!(Basis2D::SecondBasisType::Domain == Interval)) {
            std::cerr << "Second Basis Type must be Interval" << std::endl;
        }
    }
    
    switch (nr) {
        case 1:
            sing_pts_t.engine().resize(2);
            sing_pts_t = 0., 1.;
            sing_pts_x.engine().resize(2);
            sing_pts_x = 0., 1.;
            // no Peaks
            break;
        case 2: 
            sing_pts_t.engine().resize(3);
            sing_pts_t = 0., 0.5, 1.;
            sing_pts_x.engine().resize(2);
            sing_pts_x = 0., 1.;
            // no Peaks
            break;
        case 3: 
            sing_pts_t.engine().resize(3);
            sing_pts_t = 0., 2.-std::sqrt(2), 1.;
            sing_pts_x.engine().resize(2);
            sing_pts_x = 0., 1.;
            // no Peaks
        default:
            std::cerr << "Example does not exist!" << std::endl;
            exit(1);
    }
}
    
template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::sol_t(T t, int deriv_t)
{   
    switch (nr) {
        case 1:
            if (deriv_t==0){    return          std::cos(2*M_PI*t);}
            if (deriv_t==1){    return     -2*M_PI*std::sin(2*M_PI*t);}
            break;
        case 2:
            if (deriv_t==0){
                if(t < 0.5) return  t + 0.25;
                else        return -t + 1.25;
            }
            if (deriv_t==1){
                if(t < 0.5) return  1;
                else        return -1;
            }
            break;
        case 3:
            if (deriv_t==0){
                if(t < 2. - std::sqrt(2.))  return  t*t + 0.5;
                else                        return  2*(t-1)*(t-1) + 0.5;
            }
            if (deriv_t==1){
                if(t < 2.-std::sqrt(2))     return  2*t;
                else                        return  4*(t-1);
            }
            break;
        default: 
            std::cerr << "Example does not exist!" << std::endl;
            exit(1); 
    }
    return 0;
}    

template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::sol_x(T x, int deriv_x){
    switch(nr){
        case 1:
            if (deriv_x == 0)   return  -4*(x-0.5)*(x-0.5)+1;
            if (deriv_x == 1)   return    -8*(x-0.5);
            if (deriv_x == 2)   return  -8;
            break;
        case 2:
            if (deriv_x == 0)   return  8*std::pow(x-0.5, 3) - 2*x*x +1;
            if (deriv_x == 1)   return    24*(x-0.5)*(x-0.5) -4*x;
            if (deriv_x == 2)   return  48*(x-0.5) - 4;
            break;
        case 3:
            if (deriv_x == 0)   return  8*std::pow(x-0.5, 3) - 2*x*x +1;
            if (deriv_x == 1)   return    24*(x-0.5)*(x-0.5) -4*x;
            if (deriv_x == 2)   return  48*(x-0.5) - 4;
            break;
        default: 
            std::cerr << "Example does not exist!" << std::endl; 
            exit(1);
    }
    return 0;
}

template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::sol(T t, T x, int deriv_t, int deriv_x)
{
    return sol_t(t, deriv_t) * sol_x(x, deriv_x);
}

template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::sol(T t, T x)
{
    return sol_t(t, 0) * sol_x(x, 0);
}

template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::dx_sol(T t, T x)
{
    return sol_t(t,0) * sol_x(x ,1);
}
    
template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::u_t(T t)
{
    return sol_t(t,0);
}    
    
template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::u_x(T x)
{
    return sol_x(x,0);
}

// rhs_contrib_t = d/dt u_t - 0.5 * r * u_x
template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::rhs_contrib_t(T t)
{
    return sol_t(t, 1) + 0.5 * r * u_t(t);
}
    
// rhs_contrib_x = - c * d/dxx u_x + k * d/dx u_x + 0.5 * r * u_x
template < typename T, typename Basis2D>
T
SpaceTimeTensorRefSols1D<T, Basis2D>::rhs_contrib_x(T x)
{
    return - c * sol_x(x, 2) + k * sol_x(x,1) + 0.5 * r * u_x(x);
}

} // namespace lawa
