#ifndef TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_SPACETIMETENSORREFSOLS1D_H
#define TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_SPACETIMETENSORREFSOLS1D_H 1

#include <vector>

namespace lawa {
    
template <typename T, typename Basis2D>
struct SpaceTimeTensorRefSols1D {
    
    typedef flens::DenseVector<flens::Array<T> >                       DenseVectorT;
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > FullColMatrixT;
    
    // Example Number
    static int   nr;
    
    // Parameters: diffusion, convection and reaction constant
    static T c; 
    static T k; 
    static T r;
    
    // Singular Points in each dimension (dim1 = Time, dim2: Space)
    static DenseVectorT sing_pts_t, sing_pts_x;
    // Delta-Peaks in each dimension     (dim1 = Time, dim2: Space)
    static FullColMatrixT deltas_t, deltas_x; 
    
    static void setExample(int ex, T _c, T _k, T _r);
    
    static T sol_t(T t, int deriv_t);
    static T sol_x(T x, int deriv_x);
    static T sol(T t, T x, int deriv_t, int deriv_x);
    static T sol(T t, T x);
    static T dx_sol(T t, T x);
    
    // Space-Time Diffusion-Convection-Reaction Examples
    // u(t,x) = u_t(t) * u_x(x)
    // d/dt u - c * d/dxx u + k * d/dx u + r * u 
    //          = (d/dt u_t - 0.5 * r * u_x) * u_x + u_t * ( - c * d/dxx u_x + k * d/dx u_x + 0.5 * r * u_x)
    //          =: rhs_contrib_t             * u_x + u_t *  rhs_contrib_x
    static T u_t(T t);
    static T u_x(T x);
    static T rhs_contrib_t(T t);
    static T rhs_contrib_x(T x);
};
    
    
} // namespace lawa

#include <tutorials/examples/referencesolutions/spacetimetensorrefsols1d.tcc>

#endif // TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_SPACETIMETENSORREFSOLS1D_H
