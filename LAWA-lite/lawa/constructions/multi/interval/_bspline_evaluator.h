#ifndef LAWA_CONSTRUCTIONS_MULTI_INTERVAL__BSPLINE_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_MULTI_INTERVAL__BSPLINE_EVALUATOR_H 1

namespace lawa {
    
//--- linear evaluators --------------------------------------------------------
    
template <typename T>
    T
    _linear_bspline_left_evaluator0(T x, unsigned short deriv);
    
template <typename T>
    T
    _linear_bspline_inner_evaluator0(T x, unsigned short deriv);
    
template <typename T>
    T
    _linear_bspline_inner_evaluator1(T x, unsigned short deriv);
    
template <typename T>
    T
    _linear_bspline_inner_evaluator2(T x, unsigned short deriv);
    
template <typename T>
    T
    _linear_bspline_right_evaluator0(T x, unsigned short deriv);
        
//--- quadratic evaluators -----------------------------------------------------
    
template <typename T>
    T
    _quadratic_bspline_left_evaluator0(T x, unsigned short deriv);
    
template <typename T>
    T
    _quadratic_bspline_left_evaluator1(T x, unsigned short deriv);
    
template <typename T>
    T
    _quadratic_bspline_inner_evaluator0(T x, unsigned short deriv);
    
template <typename T>
    T
    _quadratic_bspline_inner_evaluator1(T x, unsigned short deriv);
    
template <typename T>
    T
    _quadratic_bspline_inner_evaluator2(T x, unsigned short deriv);
    
template <typename T>
    T
    _quadratic_bspline_inner_evaluator3(T x, unsigned short deriv);
    
template <typename T>
    T
    _quadratic_bspline_inner_evaluator4(T x, unsigned short deriv);
    
template <typename T>
    T
    _quadratic_bspline_inner_evaluator5(T x, unsigned short deriv);
    
//--- quadratic evaluators -----------------------------------------------------

template <typename T>
    T
    _cubic_bspline_inner_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_bspline_inner_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_bspline_inner_evaluator2(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_bspline_inner_evaluator3(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_bspline_inner_evaluator4(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_bspline_inner_evaluator5(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_bspline_left_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_bspline_right_evaluator1(T x, unsigned short deriv);

} // namespace lawa

#include <lawa/constructions/multi/interval/_bspline_evaluator.tcc>

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL__BSPLINE_EVALUATOR_H
