#include <cassert>
#include <iostream>

namespace lawa {
  
template <typename T>
    T
    _linear_wavelet_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_evaluator2(T x, unsigned short deriv);

//------------------------------------------------------------------------------

template <typename T>
Wavelet<T,Orthogonal,R,Multi>::Wavelet(int _d)
    : d(_d), vanishingMoments(_d)
{
    assert(d>=2);
    
    switch (d) {
        case 2: _numSplines = 3;

                _evaluator = new Evaluator[3];
                _evaluator[0] = _linear_wavelet_evaluator0;
                _evaluator[1] = _linear_wavelet_evaluator1;                
                _evaluator[2] = _linear_wavelet_evaluator2;
            
                _support = new Support<T>[3];
                _support[0] = Support<T>(-1,1);
                _support[1] = Support<T>(-1,1);
                _support[2] = Support<T>( 0,1);
            
                _singularSupport = new DenseVector<Array<T> >[3];
                _singularSupport[0].engine().resize(13,0);
                _singularSupport[0] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,
                                       0.125,0.25,0.375,0.5,0.75,1.0;
                _singularSupport[1].engine().resize(13,0);
                _singularSupport[1] = -1.0,-0.75,-0.5,-0.375,-0.25,-0.125,0.0,
                                       0.125,0.25,0.375,0.5,0.75,1.0;
                _singularSupport[2] = linspace(0.0,1.0,9);
                break;

        default: std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
}
    
template <typename T>
Wavelet<T,Orthogonal,R,Multi>::~Wavelet()
{
    delete[] _evaluator;
    delete[] _support;
    delete[] _singularSupport;
}

template <typename T>
T
Wavelet<T,Orthogonal,R,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2ih<T>(2*j*deriv+j) *
    _evaluator[type](pow2i<T>(j)*x - shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Orthogonal,R,Multi>::support(int j, int k) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2i<T>(-j) * (_support[type] + shift);    
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Orthogonal,R,Multi>::singularSupport(int j, int k) const
{
    const int typ = _type(k);
    const long shift = _shift(k);
    
    DenseVector<Array<T> > result = _singularSupport[typ];
    result += shift;
    
    return pow2i<T>(-j) * result;
}

template <typename T>
long
Wavelet<T,Orthogonal,R,Multi>::_shift(long k) const
{
    return k>=0 ? k/_numSplines : -((-k-1)/_numSplines+1);
}

template <typename T>
int
Wavelet<T,Orthogonal,R,Multi>::_type(long k) const
{
    return k>=0 ? (int) k%3 : (int) _numSplines - (int)(-k+2)%_numSplines - 1;
}

//------------------------------------------------------------------------------
    
template <typename T>
T
_linear_wavelet_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0 <= x && x < 0.125) {
            value = 2.17124059336723766167 - 33.0125730435099442245 * x;
        } else if(-0.375 <= x && x < -0.25) {
            value = -7.14234036204314274860 - 19.69261569230427143433 * x;
        } else if(0.25 <= x && x < 0.375){
            value = 2.27349781934085223456 - 5.6124494201652150829 * x;
        } else if(-0.5 <= x && x < -0.375){
            value = -1.24277494216897627173 - 3.96044123930649416268 * x;
        } else if(0.375 <= x && x < 0.5 ){
            value = 1.11118460317702247406 - 2.51294751039500238824 * x;
        } else if(-1.0 <= x && x < -0.75){
            value = -0.491630451656180539739 - 0.491630451656180539739 * x;
        } else if(0.75 <= x && x < 1){
            value = 0.096859434680319146709 - 0.096859434680319146709 * x;
        } else if(0.5 <= x && x < 0.75) {
            value = -0.48429717340159573354 + 0.67801604276223402696 * x;
        } else if(-0.75 <= x && x < -0.5){
            value = 2.458152258280902698695 + 3.44141316159326377817 *x;
        } else if(-0.125 <= x && x < 0){
            value = 2.17124059336723766167 + 17.1233461124244526904 *x;
        } else if(-0.25 <= x && x < -0.125){
            value = 2.28083109759543704076 + 18.0000701462500477231 *x;
        } else if(0.125 <= x && x < 0.25){
            value = -4.78104753844255919661 + 22.6057320109684306418* x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0 <= x && x < 0.125) {
            value = - 33.0125730435099442245 ;
        } else if(-0.375 <= x && x < -0.25) {
            value = - 19.69261569230427143433 ;
        } else if(0.25 <= x && x < 0.375){
            value = -5.6124494201652150829 ;
        } else if(-0.5 <= x && x < -0.375){
            value = -3.96044123930649416268 ;
        } else if(0.375 <= x && x < 0.5 ){
            value = -2.51294751039500238824 ;
        } else if(-1.0 <= x && x < -0.75){
            value = -0.491630451656180539739 ;
        } else if(0.75 <= x && x < 1){
            value = -0.096859434680319146709 ;
        } else if(0.5 <= x && x < 0.75) {
            value =  0.67801604276223402696 ;
        } else if(-0.75 <= x && x <-0.5){
            value =  3.44141316159326377817 ;
        } else if(-0.125 <= x && x < 0){
            value = 17.1233461124244526904 ;
        } else if(-0.25 <= x && x < -0.125){
            value = 18.0000701462500477231 ;
        } else if(0.125 <= x && x < 0.25){
            value = 22.6057320109684306418;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_evaluator1(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0.125 <= x && x< 0.25) {
            value = 6.84500149587779311472 - 33.09125988795909139990* x;
        } else if(-0.375 <= x && x < -0.25){
            value = -4.757606755902492947821 - 13.01056641923786902976 *x;
        } else if(-0.5 <= x && x < -0.375){
            value = -1.761589978102878832866 - 5.02118834510556472321*x;
        } else if(0.5 <= x && x < 0.75){
            value = 0.953647148545425993185 - 1.335106007963596390459*x;
        } else if(-1.0 <= x && x < -0.75){
            value = -0.4993361296332690191592 - 0.4993361296332690191592*x;
        } else if(0.75 <= x && x < 1.0){
            value = -0.1907294297090851986370 + 0.1907294297090851986370*x;
        } else if(-0.125 <= x && x < 0.0){
            value = 5.797244214189106508582*x;
        } else if(0.0 <= x && x < 0.125){
            value = 21.66875207906325351789* x;
        } else if(-0.75 <= x && x < -0.5){
            value = 2.496680648166345095796 + 3.495352907432883134114*x;
        } else if(0.375 <= x && x < 0.5){
            value = -2.378803377951246473911 + 5.32979504502974854373*x;
        } else if(-0.25 <= x && x < -0.125){
            value = 0.0556540975457490632366 + 6.24247699455509901447*x;
        } else if(0.25 <= x && x < 0.375){
            value = -3.523179956205757665732 + 8.38146592037511172192*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.125 <= x && x< 0.25) {
            value = - 33.09125988795909139990;
        } else if(-0.375 <= x && x < -0.25){
            value =  - 13.01056641923786902976 ;
        } else if(-0.5 <= x && x < -0.375){
            value =  - 5.02118834510556472321;
        } else if(0.5 <= x && x < 0.75){
            value = - 1.335106007963596390459;
        } else if(-1.0 <= x && x < -0.75){
            value = - 0.4993361296332690191592;
        } else if(0.75 <= x && x < 1.0){
            value = 0.1907294297090851986370;
        } else if(-0.125 <= x && x <= 0.0){//TODO: //ttttttttttest <=0.0
            value = 5.797244214189106508582;
        } else if(0.0 <= x && x < 0.125){
            value = 21.66875207906325351789;
        } else if(-0.75 <= x && x < -0.5){
            value =  3.495352907432883134114;
        } else if(0.375 <= x && x < 0.5){
            value =  5.32979504502974854373;
        } else if(-0.25 <= x && x < -0.125){
            value =  6.24247699455509901447;
        } else if(0.25 <= x && x < 0.375){
            value =  8.38146592037511172192;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_wavelet_evaluator2(T x, unsigned short deriv)
{
    double value = 0.0;
    if (deriv == 0) {
        if (0.625 <= x && x < 0.75) {
            value = 25.3256618915819297146 - 35.2183939040239842910*x;
        } else if(0.25 <= x && x < 0.375){
            value = 3.007148855493370156849 - 9.72312653910588252146*x;
        } else if(0.375 <= x && x < 0.5){
            value = 2.887238244145711420071 - 9.40336490884545922339*x;
        } else if(0.0 <= x && x < 0.125){
            value = -0.56347593702122728620*x;
        } else if (0.875 <= x && x < 1.0){
            value = -2.22797669417418531629 + 2.22797669417418531629 *x;
        } else if (0.125 <= x && x < 0.25){
            value = -0.717236204972206348034 + 5.17441370275642349806*x;
        } else if (0.75 <= x && x < 0.875){
            value = -5.94595223442177053803 + 6.47709159731428271256*x;
        } else if (0.5 <= x && x < 0.625){
            value = -22.32888385765284908917 + 41.0288792947516617951*x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1) {
        if (0.625 <= x && x < 0.75) {
            value =  - 35.2183939040239842910;
        } else if(0.25 <= x && x < 0.375){
            value = - 9.72312653910588252146;
        } else if(0.375 <= x && x < 0.5){
            value = - 9.40336490884545922339;
        } else if(0.0 <= x && x < 0.125){
            value = -0.56347593702122728620;
        } else if (0.875 <= x && x < 1.0){
            value = 2.22797669417418531629 ;
        } else if (0.125 <= x && x < 0.25){
            value = 5.17441370275642349806;
        } else if (0.75 <= x && x < 0.875){
            value = 6.47709159731428271256;
        } else if (0.5 <= x && x < 0.625){
            value = 41.0288792947516617951;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

} // namespace lawa
