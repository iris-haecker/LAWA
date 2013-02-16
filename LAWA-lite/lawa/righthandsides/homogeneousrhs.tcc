namespace lawa {
    
template <typename T>
T 
HomogeneousRHS<T>::operator()(T /*time*/, XType /*xtype*/, int /*j*/, int /*k*/) const{
    return 0;
}

} // namespace lawa

