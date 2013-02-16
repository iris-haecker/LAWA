namespace lawa {

template <typename T, typename Index>
T
NoPreconditioner<T, Index>::operator()(const Index &/*index*/) const
{
    return 1.;
}

}   // namespace lawa

