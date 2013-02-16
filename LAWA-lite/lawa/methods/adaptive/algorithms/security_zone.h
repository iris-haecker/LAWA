#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_SECURITY_ZONE_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_SECURITY_ZONE_H 1

#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa {

//Security zone for an 1d-index following Urban:2009, p.235 and KU:2010.
template <typename T, DomainType Domain, Construction Cons>
    IndexSet<Index1D>
    C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,Domain,Cons> &basis);

template <typename T, DomainType Domain, Construction Cons>
    IndexSet<Index1D>
    C(const Index1D &lambda, T c, const Basis<T,Primal,Domain,Cons> &basis);

template <typename T>
    IndexSet<Index1D>
    C_WO_XBSpline(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,R,CDF> &basis, bool only_pos=false);

template <typename T>
    IndexSet<Index1D>
    C_WO_XBSpline(const Index1D &lambda, T c, const Basis<T,Primal,R,CDF> &basis);


template <typename T, DomainType Domain, Construction Cons>
    void
    C(const Index1D &lambda, T c, const MRA<T,Primal,Domain,Cons> &mra,
      const Basis<T,Primal,Domain,Cons> &basis, IndexSet<Index1D> &ret);

// Computation of a security zone for 2d-tensor basis
template <typename T, typename Basis2D>
    IndexSet<Index2D>
    C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis);

template <typename T, typename Basis2D>
    IndexSet<Index2D>
    C_t(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis);

template <typename T, typename Basis3D>
    IndexSet<Index3D>
    C(const IndexSet<Index3D> &Lambda, T c, const Basis3D &basis);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/security_zone.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_SECURITY_ZONE_H

