// -*- C++ -*- $Id$

#if !defined(WIGNER_ECKART_H_KSDJHCUIREYT7YO7RE)
#define WIGNER_ECKART_H_KSDJHCUIREYT7YO7RE

#include "tensor.h"

namespace Tensor
{

template <typename BasisType>
class WignerEckartBasis;

// a container for representing the allowed set of projections
typedef std::set<std::pair<QuantumNumber, Projection> > QPSetType;

template <>
class WignerEckartBasis<VectorBasis>
{
   public:
      typedef VectorBasis basis_type;

      // constructs a WignerEckartBasis from a given VectorBasis and a set of allowed projections.
      // This also updates the AllowedProjections set so that pairs that don't actually occur
      // in the final basis are removed.
      WignerEckartBasis(VectorBasis const& B,
                        QPSetType& AllowedProjections,
                        SymmetryList const& FinalSL);

      // constructs a WignerEckartBasis from a given VectorBasis, taking all possible projections
      WignerEckartBasis(VectorBasis const& B,
                        SymmetryList const& FinalSL);

      basis_type const& NonAbelianBasis() const { return NonAbBasis; }
      basis_type const& AbelianBasis() const { return AbBasis; }

      // returns the subspace in the abelian basis corresponding to the 
      // given non-abelian subspace and the projection, or -1 if this state is not in the basis
      int Map(int NonAbelianSubspace, int Projection) const;

   private:
      VectorBasis NonAbBasis;
      VectorBasis AbBasis;
      // This maps the NonAbBasis subspace at a given projection,
      // onto the AbBasis subspace.  We don't attempt to merge identical quantum numbers
      // in the AbBasis here.
      std::vector<std::map<int, int> > Mapping;
};

// update the allowed projections by taking all possible combinations
// of the existing projections with the quantum number q
QPSetType UpdateAllowedProjections(QPSetType const& qp, std::set<QuantumNumber> const& q);

// given an irreducible tensor, project onto the given Projection p, using
// the Wigner Eckart mappings b1 and b2, for the left basis and right basis respectively.
template <typename T, typename B1T, typename B2T, typename S>
IrredTensor<T, B1T, B2T, S>
wigner_eckart(IrredTensor<T, B1T, B2T, S> const& x, 
              Projection const& p, 
              WignerEckartBasis<B1T> const& b1, 
              WignerEckartBasis<B1T> const& b2);

} // namespace Tensor

#include "wigner_eckart.cc"

#endif
