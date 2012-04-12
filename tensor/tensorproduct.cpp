// -*- C++ -*- $Id$

#include "tensorproduct.h"

namespace Tensor
{

// ProductBasis<VectorBasis, VectorBasis>
ProductBasis<VectorBasis, VectorBasis>::ProductBasis(VectorBasis const& B1, VectorBasis const& B2)
   : base_type(B1, B2)
{
   for (std::size_t s = 0; s < Basis_.size(); ++s)
   {
      base_type::source_type src = ReverseMapping_[s];
      Basis_.set_dim(s, B1.dim(src.first) * B2.dim(src.second));
   }
}

ProductBasis<VectorBasis, VectorBasis>::ProductBasis(VectorBasis const& B1, VectorBasis const& B2,
                                                     QuantumNumber const& Target)
   : base_type(B1, B2, Target)
{
   for (std::size_t s = 0; s < Basis_.size(); ++s)
   {
      base_type::source_type src = ReverseMapping_[s];
      Basis_.set_dim(s, B1.dim(src.first) * B2.dim(src.second));
   }
}

// ProductBasis<BasisList, VectorBasis>
ProductBasis<BasisList, VectorBasis>::ProductBasis(BasisList const& B1, VectorBasis const& B2)
   : base_type(B1, B2)
{
   for (std::size_t s = 0; s < Basis_.size(); ++s)
   {
      base_type::source_type src = ReverseMapping_[s];
      Basis_.set_dim(s, B2.dim(src.second));
   }
}

ProductBasis<BasisList, VectorBasis>::ProductBasis(BasisList const& B1, VectorBasis const& B2,
                                                   QuantumNumber const& Target)
   : base_type(B1, B2, Target)
{
   for (std::size_t s = 0; s < Basis_.size(); ++s)
   {
      base_type::source_type src = ReverseMapping_[s];
      Basis_.set_dim(s, B2.dim(src.second));
   }
}

// ProductBasis<VectorBasis, VectorBasis>
ProductBasis<VectorBasis, BasisList>::ProductBasis(VectorBasis const& B1, BasisList const& B2)
   : base_type(B1, B2)
{
   for (std::size_t s = 0; s < Basis_.size(); ++s)
   {
      base_type::source_type src = ReverseMapping_[s];
      Basis_.set_dim(s, B1.dim(src.first));
   }
}

ProductBasis<VectorBasis, BasisList>::ProductBasis(VectorBasis const& B1, BasisList const& B2,
                                                   QuantumNumber const& Target)
   : base_type(B1, B2, Target)
{
   for (std::size_t s = 0; s < Basis_.size(); ++s)
   {
      base_type::source_type src = ReverseMapping_[s];
      Basis_.set_dim(s, B1.dim(src.first));
   }
}

} // namespace Tensor
