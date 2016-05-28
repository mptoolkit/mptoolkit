// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensorproduct.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

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

ProductBasis<BasisList, BasisList>
ProductBasis<BasisList, BasisList>::MakeTriangularProjected(left_basis_type Basis1_, right_basis_type Basis2_, 
							    QuantumNumbers::QuantumNumber const& q)
{
   return ProductBasis<BasisList, BasisList>(base_type::MakeTriangularProjected(Basis1_, Basis2_, q));
}

} // namespace Tensor
