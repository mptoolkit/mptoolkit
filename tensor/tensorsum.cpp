// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensorsum.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "tensorsum.h"
#include <boost/type_traits.hpp>

namespace Tensor
{

void SumBasis<BasisList>::AddBasis(BasisList const& B)
{
   PRECONDITION(B.GetSymmetryList() == this->GetSymmetryList());
   BList_.push_back(B);
   std::vector<int> Map(B.size());
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Map[i] = this->size();
      this->push_back(B[i]);
   }
   Mapping_.push_back(Map);
}

void SumBasis<VectorBasis>::AddBasis(VectorBasis const& B)
{
   PRECONDITION(B.GetSymmetryList() == This_.GetSymmetryList());
   BList_.push_back(B);
   std::vector<int> Map(B.size());
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Map[i] = This_.size();
      This_.push_back(B[i], B.dim(i));
   }
   Mapping_.push_back(Map);
}

} // namespace Tensor
