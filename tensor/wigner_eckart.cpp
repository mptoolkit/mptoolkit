// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/wigner_eckart.cpp
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

#include "wigner_eckart.h"

using namespace QuantumNumbers;

namespace Tensor
{

WignerEckartBasis<VectorBasis>::
WignerEckartBasis(VectorBasis const& B,
                  QPSetType& AllowedProjections,
                  SymmetryList const& FinalSL)
   : NonAbBasis(B),
     AbBasis(FinalSL),
     Mapping(B.size())
{
   QPSetType UsedProjections;
   for (unsigned i = 0; i < NonAbBasis.size(); ++i)
   {
      ProjectionList p = enumerate_projections(NonAbBasis[i]);
      for (unsigned pi = 0; pi < p.size(); ++pi)
      {
         if (AllowedProjections.count(std::make_pair(NonAbBasis[i], p[pi])))
         {
            UsedProjections.insert(std::make_pair(NonAbBasis[i], p[pi]));
            Mapping[i][pi] = AbBasis.size();
            AbBasis.push_back(map_projection_to_quantum(p[pi], FinalSL), NonAbBasis.dim(i));
         }
      }
   }
   AllowedProjections = UsedProjections;
}

WignerEckartBasis<VectorBasis>::
WignerEckartBasis(VectorBasis const& B,
                  SymmetryList const& FinalSL)
   : NonAbBasis(B),
     AbBasis(FinalSL),
     Mapping(B.size())
{
   QPSetType UsedProjections;
   for (unsigned i = 0; i < NonAbBasis.size(); ++i)
   {
      ProjectionList p = enumerate_projections(NonAbBasis[i]);
      for (unsigned pi = 0; pi < p.size(); ++pi)
      {
         Mapping[i][pi] = AbBasis.size();
         AbBasis.push_back(map_projection_to_quantum(p[pi], FinalSL), NonAbBasis.dim(i));
      }
   }
}

void InsertAllowedProjections(QPSetType const& qp, QPSetType& Result, QuantumNumber const& q)
{
   ProjectionList pl = enumerate_projections(q);
   for (QPSetType::const_iterator I = qp.begin(); I != qp.end(); ++I)
   {
      QuantumNumberList ql = transform_targets(I->first, q);
      for (unsigned qi = 0; qi < ql.size(); ++qi)
      {
         for (unsigned pi = 0; pi < pl.size(); ++pi)
         {
            Projection p = sum(I->second, pl[pi]);
            if (is_possible(ql[qi], p))
               Result.insert(std::make_pair(ql[qi], p));
         }
      }
   }
}

QPSetType UpdateAllowedProjections(QPSetType const& qp, std::set<QuantumNumber> const& q)
{
   QPSetType Result;
   for (std::set<QuantumNumber>::const_iterator I = q.begin(); I != q.end(); ++I)
   {
      InsertAllowedProjections(qp, Result, *I);
   }
   return Result;
}

} // namespace Tensor
 
