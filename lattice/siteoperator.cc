// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/siteoperator.cc
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

inline
SiteOperator SiteOperator::Identity(SiteBasis const& Basis)
{
   SiteOperator Result(Basis, QuantumNumber(Basis.GetSymmetryList()), LatticeCommute::Bosonic);
   for (std::size_t i = 0; i < Basis.size(); ++i)
   {
      Result(i,i) = 1.0;
   }
   return Result;
}

inline
SiteOperator SiteOperator::Identity(SiteBasis const& Basis1, SiteBasis const& Basis2)
{
   CHECK_EQUAL(Basis1, Basis2);
   return SiteOperator::Identity(Basis1);
}

inline
SiteOperator
MakeIdentityFrom(SiteOperator const& x)
{
   CHECK_EQUAL(x.Basis1(), x.Basis2());
   return SiteOperator::Identity(x.Basis1());
}

inline
void CoerceSymmetryListInPlace(SiteOperator& s, SymmetryList const& sl)
{
   CoerceSymmetryListInPlace(s.Basis_, sl);
   CoerceSymmetryListInPlace(static_cast<IrredTensor<std::complex<double> >&>(s), sl);
}

inline
SiteOperator prod(SiteOperator const& x, SiteOperator const& y, QuantumNumber Trans)
{
   DEBUG_PRECONDITION_EQUAL(x.Basis(), y.Basis());
   return SiteOperator(x.Basis(), prod(x.base(), y.base(), Trans),
                       x.Commute() * y.Commute());
}

inline
SiteOperator adjoint(SiteOperator const& x)
{
   return SiteOperator(x.Basis(), adjoint(x.base()), x.Commute());
}

inline
SiteOperator inv_adjoint(SiteOperator const& x)
{
   return SiteOperator(x.Basis(), inv_adjoint(x.base()), x.Commute());
}

inline
SiteOperator
tensor_prod(SiteOperator const& S1, SiteOperator const& S2,
            SiteProductBasis const& SPBasis,
            QuantumNumber const& q)
{
   return SiteOperator(SPBasis.Basis(),
                       tensor_prod(S1.base(), S2.base(),
                                   SPBasis.PBasis(), SPBasis.PBasis(),q),
                       S1.Commute() * S2.Commute());
}

inline
SiteOperator
tensor_prod(SiteOperator const& S1, SiteOperator const& S2, QuantumNumber const& q)
{
   SiteProductBasis SP(S1.Basis1(), S2.Basis1());
   return tensor_prod(S1, S2, SP, q);
}

inline
SiteOperator
tensor_prod(SiteOperator const& S1, SiteOperator const& S2)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(S1.TransformsAs(), S2.TransformsAs());
   CHECK_EQUAL(ql.size(), 1);
   return tensor_prod(S1, S2, ql[0]);
}

inline
SiteOperator
cross(SiteOperator const& x, SiteOperator const& y)
{
   CHECK(cross_product_exists(x.TransformsAs(), y.TransformsAs()))
      ("Cross product does not exist for these operators")
      (x.TransformsAs())(y.TransformsAs());
   return cross_product_factor(x.TransformsAs(), y.TransformsAs())
      * prod(x, y, cross_product_transforms_as(x.TransformsAs(), y.TransformsAs()));
}
