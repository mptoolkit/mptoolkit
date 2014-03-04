// -*- C++ -*- $Id$

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

