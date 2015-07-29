// -*- C++ -*- $Id$

#include "lattice/latticesite.h"




LatticeSite SiteTensorProduct(LatticeSite const& A, LatticeSite const& B)
{
   SiteBasis ABasis = A.Basis1();
   SiteBasis BBasis = B.Basis1();
   SiteProductBasis Basis(ABasis, BBasis);

   LatticeSite Result;

   for (LatticeSite::const_iterator ai = A.begin(); ai != A.end(); ++ai)
   {
      for (LatticeSite::const_iterator bi = B.begin(); bi != B.end(); ++bi)
      {
	 if (!(ai->first == "I" || ai->first == "P" || bi->first == "I" || bi->first == "P"))
	    continue;

	 QuantumNumbers::QuantumNumberList ql = transform_targets(ai->second.TransformsAs(),
						  bi->second.TransformsAs());
	 for (std::size_t q = 0; q < ql.size(); ++q)
	 {
	    SiteOperator ProdOp = tensor_prod(ai->second, bi->second, Basis, ql[q]);
	    std::string Name = ai->first + "*" + bi->first;
	    if (ql.size() != 1)
	       Name = Name + '[' + boost::lexical_cast<std::string>(ql[q]) + ']';
	    if (Name == "I*I")
	       Result["I"] = ProdOp;  // make sure we have the identity operator
	    else if (Name == "P*P") 
	       Result["P"] = ProdOp;  // make sure we have the parity operator
	    Result[Name] = ProdOp;
	 }
      }
   }
   return Result;
}
