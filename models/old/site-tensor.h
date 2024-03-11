// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/old/site-tensor.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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
