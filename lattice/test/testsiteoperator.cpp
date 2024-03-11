// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/test/testsiteoperator.cpp
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

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"

using namespace QuantumNumbers;

int main()
{
   SymmetryList Symmetry("N:U(1),S:SU(2)");
   // helper function object for creating quantum numbers at compile time
   QNConstructor<U1,SU2> QN(Symmetry);

   SiteBasis HubbardBasis(Symmetry);
   HubbardBasis.push_back("Empty",  QN(0, 0));
   HubbardBasis.push_back("Single", QN(1, 0.5));
   HubbardBasis.push_back("Double", QN(2, 0));

   // Define the annihilation operator
   SiteOperator C(HubbardBasis, QN(-1, 0.5));
   C("Empty",  "Single")    = sqrt(2.0);
   C("Single", "Double")    = 1.0;

   // and the creation operator
   SiteOperator CH(HubbardBasis, QN(1, 0.5));
   CH("Single", "Empty")     = 1.0;
   CH("Double", "Single")    = -sqrt(2.0);

   // check that CH = adjoint(C)
   CHECK(equal(adjoint(C), CH))(adjoint(C))(CH);

   // check that C * adoint(C) is what we expect (ie. N-2)
   SiteOperator CCH(HubbardBasis, QN(0, 0));
   CCH("Empty",  "Empty")  = -2.0;
   CCH("Single", "Single") = -1.0;

   CHECK(equal(CCH, sqrt(2.0) * prod(C, adjoint(C), QN(0,0))))(CCH)(sqrt(2.0)
                                                                    * prod(C, adjoint(C), QN(0,0)));

   // check that adjoint(C) * C is what we expect (ie. N)
   SiteOperator N(HubbardBasis, QN(0, 0));
   N("Single", "Single")  = 1.0;
   N("Double", "Double")  = 2.0;

   CHECK(equal(N, sqrt(2.0) * prod(adjoint(C), C, QN(0,0))))(N)(sqrt(2.0)
                                                                * prod(adjoint(C), C, QN(0,0)));

   // the 'inverse adjoint' should give the opposite sign for these spin-1/2 operators.

   CHECK(equal(inv_adjoint(C), -CH))(inv_adjoint(C))(-CH);
   CHECK(equal(CCH, -sqrt(2.0)
               * prod(C, inv_adjoint(C), QN(0,0))))(CCH)(-sqrt(2.0)
                                                         * prod(C, inv_adjoint(C), QN(0,0)));

   CHECK(equal(N, -sqrt(2.0)
               * prod(inv_adjoint(C), C, QN(0,0))))(N)(-sqrt(2.0)
                                                       * prod(inv_adjoint(C), C, QN(0,0)));

   TRACE(N);
   TRACE(CH);
   show_projections(std::cout, CH);


   SiteProductBasis SPB(HubbardBasis, HubbardBasis);
   SiteOperator CpCH = tensor_prod(C, CH, SPB, QN(0,1));
   SiteOperator CpCHCheck(CpCH.Basis(), CpCH.TransformsAs());
   std::vector<std::pair<SiteOperator, SiteOperator> > Decomp = decompose_tensor_prod(CpCH, SPB);
   for (std::size_t i = 0; i < Decomp.size(); ++i)
   {
      CpCHCheck += tensor_prod(Decomp[i].first, Decomp[i].second, SPB, CpCH.TransformsAs());
   }
   CHECK_CLOSE(CpCH, CpCHCheck);
}
