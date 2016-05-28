// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/old/factorized-u1u1.h
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
//
// factorized-u1u1.h
//
// Site blocks for bosonic and fermionic models factorized into 
// (0 up) \otimes (0 down) sites.
//

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"




inline
LatticeSite FactorizedSite(int MyN, half_int MySz, 
                         std::string const& Sym1 = "N", 
                         std::string const& Sym2 = "Sz")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CH, C, P, R, N, I, Sz;
   LatticeSite Site;

   Basis.push_back("0", QN(0,0));
   Basis.push_back("1", QN(MyN,MySz));

   CH = SiteOperator(Basis, QN(MyN,MySz), LatticeCommute(LatticeCommute::Fermionic) * MyN);
   C = SiteOperator(Basis, QN(-MyN,-MySz), LatticeCommute(LatticeCommute::Fermionic) * MyN);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // annihilate particle
   C("0", "1") = std::sqrt(double(MyN));
   
   // create particle
   CH = adjoint(C);

   // identity
   I("0", "0") = 1;
   I("1", "1") = 1;

   // parity = (-1)^N
   P("0", "0") = 1;
   P("1", "1") = minus1pow(MyN);

   R = I;

   // particle number  
   N = prod(CH, C, QN(0,0));

   // z-component of spin
   Sz("1", "1") = MySz.to_double();

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["N"] = N;
   Site["Sz"] = Sz;
   Site["CH"] = CH;
   Site["C"] = C;
   return Site;
}
