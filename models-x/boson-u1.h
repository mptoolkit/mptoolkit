// -*- C++ -*- $Id$

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"




void SetMatElementU1(SiteOperator& s, int n1, int n2, double x)
{
   std::string q1 = boost::lexical_cast<std::string>(n1);
   std::string q2 = boost::lexical_cast<std::string>(n2);

   int l1 = s.Basis1().LookupOrNeg(q1);
   int l2 = s.Basis2().LookupOrNeg(q2);

   if (l1 >= 0 && l2 >= 0)
      s(q1, q2) = x;
}

inline
LatticeSite BosonU1(int MaxN, std::string const& Sym1 = "N")
{
   SymmetryList Symmetry(Sym1+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   // Setup the site basis
   for (int n = 0; n <= MaxN; ++n)
   {
      std::string q = boost::lexical_cast<std::string>(n);
      Basis.push_back(q, QN(n));
   }

   LatticeSite Site;

   // identity
   SiteOperator& I = Site["I"];
   I = SiteOperator::Identity(Basis);

   // creation
   SiteOperator& BH = Site["BH"];
   BH = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
     SetMatElementU1(BH, n+1, n, std::sqrt(double(n + 1)));
   }

   // annihilation
   SiteOperator& B = Site["B"];
   B = adjoint(BH);

   // parity, is identity
   SiteOperator& P = Site["P"];
   P = I;

   // number operator
   SiteOperator& N = Site["N"];
   N = prod(BH, B, QN(0));

   // N2 is N*(N-1)
   SiteOperator& N2 = Site["N2"];
   N2 = prod(N, N-I, QN(0));

   // Reflection operator is identity
   SiteOperator& R = Site["R"];
   R = I;

   // U is (-1)^N
   SiteOperator& U = Site["U"];
   U = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
      SetMatElementU1(U, n, n, minus1pow(n));
   }

   // projections onto each state
   for (int n = 0; n <= MaxN; ++n)
   {
      SiteOperator X(Basis, QN(0), LatticeCommute::Bosonic);
      SetMatElementU1(X, n, n, 1);
      std::string OpName = std::string("P_")
	+ boost::lexical_cast<std::string>(n);
      Site[OpName] = X;
   }

   return Site;
}