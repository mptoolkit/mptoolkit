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
LatticeSite CreateBoseHubbardSpinlessU1Site(int MaxN, std::string const& Sym1 = "N")
{
   SymmetryList Symmetry(Sym1+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator B, BH, P, R, N, N2, Q, I, U;
   LatticeSite Site;

   // Setup the site basis
   for (int n = 0; n <= MaxN; ++n)
   {
      std::string q = boost::lexical_cast<std::string>(n);
      Basis.push_back(q, QN(n));
   }

   BH = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
     SetMatElementU1(BH, n+1, n, std::sqrt(double(n + 1)));
   }

   Site["BH"] = BH;

   B = adjoint(BH);
   Site["B"] = B;
   I = SiteOperator::Identity(Basis);
   P = I;
   Site["I"] = I;
   N = prod(BH, B, QN(0));
   Site["N"] = N;
   N2 = prod(N, N-I, QN(0));
   Site["N2"] = N2;
   R = I;
   Site["R"] = R;
   Site["P"] = P;

   U = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
      SetMatElementU1(U, n, n, minus1pow(n));
   }
   Site["U"] = U;

   Site["BN"] = B*N;
   Site["BN2"] = B*N*N;
   Site["BN3"] = B*N*N*N;
   Site["BU"] = B*U;

   Site["BHN"] = BH*N;
   Site["BHN2"] = BH*N*N;
   Site["BHN3"] = BH*N*N*N;
   Site["BHU"] = BH*U;

   // projections onto each state
   for (int n = 0; n <= MaxN; ++n)
   {
      SiteOperator X(Basis, QN(0), LatticeCommute::Bosonic);
      SetMatElementU1(X, n, n, 1);
      std::string OpName = std::string("P_")
	+ boost::lexical_cast<std::string>(n);
      Site[OpName] = X;
   }

   DEBUG_TRACE(BH)(B)(I)(N)(N2)(U);

   return Site;
}
