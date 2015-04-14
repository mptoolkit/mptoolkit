// -*- C++ -*- $Id$

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"




inline
LatticeSite CreateU1SpinSite(half_int Spin, std::string const& Sym = "Sz")
{
   SymmetryList Symmetry(Sym+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator Sp, Sm, Sz, mSz, R, P, I, Spp, Smm, Sz2;
   LatticeSite Site;

   std::map<half_int, std::string> SpinBasis;
   for (half_int s = -Spin; s <= Spin; ++s)
   {
      SpinBasis[s] = boost::lexical_cast<std::string>(s);
      Basis.push_back(SpinBasis[s], QN(s));
   }

   Sp = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sz2 = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   mSz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Spp = SiteOperator(Basis, QN(to_int(2*Spin)), LatticeCommute::Bosonic);
   Smm = SiteOperator(Basis, QN(to_int(-2*Spin)), LatticeCommute::Bosonic);

   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   // for the mSz operator, we want (-1)^Sz, but for half-integer spin
   // this is imaginary, so we want to multiply it by i in that case
   half_int Fudge = Spin.is_integral() ? 0.0 : 0.5;
   for (half_int s = -Spin; s <= Spin; ++s)
   {
      I(SpinBasis[s], SpinBasis[s]) = 1.0;
      P(SpinBasis[s], SpinBasis[s]) = 1.0;
      R(SpinBasis[s], SpinBasis[s]) = 1.0;
      Sz(SpinBasis[s], SpinBasis[s]) = s.to_double();
      mSz(SpinBasis[s], SpinBasis[s]) = minus1pow((s+Fudge).to_int());
   }

   // Sp and Sm operators
   for (half_int s = -Spin; s < Spin; ++s)
   {
      Sp(SpinBasis[s+1], SpinBasis[s]) = std::sqrt((Spin - s) * (Spin + s + 1));
   }
   Sm = adjoint(Sp);

   // Spp and Smm operators - maximal spin flips |-Spin><Spin| and |Spin><-Spin|
   // with amplitude 1.0
   Spp(SpinBasis[Spin], SpinBasis[-Spin]) = 1.0;
   Smm(SpinBasis[-Spin], SpinBasis[Spin]) = 1.0;

   Sz2 = Sz * Sz;
   
   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sz"] = Sz;
   Site["mSz"] = mSz;
   Site["Sz2"] = prod(Sz, Sz, QN(0));
   Site["Spp"] = Spp;
   Site["Smm"] = Smm;
   Site["Sz2"] = Sz2;
   return Site;
}
