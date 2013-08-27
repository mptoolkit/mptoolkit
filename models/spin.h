// -*- C++ -*- $Id$

#include "siteoperator/latticesite.h"
#include "quantumnumbers/u1.h"

inline
LatticeSite CreateSpinSite(half_int Spin)
{
   SymmetryList Symmetry("S:Null");
   SiteBasis Basis(Symmetry);
   SiteOperator Sp, Sm, Sz, Sx, Sy, mSz, mSx, mSy, R, P, I;
   LatticeSite Site;

   QuantumNumbers::QuantumNumber q(Symmetry); // no symmetries, only one quantum number

   std::map<half_int, std::string> SpinBasis;
   for (half_int s = -Spin; s <= Spin; ++s)
   {
      SpinBasis[s] = boost::lexical_cast<std::string>(s);
      Basis.push_back(SpinBasis[s], q);
   }

   Sp = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sx = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sy = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   mSz = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   mSy = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   mSz = SiteOperator(Basis, q, LatticeCommute::Bosonic);

   P = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   R = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   I = SiteOperator(Basis, q, LatticeCommute::Bosonic);

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

   for (half_int s = -Spin; s < Spin; ++s)
   {
      Sp(SpinBasis[s+1], SpinBasis[s]) = std::sqrt((Spin - s) * (Spin + s + 1));
   }

   Sm = adjoint(Sp);

   Sx = 0.5 * (Sp + Sm);
   Sy = std::complex<double>(0.0, -0.5) * (Sp - Sm);

   // the mSx and mSz operators.  I don't know the general formula, so
   // we'll just cover the cases I need
   if (Spin == 1)
   {
      mSx = I - 2.0*Sx*Sx;
      mSy = I - 2.0*Sy*Sy;
   }
   else if (Spin == 2)
   {
      mSx = I - (8.0/3.0)*Sx*Sx + (2.0/3.0)*Sx*Sx*Sx*Sx;
      mSy = I - (8.0/3.0)*Sy*Sy + (2.0/3.0)*Sy*Sy*Sy*Sy;
   }
   
   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sx"] = Sx;
   Site["Sy"] = Sy;
   Site["Sz"] = Sz;
   Site["mSz"] = mSz;
   Site["mSy"] = mSy;
   Site["mSx"] = mSx;
   Site["Sz2"] = prod(Sz, Sz, q);
   return Site;
}
