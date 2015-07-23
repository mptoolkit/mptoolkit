// -*- C++ -*- $Id$

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"

inline
LatticeSite SpinU1(half_int Spin, std::string const& Sym = "Sz")
{
   SymmetryList Symmetry(Sym+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator Sp, Sm, Sz, R, P, I;
   LatticeSite Site("U(1) Spin "+to_string_fraction(Spin));

   std::map<half_int, std::string> SpinBasis;
   for (half_int s = -Spin; s <= Spin; ++s)
   {
      SpinBasis[s] = boost::lexical_cast<std::string>(s);
      Basis.push_back(SpinBasis[s], QN(s));
   }

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("Sp"  , "raising operator")
      ("Sm"  , "lowering operator")
      ("Sz"  , "z-component of spin")
      ;

   Sp = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   for (half_int s = -Spin; s <= Spin; ++s)
   {
      I(SpinBasis[s], SpinBasis[s]) = 1.0;
      P(SpinBasis[s], SpinBasis[s]) = 1.0;
      R(SpinBasis[s], SpinBasis[s]) = 1.0;
      Sz(SpinBasis[s], SpinBasis[s]) = s.to_double();
   }

   // Sp and Sm operators
   for (half_int s = -Spin; s < Spin; ++s)
   {
      Sp(SpinBasis[s+1], SpinBasis[s]) = std::sqrt((Spin - s) * (Spin + s + 1));
   }
   Sm = adjoint(Sp);

   // Example of defining a named constant (argument)
   Site.arg("Spin") = Spin.to_double();

   // Example of defining a function.  The first parameter has a default value
   Site.func("Uz")(arg("theta") = math_const::pi) = "exp(theta*i*Sz)";

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sz"] = Sz;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
