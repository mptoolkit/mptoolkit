// -*- C++ -*- $Id$

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/z2.h"
#include "siteoperator/block.h"

typedef Block<SiteOperator> SiteBlock;

std::string StateName(half_int s, bool Symmetric)
{
   std::string Result = boost::lexical_cast<std::string>(s);
   Result += Symmetric ? 's' : 'a';
   return Result;
}

inline
SiteBlock CreateZ2SpinSite(half_int Spin)
{
   SymmetryList Symmetry("Z:Z2");
   QuantumNumbers::QNConstructor<QuantumNumbers::Z2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator Sz, Sx, Sy, mSz, R, P, I, Z;
   SiteBlock Block;

   PRECONDITION(Spin >= 0)("Spin must be >= 0")(Spin);

   half_int Base = 0.5; // lowest possible quantum number
   if (Spin.is_integral())
   {
      Base = 1;
      Basis.push_back(StateName(0,true), QN(1));
   }
   for (half_int s = Base; s <= Spin; ++s)
   {
      Basis.push_back(StateName(s, true), QN(1));
      Basis.push_back(StateName(s, false), QN(-1));
   }

   Sx = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   Sy = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);
   P = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   Z = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   I = SiteOperator::Identity(Basis);

   P = I; // bosonic, so the parity operator is identity
   R = I; // spatial reflection


   // Set the matrix elements for Sz
   // We only set half the matrix elements of Sx,Sy, and at the end
   // make them hermitian
   // Since Base > 0 we do the s=0 components later (these are special anyway)
   for (half_int s = Base; s <= Spin; ++s)
   {
      Z(StateName(s, true), StateName(s, true)) = 1;

      Z(StateName(s, false), StateName(s, false)) = -1;
      
      Sz(StateName(s, true), StateName(s, false)) = s.to_double();
      Sz(StateName(s, false), StateName(s, true)) = s.to_double();
           
      if (s < Spin)
      {
         Sx(StateName(s, false), StateName(s+1, false)) 
            = 0.5 * sqrt((Spin-s)*(Spin+s+1));
         Sx(StateName(s, true), StateName(s+1, true)) 
            = 0.5 * sqrt((Spin-s)*(Spin+s+1));

         Sy(StateName(s, false), StateName(s+1, true)) 
            = std::complex<double>(0.0, -0.5 * sqrt((Spin-s)*(Spin+s+1)));
         Sy(StateName(s, true), StateName(s+1, false)) 
            = std::complex<double>(0.0, -0.5 * sqrt((Spin-s)*(Spin+s+1)));
      }
   }

   // do the s == 0 part
   if (Spin.is_integral())
   {
      Z("0s", "0s") = 1;
      if (Spin > 0)
      {
         Sx("0s", "1s") 
            = sqrt(Spin*(Spin+1)/2);

         Sy("0s", "1a") 
            = std::complex<double>(0.0, -sqrt(Spin*(Spin+1)/2));
      }
   }

   Sx = Sx + adjoint(Sx);
   Sy = Sy + adjoint(Sy);

   if (!Spin.is_integral())
   {
      // Some of the matrix elements of Sx,Sy for s=0.5 are diagonal in m
      Sx(StateName(0.5, true), StateName(0.5, true)) 
         = 0.5 * (Spin+0.5);
      Sx(StateName(0.5, false), StateName(0.5, false)) 
         = -0.5 * (Spin+0.5);

      Sy(StateName(0.5, true), StateName(0.5, false)) 
         = std::complex<double>(0.0, -0.5 * (Spin+0.5));
      Sy(StateName(0.5, false), StateName(0.5, true)) 
         = std::complex<double>(0.0, 0.5 * (Spin+0.5));
   }

#if 0
   Sx("s", "s") =  0.5;
   Sx("a", "a") = -0.5;

   Sz("s", "a") =  0.5;
   Sz("a", "s") =  0.5;

   Sy("s", "a") = std::complex<double>(0.0, -0.5);
   Sy("a", "s") = std::complex<double>(0.0,  0.5);

   Z("s", "s") =  1.0;
   Z("a", "a") = -1.0;
#endif

   Block["I"] = I;
   Block["P"] = P;
   Block["R"] = R;
   Block["Sx"] = Sx;
   Block["Sy"] = Sy;
   Block["Sz"] = Sz;
   Block["Sz2"] = Sz * Sz;
   Block["Z"] = Z;
   return Block;
}
