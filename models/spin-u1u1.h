// -*- C++ -*- $Id: spin-u1.h 1011 2009-04-22 14:36:31Z ianmcc $
// SU(3) spin chain in the fundamental representation,
// broken down to U(1)xU(1)

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/u1.h"
#include "siteoperator/block.h"

typedef Block<SiteOperator> SiteBlock;

inline
SiteBlock CreateU1U1SpinSite()
{
   SymmetryList Symmetry("Sz:U(1),Qz:U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::U1> 
      QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator L3, L8;
   SiteOperator Tp, Tm, Up, Um, Vp, Vm, Sz, Qz;

   SiteOperator R, P, I, U;
   SiteBlock Block;

   // Basis is labelled by Sz = L3/2 and Qz = L8 * sqrt(3)/2
   // we label the local basis as 1,-1,0, the eigenvalue of L3
   Basis.push_back("1", QN(0.5, 0.5));
   Basis.push_back("-1", QN(-0.5, 0.5));
   Basis.push_back("0", QN(0, -1));

   L3 = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   L8 = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   Sz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Qz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // T+ = L1 + i*L2
   // T- = L1 - i*L2
   Tp = SiteOperator(Basis, QN(1,0), LatticeCommute::Bosonic);
   Tm = SiteOperator(Basis, QN(-1,0), LatticeCommute::Bosonic);

   // V+ = L4 + i*L5
   // V- = L4 - i*L5
   Vp = SiteOperator(Basis, QN(0.5,1.5), LatticeCommute::Bosonic);
   Vm = SiteOperator(Basis, QN(-0.5,-1.5), LatticeCommute::Bosonic);

   // U+ = L6 + i*L7
   // U- = L6 - i*L7
   Up = SiteOperator(Basis, QN(-0.5,1.5), LatticeCommute::Bosonic);
   Um = SiteOperator(Basis, QN(0.5,-1.5), LatticeCommute::Bosonic);

   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   U = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   L3("1", "1")   =  1.0;
   L3("-1", "-1") = -1.0;

   L8("1", "1")   =  1.0 / std::sqrt(3.0);
   L8("-1", "-1") =  1.0 / std::sqrt(3.0);
   L8("0", "0")   = -2.0 / std::sqrt(3.0);

   Sz("1", "1")   =  0.5;
   Sz("-1", "-1") = -0.5;
   
   Qz("1", "1")   =  0.5;
   Qz("-1", "-1") =  0.5;
   Qz("0", "0")   = -1.0;

   Tp("1", "-1")  =  2.0;
   Tm("-1", "1")  =  2.0;

   Vp("1", "0")   =  2.0;
   Vm("0", "1")   =  2.0;

   Up("-1", "0")  =  2.0;
   Um("0", "-1")  =  2.0;

   I("1", "1")    =  1.0;
   I("-1", "-1")  =  1.0;
   I("0", "0")    =  1.0;

   P = I;
   R = I;

   U("1", "1")    = -1.0;
   U("-1", "-1")  = -1.0;
   U("0", "0")    =  1.0;
   
   Block["I"] = I;
   Block["P"] = P;
   Block["R"] = R;
   Block["U"] = U;
   Block["L3"] = L3;
   Block["L8"] = L8;
   Block["Sz"] = Sz;
   Block["Qz"] = Qz;
   Block["Tp"] = Tp;
   Block["Tm"] = Tm;
   Block["Vp"] = Vp;
   Block["Vm"] = Vm;
   Block["Up"] = Up;
   Block["Um"] = Um;
   return Block;
}
