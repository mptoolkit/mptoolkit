// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// junk/kondo-exactdiag.cpp
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


#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1u1.h"
#include "models/spin-u1.h"
#include "tensor/regularize.h"
#include "tensor/tensor_exponential.h"

typedef IrredTensor<std::complex<double> > OpType;

OpType CHup, CHdown, Cup, Cdown, CHupP, CHdownP, CupP, CdownP, BandI;
OpType BandSz, BandSp, BandSm, SpinSp, SpinSm, SpinSz, SpinI;

#define SIGN -

OpType MakeHopping(OpType const& LeftIdent)
{
   return tensor_prod(tensor_prod(LeftIdent, CHupP), Cup)
      + tensor_prod(tensor_prod(LeftIdent, CHdownP), Cdown)
      SIGN tensor_prod(tensor_prod(LeftIdent, CupP), CHup)
      SIGN tensor_prod(tensor_prod(LeftIdent, CdownP), CHdown);
}

OpType MakeHopping(OpType const& LeftIdent, OpType const& SpinIdent)
{
   return 
      tensor_prod(tensor_prod(tensor_prod(tensor_prod(LeftIdent, CHupP), SpinIdent), SpinIdent), Cup)
      + tensor_prod(tensor_prod(tensor_prod(tensor_prod(LeftIdent, CHdownP), SpinIdent), SpinIdent), Cdown)
      SIGN tensor_prod(tensor_prod(tensor_prod(tensor_prod(LeftIdent, CupP), SpinIdent), SpinIdent), CHup)
      SIGN tensor_prod(tensor_prod(tensor_prod(tensor_prod(LeftIdent, CdownP), SpinIdent), SpinIdent), CHdown);
}

OpType ApplyIdent(OpType const& Op, OpType const& Ident, int n = 1)
{
   OpType Result = Op;
   for (int i = 0; i < n; ++i)
   {
      Result = tensor_prod(Result, Ident);
   }
   return Result;
}

MatrixOperator MakeMatrixOp(OpType const& Op)
{
   VectorBasis B1(Op.Basis1());
   VectorBasis B2(Op.Basis2());
   MatrixOperator Result(B1, B2, Op.TransformsAs());
   Result.data() = scalar(LinearAlgebra::Matrix<double>(1,1,1.0)) * Op.data();
   return Result;
}

int main()
{
   SiteBlock Band = CreateU1HubbardSite();
   SiteBlock Spin = CreateU1SpinSite(0.5);

   QuantumNumbers::SymmetryList MySL("N:U(1),Sz:U(1)");
   Band.CoerceSymmetryList(MySL);
   Spin.CoerceSymmetryList(MySL);

   QuantumNumber Ident(MySL, "0,0");

   // operators in the fermion band
   CHup = Band["CHup"];
   CHupP = prod(OpType(Band["CHup"]), OpType(Band["P"]), QuantumNumber(MySL, "1,0.5"));
   CHdown = Band["CHdown"];
   CHdownP = prod(OpType(Band["CHdown"]), OpType(Band["P"]), QuantumNumber(MySL, "1,-0.5"));

   Cup = Band["Cup"];
   CupP = prod(OpType(Band["Cup"]), OpType(Band["P"]), QuantumNumber(MySL, "-1,-0.5"));
   Cdown = Band["Cdown"];
   CdownP = prod(OpType(Band["Cdown"]), OpType(Band["P"]), QuantumNumber(MySL, "-1,0.5"));

   BandSp = Band["Sp"];
   BandSm = Band["Sm"];
   BandSz = Band["Sz"];

   BandI = Band["I"];

   // operators in the spin site
   SpinSp = Spin["Sp"];
   SpinSm = Spin["Sm"];
   SpinSz = Spin["Sz"];

   SpinI = Spin["I"];

   // assemble the Hamiltonian
   OpType Hop1 = tensor_prod(CHupP, Cup, Ident)
      + tensor_prod(CHdownP, Cdown, Ident)
      SIGN tensor_prod(CupP, CHup, Ident)
      SIGN tensor_prod(CdownP, CHdown, Ident);

   OpType LeftIdent = ApplyIdent(BandI, BandI);
   OpType MidIdent = ApplyIdent(LeftIdent, SpinI, 2);

   // incorporate the 2 spin sites
   OpType Hopping = ApplyIdent(ApplyIdent(Hop1, SpinI, 2), BandI)
      + MakeHopping(BandI, SpinI);

   // and the remaining site
   Hopping = ApplyIdent(Hopping, BandI)
      + MakeHopping(MidIdent);

   // now the spin interactions
   OpType SpinJ1 = 0.5 * (tensor_prod(tensor_prod(BandI, BandSp), SpinSm) 
			  + tensor_prod(tensor_prod(BandI, BandSm), SpinSp))
      + tensor_prod(tensor_prod(BandI, BandSz), SpinSz);
   SpinJ1 = ApplyIdent(SpinJ1, SpinI);
   SpinJ1 = ApplyIdent(SpinJ1, BandI, 2);

   OpType SpinJ2 = ApplyIdent(LeftIdent, SpinI);
   SpinJ2 = 0.5 * (tensor_prod(tensor_prod(SpinJ2, SpinSp), BandSm) 
		   + tensor_prod(tensor_prod(SpinJ2, SpinSm), BandSp))
      + tensor_prod(tensor_prod(SpinJ2, SpinSz), BandSz);
   SpinJ2 = ApplyIdent(SpinJ2, BandI);
   
   OpType SzSz = tensor_prod(tensor_prod(LeftIdent, SpinSz), SpinSz);
   SzSz = ApplyIdent(SzSz, BandI, 2);

   OpType SS = 0.5 * (tensor_prod(tensor_prod(LeftIdent, SpinSp), SpinSm) 
		      + tensor_prod(tensor_prod(LeftIdent, SpinSm), SpinSp))
      + tensor_prod(tensor_prod(LeftIdent, SpinSz), SpinSz);
   SS = ApplyIdent(SS, BandI, 2);

   OpType Hamiltonian = -0.5*Hopping + 0.8*(SpinJ1+SpinJ2);

   TRACE(norm_frob_sq(Hamiltonian-adjoint(Hamiltonian)));

   // add a constant to make the trace non-zero
   //Hamiltonian += ApplyIdent(MidIdent, BandI, 2);

//   MatrixOperator H = MakeMatrixOp(Hamiltonian);

//   MatrixOperator Reg = Regularize(H.Basis1());
//   MatrixOperator HReg = triple_prod(Reg, H, herm(Reg));

#if 0
   for (double Beta = 1; Beta <= 1000; Beta += 1)
   {
//      MatrixOperator m = HReg;
//      m *= Beta;
      OpType m = Exponentiate(-Beta*Hamiltonian);
      double z = trace(m).real();
      double mag = trace(prod(m,SS,Ident)).real();
      double magz = trace(prod(m,SzSz,Ident)).real();
      TRACE(Beta)(z)(mag/z)(magz/z);
   }
#endif

#if 0
   std::cout << '{';
   bool first = true;
   for (const_iterator<OpType>::type I = iterate(Hamiltonian); I; ++I)
   {
      for (const_inner_iterator<OpType>::type J = iterate(I); J; ++J)
      {
	 if (first)
	    first = false;
	 else
	    std::cout << ",\n";

	 std::cout << "Rule[{" << (J.index1()+1) << ',' << (J.index2()+1) << "}, " << J->real() << "]";
      }
   }
   std::cout << "}\n";
#endif

#if 1
   double Beta = 1;
   OpType Ex = Exponentiate(-Beta*Hamiltonian);
   Ex *= 1.0 / norm_frob(Ex);

   std::cout.precision(12);
   OpType m = Ex;
   for (int i = 1; i <= 500; ++i)
   {
      double z = trace(m).real();
      double mag = trace(prod(m,SS,Ident)).real();
      std::cout << (Beta*i) << "  " << z << "  " << (mag/z) << '\n';
      m = prod(Ex, m, Ident);
      m *= 1.0 / norm_frob(m);
   }
#endif
}
