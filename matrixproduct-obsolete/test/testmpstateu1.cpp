// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/test/testmpstateu1.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/mpstate.h"
#include "quantumnumbers/u1.h"

int main()
{
   // local site basis
   SymmetryList Symmetry("Sz:U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   BasisList B(Symmetry);
   B.push_back(QN(-1));
   B.push_back(QN(0));
   B.push_back(QN(1));

   // auxiliary basis
   VectorBasis VB(Symmetry);
   VB.push_back(QN(-0.5), 1);
   VB.push_back(QN(0.5), 1);

   MatrixOperator AuxIdent(VB, VB, QN(0));
   set_element(AuxIdent.data(), 0, 0, LinearAlgebra::identity_matrix<double>(1));
   set_element(AuxIdent.data(), 1, 1, LinearAlgebra::identity_matrix<double>(1));

   MPStateComponent x(B, VB, VB);

   MPStateComponent FB1 = MPStateComponent::ConstructFullBasis1(B, VB);

   MatrixOperator FullIdent1(FB1.Basis1(), FB1.Basis1(), QN(0));
   for (std::size_t i = 0; i < FB1.Basis1().size(); ++i)
   {
      set_element(FullIdent1.data(), i, i, 
		  LinearAlgebra::identity_matrix<double>(1));
   }

   CHECK_CLOSE(scalar_prod(FB1, herm(FB1)), FullIdent1);

   MPStateComponent FB2 = MPStateComponent::ConstructFullBasis2(VB, B);

   // FB1.Basis1() and FB2.Basis2() are different here because
   // they are constructed as a tensor product x \otimes y
   // and y \otimes x.  So we cannot compare them and must construct
   // another identity matrix.

   MatrixOperator FullIdent2(FB2.Basis2(), FB2.Basis2(), QN(0));
   for (std::size_t i = 0; i < FB2.Basis2().size(); ++i)
   {
      set_element(FullIdent2.data(), i, i, 
		  LinearAlgebra::identity_matrix<double>(1));
   }

   CHECK_CLOSE(scalar_prod(herm(FB2), FB2), FullIdent2);


   SimpleOperator SiteIdent(B, B, QN(0));
   SiteIdent(0,0) = 1;
   SiteIdent(1,1) = 1;
   SiteIdent(2,2) = 1;

   CHECK_CLOSE(operator_prod(SiteIdent, FB1, herm(FB1)), FullIdent1);
   CHECK_CLOSE(operator_prod(SiteIdent, herm(FB2), FB2), FullIdent2);

   CHECK_CLOSE(operator_prod(SiteIdent, FB1, AuxIdent, herm(FB1)), FullIdent1);
   CHECK_CLOSE(operator_prod(SiteIdent, herm(FB2), AuxIdent, FB2), FullIdent2);
}
