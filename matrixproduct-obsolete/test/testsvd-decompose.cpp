// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/test/testsvd-decompose.cpp
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

#include "matrixproduct/mpstate.h"
#include "quantumnumbers/su2.h"
#include "quantumnumbers/u1.h"
#include "matrixproduct/density.h"

using namespace QuantumNumbers;
using namespace Tensor;

int main()
{
   // local site basis
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   BasisList Aux(Symmetry);
   Aux.push_back(QN(1));

   // auxiliary basis
   VectorBasis VB(Symmetry);
   VB.push_back(QN(0.5), 1);

   MatrixOperator AuxIdent(VB, VB, QN(0));
   set_element(AuxIdent.data(), 0, 0, LinearAlgebra::identity_matrix<double>(1));

   MPStateComponent x(Aux, VB, VB);

   MPStateComponent FB1 = MPStateComponent::ConstructFullBasis1(Aux, VB);

   MatrixOperator FullIdent1(FB1.Basis1(), FB1.Basis1(), QN(0));
   for (std::size_t i = 0; i < FB1.Basis1().size(); ++i)
   {
      set_element(FullIdent1.data(), i, i, 
		  LinearAlgebra::identity_matrix<double>(1));
   }

   CHECK_CLOSE(scalar_prod(FB1, herm(FB1)), FullIdent1);

   MPStateComponent FB2 =  MPStateComponent::ConstructFullBasis2(VB, Aux);

   MatrixOperator FullIdent2(FB2.Basis2(), FB2.Basis2(), QN(0));
   for (std::size_t i = 0; i < FB2.Basis2().size(); ++i)
   {
      set_element(FullIdent2.data(), i, i, 
		  LinearAlgebra::identity_matrix<double>(1));
   }

   CHECK_CLOSE(scalar_prod(herm(FB2), FB2), FullIdent2);

   Tensor::ProductBasis<BasisList, BasisList> PB(FB2.LocalBasis(), FB1.LocalBasis());

   MPStateComponent Combined(PB.Basis(), FB2.Basis1(), FB1.Basis2());
   for (unsigned i = 0; i < FB2.LocalBasis().size(); ++i)
   {
      for (unsigned j = 0; j < FB1.LocalBasis().size(); ++j)
      {
	 ProductBasis<BasisList, BasisList>::const_iterator IEnd = PB.end(i,j);
	 ProductBasis<BasisList, BasisList>::const_iterator I = PB.begin(i,j);
	 for ( ; I != IEnd; ++I)
	 {
	    Combined[*I] = prod(FB2[i], FB1[j], PB.Basis()[*I]);
	 }
      }
   }
   
   SingularDecomposition<MPStateComponent, MPStateComponent> SL(Combined, PB);
   for (SingularDecomposition<MPStateComponent, MPStateComponent>::const_iterator I = 
	   SL.begin(); I != SL.end(); ++I)
   {
      TRACE(I->Eigenvalue)(SL.Lookup(I->Subspace));
   }


   MPStateComponent A,B;
   MatrixOperator C;
   SL.ConstructMatrices(SL.begin(), SL.end(), A, C, B);

   TRACE(A)(C)(B);

   TRACE(scalar_prod(herm(A), A))(scalar_prod(B, herm(B)));

   MPStateComponent CC = local_tensor_prod(prod(A, C), B);
   TRACE(CC)(Combined);
   TRACE(norm_frob(CC-Combined));
   TRACE(FB2)(FB1);
}
