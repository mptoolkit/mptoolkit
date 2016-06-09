// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-canonical.cpp
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

#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpstate.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp/copyright.h"
#include "common/environment.h"

using QuantumNumbers::QuantumNumberList;
using LinearAlgebra::Matrix;

std::pair<MPStateComponent, MPStateComponent>
ConstructMaxEntangledState(BasisList const& b1, BasisList const& b2)
{
   CHECK_EQUAL(b1.size(), b2.size())("The real/aux pair of local sites has a different size basis");

   // Start with a right basis that is the identity quantum number
   VectorBasis rb(b1.GetSymmetryList());
   rb.push_back(QuantumNumbers::QuantumNumber(rb.GetSymmetryList()), 1);

   // The intermediate basis is the same as b2
   VectorBasis ib(b2);

   // The left basis is the tensor of b1[i] * b2[i]
   VectorBasis lb(b1.GetSymmetryList());
   for (unsigned i = 0; i < b1.size(); ++i)
   {
      QuantumNumberList ql = transform_targets(b2[i], b1[i]);
      for (unsigned q = 0; q < ql.size(); ++q)
         lb.push_back(ql[q], 1);
   }

   // Now make another pass and construct the MPS
   MPStateComponent rA(b2, ib, rb);
   for (unsigned i = 0; i < b2.size(); ++i)
   {
      rA[i](i, 0) = LinearAlgebra::Matrix<std::complex<double> >(1,1,1.0);
   }

   MPStateComponent lA(b1, lb, ib);
   unsigned j = 0;
   for (unsigned i = 0; i < b1.size(); ++i)
   {
      QuantumNumberList ql = transform_targets(b2[i], b1[i]);
      for (unsigned q = 0; q < ql.size(); ++q)
      {
         // TODO: The coefficient for the SU(2) case is just a guess
         lA[i](j,i) = LinearAlgebra::Matrix<std::complex<double> >
            (1,1, double(degree(ql[q])) / (degree(b2[i]) * degree(b1[i])));
         ++j;
      }
   }

   // Finally, it is possible that the same quantum number is repeated more than once
   // in the left basis.  Collapse it.
   MatrixOperator U = CollapseBasis(lA.Basis1());
   lA = prod(U, lA);

   return std::make_pair(lA, rA);
}

LinearWavefunction ConstructCanonicalWavefunction(Lattice const& L)
{
   std::vector<BasisList> AllBasis(L.size());
   for (int i = 0; i < L.size(); ++i)
   {
      AllBasis[i] = L[i+1].Basis1().Basis();
   }

   LinearWavefunction Result(L.GetSymmetryList());

   VectorBasis CurrentBasis(L.GetSymmetryList());
   CurrentBasis.push_back(QuantumNumbers::QuantumNumber(L.GetSymmetryList()), 1);
   MatrixOperator U = MatrixOperator::make_identity(CurrentBasis);
   Lattice::const_iterator I = L.end();
   for (int i = AllBasis.size()-1; i >= 0; i -=2)
   {
      std::cout << "Working (" << ((i+1)/2) << ") ..." << std::endl;
      DEBUG_TRACE("here")(CurrentBasis);
      std::pair<MPStateComponent, MPStateComponent> APair = 
         ConstructMaxEntangledState(AllBasis[i-1], AllBasis[i]);
      // make the tensor product with the identity operator in the CurrentBasis
      MatrixOperator Ident = MatrixOperator::make_identity(CurrentBasis);
      ProductBasis<VectorBasis, VectorBasis> BR = make_product_basis(CurrentBasis, APair.second.Basis2());
      ProductBasis<VectorBasis, VectorBasis> BI = make_product_basis(CurrentBasis, APair.first.Basis2());
      ProductBasis<VectorBasis, VectorBasis> BL = make_product_basis(CurrentBasis, APair.first.Basis1());

      U = CollapseBasis(BL.Basis());
      MPStateComponent NewA2(APair.second.SiteBasis(), BI.Basis(), BR.Basis());
      MPStateComponent NewA1(APair.first.SiteBasis(), U.Basis1(), BI.Basis());
      for (unsigned j = 0; j < APair.first.SiteBasis().size(); ++j)
      {
         NewA2[j] = tensor_prod(Ident, APair.second[j], BI, BR);
         NewA1[j] = U * tensor_prod(Ident, APair.first[j], BL, BI);
      }
      // Now add the A matrices to the wavefunction
      Result.push_front(NewA2);
      Result.push_front(NewA1);
      // update CurrentBasis
      CurrentBasis = NewA1.Basis1();
   }

   return Result;
}

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-canonical <lattice> <quantum number> <outfile>\n";
      return 1;
   }

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(argv[3], 1, PageSize, CacheSize);
   pvalue_ptr<OperatorList> OpList = pheap::ImportHeap(argv[1]);

   QuantumNumber Q(OpList->GetSymmetryList(), std::string(argv[2]));

   MPWavefunction Psi = ConstructCanonicalWavefunction(OpList->GetLattice());
   std::cout << "Projecting onto quantum number " << Q << std::endl;
   project(Psi, Q);

   std::cout.precision(20);
   double Normalization =  norm_2_sq(Psi);
   std::cout << "\nNormalization factor is " << Normalization << '\n';
   Psi.normalize();

   Psi.Attributes()["Beta"] = 0;
   Psi.Attributes()["Normalization"] = Normalization;

   pvalue_ptr<MPWavefunction> Ret = new MPWavefunction(Psi);
   pheap::ShutdownPersistent(Ret);
}
