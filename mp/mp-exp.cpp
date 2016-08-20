// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-exp.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/linearoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "linearalgebra/matrix_utility.h"

struct ExpProd
{
   typedef std::complex<double> result_type;
   typedef SimpleOperator const& first_argument_type;
   typedef SimpleOperator const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.Basis1().size(), 1);
      DEBUG_PRECONDITION_EQUAL(y.Basis2().size(), 1);
      return trace(prod(x, y, QuantumNumbers::QuantumNumber(x.GetSymmetryList())));
   }
};

void TwoSiteScalarApply(MPOpComponent& A, MPOpComponent& B, std::complex<double> x)
{
   PRECONDITION_EQUAL(A.Basis1().size(), 1);
   PRECONDITION_EQUAL(B.Basis2().size(), 1);

   QuantumNumber Ident(A.GetSymmetryList());
   PRECONDITION_EQUAL(A.Basis1()[0], Ident);
   PRECONDITION_EQUAL(B.Basis2()[0], Ident);

   SimpleOperator C = ExpandBasis2(A);
   A = prod(A, C);

   C = ExpandBasis1(B);
   B = prod(C, B);

   ProductBasis<BasisList, BasisList> PB(A.SiteBasis(), B.SiteBasis());

   SimpleOperator TP;
   for (MPOpComponent::const_iterator I = A.begin(); I != A.end(); ++I)
   {
      TP += tensor_prod(I->second, B[I->first], PB, PB, Ident, ExpProd());
   }

   // now we apply the exponential to TP
   LinearAlgebra::Matrix<std::complex<double> > M(TP.data());

   LinearAlgebra::Vector<std::complex<double> > EVec = DiagonalizeHermitian(M);
   for (std::size_t i = 0; i < EVec.size(); ++i)
   {
      EVec[i] = exp(EVec[i] * x);
   }

   zero_all(TP.data());
   M = transpose(M) * diagonal_matrix(EVec) * M;
   TRACE(M);
   for (std::size_t i = 0; i < size1(M); ++i)
   {
      for (std::size_t j = 0; j < size2(M); ++j)
      {
         if (norm_2(M(i,j)) > 1E-10) set_element(TP.data(), i, j, M(i,j));
      }
   }

   ProductBasis<BasisList, BasisList> alpha(A.SiteBasis(), adjoint(A.SiteBasis()));
   ProductBasis<BasisList, BasisList> beta(B.SiteBasis(), adjoint(B.SiteBasis()));

   SimpleOperator Partial = partial_transpose(TP, PB, PB, alpha, beta);

   M = Partial.data();
   LinearAlgebra::Matrix<std::complex<double> > U, Vt;
   LinearAlgebra::Vector<double> D;
   LinearAlgebra::SingularValueDecomposition(M, U, D, Vt);

   BasisList AuxBasis(TP.GetSymmetryList());

   for (std::size_t i = 0; i < size(D) && D[i] > 1E-10; ++i)
   {

      int iMax = 0;
      double Max = LinearAlgebra::norm_2(U(iMax,i));

      for (std::size_t k = 1; k < size1(U); ++k)
      {
         if (LinearAlgebra::norm_2(U(k,i)) > Max)
         {
            iMax = k;
            Max = LinearAlgebra::norm_2(U(k,i));
         }
      }

      AuxBasis.push_back(alpha[iMax]);
   }

   A = MPOpComponent(A.SiteBasis(), A.Basis1(), AuxBasis);
   B = MPOpComponent(B.SiteBasis(), AuxBasis, B.Basis2());

   for (std::size_t i = 0; i < size(D); ++i)
   {
      if (D[i] < 1E-10) continue;

      for (std::size_t k = 0; k < size1(U); ++k)
      {
         if (LinearAlgebra::norm_2(U(k,i)) > 1E-10)
         {
            int l,m;
            std::tie(l,m) = alpha.rmap(k);

            A[alpha[k]].data()(l,m)(1,k) = U(k,i) * D[i];
         }
      }

      for (std::size_t k = 0; k < size2(Vt); ++k)
      {
         if (LinearAlgebra::norm_2(Vt(i,k)) > 1E-10)
         {
            int l,m;
            std::tie(l,m) = beta.rmap(k);

            B[beta[k]].data()(l,m)(k,1) = Vt(i,k) / std::sqrt(double(degree(B.SiteBasis()[l])));
         }
      }
   }
}

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-exp <lattice> <operator> <real> <imag> <output-operator>\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);
   std::string Operator = argv[2];

   SplitOperator Op = (*System)[Operator];
   std::string OutOp = argv[5];

   double Real = boost::lexical_cast<double>(argv[3]);
   double Imag = boost::lexical_cast<double>(argv[4]);

   while (Op.RightSize() > 1)
   {
      if (Op.Left().Basis1().size() == 1 && is_scalar(Op.Left().Basis1()[0])
          && Op.Right().Basis2().size() == 1 && is_scalar(Op.Right().Basis2()[0]))
      {
         Op.Left() = prod(Op.Left(), Op.Center());
         TwoSiteScalarApply(Op.Left(), Op.Right(), std::complex<double>(Real, Imag));
         Op.Center() = TruncateBasis2(Op.Left());
         Op.RotateRight();
      }
      Op.RotateRight();
   }

   Op.Right() = prod(Op.Center(), Op.Right());
   Op.Center() = TruncateBasis1(Op.Right());

   while (Op.LeftSize() > 1)
   {
      Op.RotateLeft();
   }

   (*System.mutate())[OutOp] = Op.AsLinearOperator();
   pheap::ShutdownPersistent(System);
}
