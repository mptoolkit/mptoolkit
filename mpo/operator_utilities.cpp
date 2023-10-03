// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/operator_utilities.cpp
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

#include "operator_utilities.h"

int
linear_dimension(SimpleOperator const& c)
{
   int n = 0;
   for (unsigned i = 0; i < c.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < c.Basis2().size(); ++j)
      {
         if (is_transform_target(c.Basis2()[j], c.TransformsAs(), c.Basis1()[i]))
         {
            ++n;
         }
      }
   }
   return n;
}

LinearAlgebra::Vector<std::complex<double>>
linearize(SimpleOperator const& c)
{
   int TotalSize = linear_dimension(c);
   LinearAlgebra::Vector<std::complex<double>> Result(TotalSize);
   int n = 0;
   for (unsigned i = 0; i < c.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < c.Basis2().size(); ++j)
      {
         if (is_transform_target(c.Basis2()[j], c.TransformsAs(), c.Basis1()[i]))
         {
            Result[n++] = c(i,j);
         }
      }
   }
   CHECK_EQUAL(n, TotalSize);
   return Result;
}

LinearAlgebra::Vector<std::complex<double>>
linearize(SimpleRedOperator const& c, QuantumNumberList const& QList)
{
   int TotalSize = linear_dimension(c, QList);
   LinearAlgebra::Vector<std::complex<double>> Result(TotalSize);
   int i = 0;
   for (auto q : QList)
   {
      SimpleOperator x = project(c, q);
      int n = linear_dimension(x);
      Result[LinearAlgebra::range(i, i+n)] = linearize(x);
      i += n;
   }
   return Result;
}

int
linear_dimension(SimpleRedOperator const& c, QuantumNumberList const& QList)
{
   int n = 0;
   for (auto q : QList)
   {
      n += linear_dimension(project(c,q));
   }
   return n;
}

SimpleOperator CollapseBasis(BasisList const& b)
{
   std::set<QuantumNumber> QN(b.begin(), b.end());
   BasisList NewB(b.GetSymmetryList(), QN.begin(), QN.end());
   SimpleOperator C(NewB, b);
   for (unsigned j = 0; j < b.size(); ++j)
   {
      unsigned i = std::find(NewB.begin(), NewB.end(), b[j]) - NewB.begin();
      C(i,j) = 1.0;
   }
   return C;
}

SimpleOperator ProjectBasis(BasisList const& b, QuantumNumbers::QuantumNumber const& q)
{
   BasisList NewB(q);
   SimpleOperator Result(NewB, b);
   for (unsigned j = 0; j < b.size(); ++j)
   {
      if (b[j] == q)
         Result(0, j) = 1.0;
   }
   return Result;
}

std::complex<double> PropIdent(SimpleOperator const& X, double UnityEpsilon)
{
   if (X.Basis1() != X.Basis2())
      return 0.0;
   SimpleOperator Ident = SimpleOperator::make_identity(X.Basis1());
   std::complex<double> x = inner_prod(Ident, X) / double(X.Basis1().total_degree());
   if (norm_frob_sq(X-x*Ident) > UnityEpsilon*UnityEpsilon)
      x = 0.0;
   return x;
}

SimpleOperator
swap_gate(BasisList const& B1, BasisList const& B2,
          ProductBasis<BasisList, BasisList> const& Basis_21,
          ProductBasis<BasisList, BasisList> const& Basis_12)
{
   QuantumNumbers::QuantumNumber Ident(B1.GetSymmetryList());

   SimpleOperator Result(Basis_21.Basis(), Basis_12.Basis(), Ident);
   for (unsigned i = 0; i < Basis_12.size(); ++i)
   {
      std::pair<int,int> x = Basis_12.rmap(i);

      // Find the corresponding element in Basis1 with swapped indices
      for (ProductBasis<BasisList>::const_iterator I = Basis_21.begin(x.second, x.first); I != Basis_21.end(x.second, x.first); ++I)
      {
         if (Basis_21[*I] == Basis_12[i])
         {
            // This phase factor works for the xxx spin chain
            Result(*I, i) = conj_phase(Basis_12[i], adjoint(Basis_12.Left()[x.first]), Basis_12.Right()[x.second]);
         }
      }
   }
   //TRACE(Result)(Result.Basis1())(Result.Basis2());
   return Result;
}

SimpleOperator
swap_gate_fermion(BasisList const& B1, LinearAlgebra::Vector<double> const& Parity1,
                  BasisList const& B2, LinearAlgebra::Vector<double> const& Parity2,
                  ProductBasis<BasisList, BasisList> const& Basis_21,
                  ProductBasis<BasisList, BasisList> const& Basis_12)
{
   DEBUG_CHECK_EQUAL(Parity1.size(), B1.size());
   DEBUG_CHECK_EQUAL(Parity2.size(), B2.size());
   QuantumNumbers::QuantumNumber Ident(B1.GetSymmetryList());

   SimpleOperator Result(Basis_21.Basis(), Basis_12.Basis(), Ident);
   for (unsigned i = 0; i < Basis_12.size(); ++i)
   {
      std::pair<int,int> x = Basis_12.rmap(i);

      // Find the corresponding element in Basis1 with swapped indices
      for (ProductBasis<BasisList>::const_iterator I = Basis_21.begin(x.second, x.first); I != Basis_21.end(x.second, x.first); ++I)
      {
         if (Basis_21[*I] == Basis_12[i])
         {
            // This phase factor works for the xxx spin chain
            Result(*I, i) = conj_phase(Basis_12[i], adjoint(Basis_12.Left()[x.first]), Basis_12.Right()[x.second])
               * Parity1[x.first] * Parity2[x.second];
         }
      }
   }
   //TRACE(Result)(Result.Basis1())(Result.Basis2());
   return Result;
}
