// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// wavefunction/momentum_operations.cpp
//
// Copyright (C) 2014-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "momentum_operations.h"
#include "operator_actions.h"

void
ScalePoly(KMatrixPolyType& K, std::complex<double> x)
{
   for (auto& k : K)
   {
      k.second *= x;
   }
}

MatrixPolyType
delta_shift(MatrixPolyType const& In, QuantumNumber const& QShift)
{
   MatrixPolyType Result(In);
   for (MatrixPolyType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

// do finite momentum at the same time?
std::vector<MatrixPolyType>
inject_left(std::vector<MatrixPolyType> const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());

   std::vector<MatrixPolyType> Result(Op.Basis2().size());
   int MaxDegree = 0;
   for (unsigned i = 0; i < In.size(); ++i)
      MaxDegree = std::max(In[i].degree(), MaxDegree);

   for (int Degree = 0; Degree <= MaxDegree; ++Degree)
   {
      StateComponent E(Op.Basis1(), Psi1.Basis1(), Psi2.Basis1());
      for (unsigned k = 0; k < E.size(); ++k)
      {
         E[k] = In[k][Degree];
      }

      E = inject_left(E, Psi1, Op, Psi2);

      CHECK_EQUAL(E.size(), Result.size());
      for (unsigned i = 0; i < Result.size(); ++i)
      {
         Result[i][Degree] += E[i];
      }
   }
   return Result;
}

MatrixPolyType
inject_left(MatrixPolyType const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   std::vector<MatrixPolyType> Vec(1, In);
   Vec = inject_left(Vec, Psi1, Op, Psi2);
   CHECK_EQUAL(Vec.size(), 1);
   return Vec[0];
}

#if 0
// Calculates result' = C[column] = sum_{j > Column} E[j] * T_Op(j, Column)
// Assumes that E[j] is defined, for j > Column
MatrixPolyType
MultiplyLeft(std::vector<MatrixPolyType> const& E,
             BasicTriangularMPO const& Op,
             LinearWavefunction const& Psi,
             QuantumNumber const& QShift, int Column)
{
   CHECK_EQUAL(Op.size(), Psi.size());
   MatrixPolyType Result;

   // replace this stuff with the inject_left implementation, and extract_column() in TriangularOperator

   GenericMPO OpCol = extract_column(Op, Column);
   std::vector<MatrixPolyType> C = inject_left(E, Psi, OpCol, Psi);
   return delta_shift(C[Column], QShift);
}
#endif

// Here the E-matrix is assumed to not require conjugation, so
// the result is implented as inner_prod(E, Rho)
ComplexPolyType
ExtractOverlap(MatrixPolyType const& E, MatrixOperator const& Rho)
{
   ComplexPolyType Overlap;
   for (auto const& x : E)
   {
      if (x.second.TransformsAs() == Rho.TransformsAs())
	 Overlap[x.first] = inner_prod(Rho, x.second);
   }
   return Overlap;
}

KComplexPolyType
ExtractOverlap(KMatrixPolyType const& E, MatrixOperator const& Rho)
{
   KComplexPolyType Result;
   for (auto const& x : E)
      Result[x.first] = ExtractOverlap(x.second, Rho);
   return Result;
}

//
// Momentum
//

KMatrixPolyType
delta_shift(KMatrixPolyType const& In, QuantumNumber const& QShift)
{
   KMatrixPolyType Result(In);
   for (KMatrixPolyType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

std::vector<MatrixPolyType>
delta_shift(std::vector<MatrixPolyType> const& In, QuantumNumber const& QShift)
{
   std::vector<MatrixPolyType> Result(In);
   for (std::vector<MatrixPolyType>::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = delta_shift(*I, QShift);
   }
   return Result;
}

std::vector<KMatrixPolyType>
delta_shift(std::vector<KMatrixPolyType> const& In, QuantumNumber const& QShift)
{
   std::vector<KMatrixPolyType> Result(In);
   for (std::vector<KMatrixPolyType>::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = delta_shift(*I, QShift);
   }
   return Result;
}

void
add_triple_prod(MatrixPolyType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                MatrixPolyType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp)
{
   // loop over degrees of the polynomial
   for (MatrixPolyType::const_iterator D = E.begin(); D != E.end(); ++D)
   {
       add_triple_prod(Result[D->first], Factor, x, D->second, y, qxy, qEp);
   }
}

void
add_triple_prod(KMatrixPolyType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                KMatrixPolyType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp)
{
   // loop over momenta
   for (KMatrixPolyType::const_iterator K = E.begin(); K != E.end(); ++K)
   {
      add_triple_prod(Result[K->first], Factor, x, K->second, y, qxy, qEp);
   }
}

void
add_triple_prod(MatrixPolyType& Result, std::complex<double> Factor,
              MatrixOperator const& x,
              MatrixPolyType const& E,
              HermitianProxy<MatrixOperator> const& y,
              QuantumNumber const& qxy,
              QuantumNumber const& qEp)
{
   // loop over degrees of the polynomial
   for (MatrixPolyType::const_iterator D = E.begin(); D != E.end(); ++D)
   {
      add_triple_prod(Result[D->first], Factor, x, D->second, y, qxy, qEp);
   }
}

void
add_triple_prod(KMatrixPolyType& Result, std::complex<double> Factor,
              MatrixOperator const& x,
              KMatrixPolyType const& E,
              HermitianProxy<MatrixOperator> const& y,
              QuantumNumber const& qxy,
              QuantumNumber const& qEp)
{
   // loop over degrees of the polynomial
   for (KMatrixPolyType::const_iterator K = E.begin(); K != E.end(); ++K)
   {
      add_triple_prod(Result[K->first], Factor, x, K->second, y, qxy, qEp);
   }
}


std::vector<KMatrixPolyType>
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   std::vector<KMatrixPolyType> const& E,
                   StateComponent const& B)
{
   //   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().Basis1());
   DEBUG_PRECONDITION_EQUAL(M.LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis1().size(), E.size());

   std::vector<KMatrixPolyType> Result(M.Basis2().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index2()], *S,
                                  herm(A.base()[S.index1()]),
                                  E[J.index1()],
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<KMatrixPolyType>
contract_from_right(HermitianProxy<OperatorComponent> const& M,
                    StateComponent const& A,
                    std::vector<KMatrixPolyType> const& F,
                    HermitianProxy<StateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2().size(), F.size());

   std::vector<KMatrixPolyType> Result(M.base().Basis1().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index1()], herm(*S),
                                  A[S.index1()],
                                  F[J.index2()],
                                  herm(B.base()[S.index2()]),
                                  k->TransformsAs(),
                                  M.base().Basis1()[J.index1()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<KMatrixPolyType>
inject_left(std::vector<KMatrixPolyType> const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   PRECONDITION_EQUAL(Op.Basis1().size(), In.size());

   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   GenericMPO::const_iterator OpIter = Op.begin();

   std::vector<KMatrixPolyType> E;
   std::vector<KMatrixPolyType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      //      std::vector<KMatrixPolyType>(OpIter->Basis2().size()).swap(Result);

      Result = contract_from_left(*OpIter, herm(*I1), E, *I2);

      ++I1; ++I2; ++OpIter;
   }
   return Result;
}

std::vector<MatrixPolyType>
inject_left_mask(std::vector<MatrixPolyType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   //   PRECONDITION_EQUAL(Op.Basis1().size(), In.size());

   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   GenericMPO::const_iterator OpIter = Op.begin();
   std::vector<std::vector<int>>::const_iterator MaskIter = Mask.begin();

   std::vector<MatrixPolyType> E;
   std::vector<MatrixPolyType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      Result = contract_from_left_mask(*OpIter, herm(*I1), E, *I2, *(MaskIter+1), *MaskIter);

      ++I1; ++I2; ++OpIter; ++MaskIter;
   }
   return delta_shift(Result, QShift);
}

std::vector<KMatrixPolyType>
inject_left_mask(std::vector<KMatrixPolyType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   //   PRECONDITION_EQUAL(Op.Basis1().size(), In.size());

   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   GenericMPO::const_iterator OpIter = Op.begin();
   std::vector<std::vector<int>>::const_iterator MaskIter = Mask.begin();

   std::vector<KMatrixPolyType> E;
   std::vector<KMatrixPolyType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      Result = contract_from_left_mask(*OpIter, herm(*I1), E, *I2, *(MaskIter+1), *MaskIter);

      ++I1; ++I2; ++OpIter; ++MaskIter;
   }
   return delta_shift(Result, QShift);
}

std::vector<MatrixPolyType>
inject_right_mask(std::vector<MatrixPolyType> const& In,
                  LinearWavefunction const& Psi1,
                  QuantumNumber const& QShift,
                  GenericMPO const& Op,
                  LinearWavefunction const& Psi2,
                  std::vector<std::vector<int>> const& Mask)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());

   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   GenericMPO::const_iterator OpIter = Op.end();
   std::vector<std::vector<int> >::const_iterator MaskIter = Mask.end();

   std::vector<MatrixPolyType> F;
   std::vector<MatrixPolyType> Result(In);

   while (OpIter != Op.begin())
   {
      std::swap(F, Result);

      --I1; --I2; --OpIter; --MaskIter;

      Result = contract_from_right_mask(herm(*OpIter), *I1, F, herm(*I2), *(MaskIter-1), *MaskIter);

   }
   return delta_shift(Result, adjoint(QShift));
}

std::vector<KMatrixPolyType>
inject_right_mask(std::vector<KMatrixPolyType> const& In,
                  LinearWavefunction const& Psi1,
                  QuantumNumber const& QShift,
                  GenericMPO const& Op,
                  LinearWavefunction const& Psi2,
                  std::vector<std::vector<int>> const& Mask)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());

   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   GenericMPO::const_iterator OpIter = Op.end();
   std::vector<std::vector<int> >::const_iterator MaskIter = Mask.end();

   std::vector<KMatrixPolyType> F;
   std::vector<KMatrixPolyType> Result(In);

   while (OpIter != Op.begin())
   {
      std::swap(F, Result);

      --I1; --I2; --OpIter; --MaskIter;

      Result = contract_from_right_mask(herm(*OpIter), *I1, F, herm(*I2), *(MaskIter-1), *MaskIter);

   }
   return delta_shift(Result, adjoint(QShift));
}
