// -*- C++ -*- $Id$

#include "momentum_operations.h"
#include "operator_actions.h"

PolynomialType
delta_shift(PolynomialType const& In, QuantumNumber const& QShift)
{
   PolynomialType Result(In);
   for (PolynomialType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

// do finite momentum at the same time?
std::vector<PolynomialType>
inject_left(std::vector<PolynomialType> const& In, 
            LinearWavefunction const& Psi1, 
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());

   std::vector<PolynomialType> Result(Op.Basis2().size());
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

PolynomialType
inject_left(PolynomialType const& In, 
            LinearWavefunction const& Psi1, 
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   std::vector<PolynomialType> Vec(1, In);
   Vec = inject_left(Vec, Psi1, Op, Psi2);
   CHECK_EQUAL(Vec.size(), 1);
   return Vec[0];
}

// Calculates result' = C[column] = sum_{j > Column} E[j] * T_Op(j, Column)
// Assumes that E[j] is defined, for j > Column
PolynomialType
MultiplyLeft(std::vector<PolynomialType> const& E, 
             TriangularMPO const& Op, 
             LinearWavefunction const& Psi, 
             QuantumNumber const& QShift, int Column)
{
   CHECK_EQUAL(Op.size(), Psi.size());
   PolynomialType Result;

   // replace this stuff with the inject_left implementation, and extract_column() in TriangularOperator

   GenericMPO OpCol = extract_column(Op, Column);
   std::vector<PolynomialType> C = inject_left(E, Psi, OpCol, Psi);
   return delta_shift(C[Column], QShift);
}

Polynomial<std::complex<double> >
ExtractOverlap(Polynomial<MatrixOperator> const& E, MatrixOperator const& Rho)
{
   Polynomial<std::complex<double> > Overlap;
   for (Polynomial<MatrixOperator>::const_iterator I = E.begin(); I != E.end(); ++I)
   {
      Overlap[I->first] = inner_prod(I->second, Rho);
   }
   return Overlap;
}

//
// Momentum
//

MomentumPolynomialType
delta_shift(MomentumPolynomialType const& In, QuantumNumber const& QShift)
{
   MomentumPolynomialType Result(In);
   for (MomentumPolynomialType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

std::vector<MomentumPolynomialType>
delta_shift(std::vector<MomentumPolynomialType> const& In, QuantumNumber const& QShift)
{
   std::vector<MomentumPolynomialType> Result(In);
   for (std::vector<MomentumPolynomialType>::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = delta_shift(*I, QShift);
   }
   return Result;
}

void
add_triple_prod(MomentumPolynomialType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                MomentumPolynomialType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp)
{
   // loop over momenta
   for (MomentumPolynomialType::const_iterator K = E.begin(); K != E.end(); ++K)
   {
      // loop over degrees of the polynomial
      for (PolynomialType::const_iterator D = K->second.begin(); D != K->second.end(); ++D)
      {
         Result[K->first][D->first] += Factor * triple_prod(x, D->second, y, qxy, qEp);
      }
   }
}

std::vector<MomentumPolynomialType>
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              std::vector<MomentumPolynomialType> const& E, 
              StateComponent const& B)
{
   //   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().Basis1());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1().size(), E.size());
   
   std::vector<MomentumPolynomialType> Result(M.base().Basis2().size());

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
                  add_triple_prod(Result[J.index2()], herm(*S), 
                                  herm(A.base()[S.index1()]), 
                                  E[J.index1()], 
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.base().Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<MomentumPolynomialType>
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              std::vector<MomentumPolynomialType> const& E, 
              StateComponent const& B,
              std::vector<int> const& OutMask,
              std::vector<int> const& InMask)
{
   std::vector<MomentumPolynomialType> Result(M.base().Basis2().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // skip over masked components
      if (!InMask[I.index()])
         continue;

      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // skip over masked components
         if (!OutMask[J.index2()])
            continue;

         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index2()], herm(*S), 
                                  herm(A.base()[S.index1()]), 
                                  E[J.index1()], 
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.base().Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<MomentumPolynomialType>
inject_left(std::vector<MomentumPolynomialType> const& In, 
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

   std::vector<MomentumPolynomialType> E;
   std::vector<MomentumPolynomialType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      //      std::vector<MomentumPolynomialType>(OpIter->Basis2().size()).swap(Result);

      Result = operator_prod(herm(*OpIter), herm(*I1), E, *I2);

      ++I1; ++I2; ++OpIter;
   }
   return Result;
}

std::vector<MomentumPolynomialType>
inject_left_mask(std::vector<MomentumPolynomialType> const& In, 
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
   std::vector<std::vector<int> >::const_iterator MaskIter = Mask.begin();

   std::vector<MomentumPolynomialType> E;
   std::vector<MomentumPolynomialType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      //      std::vector<MomentumPolynomialType>(OpIter->Basis2().size()).swap(Result);

      Result = operator_prod(herm(*OpIter), herm(*I1), E, *I2, *(MaskIter+1), *MaskIter);

      ++I1; ++I2; ++OpIter; ++MaskIter;
   }
   return delta_shift(Result, QShift);
}
