// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// wavefunction/operator_actions.cpp
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "operator_actions.h"
#include "mpo/operator_component.h"

MatrixOperator
inject_left(MatrixOperator const& m,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   DEBUG_CHECK_EQUAL(m.Basis1(), Psi1.Basis1());
   DEBUG_CHECK_EQUAL(m.Basis2(), Psi2.Basis1());
   if (Op.is_null())
      return MatrixOperator();
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   MatrixOperator Result = m;
   StateComponent E(Op.Basis1(), m.Basis1(), m.Basis2());
   E[0] = m;
   E.debug_check_structure();
   E = inject_left(E, Psi1, Op, Psi2);
   return E[0];
}


MatrixOperator
inject_left(MatrixOperator const& m,
            GenericMPO const& Op,
            LinearWavefunction const& Psi)
{
   return inject_left(m, Psi, Op, Psi);
}

MatrixOperator
inject_left_qshift(MatrixOperator const& m,
                   GenericMPO const& Op,
                   LinearWavefunction const& Psi,
                   QuantumNumber const& QShift)
{
   return inject_left(delta_shift(m, QShift), Op, Psi);
}

StateComponent
inject_left(StateComponent const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   StateComponent Result = In;
   In.debug_check_structure();
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   GenericMPO::const_iterator OpIter = Op.begin();

   while (OpIter != Op.end())
   {
      Result = contract_from_left(*OpIter, herm(*I1), Result, *I2);
      ++I1; ++I2; ++OpIter;
   }
   DEBUG_CHECK(I1 == Psi1.end());
   DEBUG_CHECK(I2 == Psi2.end());
   return Result;
}

MatrixOperator
inject_left(MatrixOperator const& m,
            LinearWavefunction const& Psi1,
            QuantumNumbers::QuantumNumber const& QShift,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   DEBUG_CHECK_EQUAL(m.Basis1(), Psi1.Basis2());
   DEBUG_CHECK_EQUAL(m.Basis2(), Psi2.Basis2());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   MatrixOperator mm = delta_shift(m, QShift);
   StateComponent E(Op.Basis1(), mm.Basis1(), mm.Basis2());
   E[0] = mm;
   E.debug_check_structure();
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   GenericMPO::const_iterator OpIter = Op.begin();
   while (OpIter != Op.end())
   {
      if (I1 == Psi1.end())
      {
         I1 = Psi1.begin();
         I2 = Psi2.begin();
         E = delta_shift(E, QShift);
      }
      E = contract_from_left(*OpIter, herm(*I1), E, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return E[0]; //delta_shift(E[0], QShift);
}

MatrixOperator
inject_right(MatrixOperator const& m,
             GenericMPO const& Op,
             LinearWavefunction const& Psi)
{
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   PRECONDITION_EQUAL(Psi.size(), Op.size());
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   StateComponent E(Op.Basis2(), m.Basis1(), m.Basis2());
   E[0] = m;
   E.debug_check_structure();
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I = Psi.end();
   GenericMPO::const_iterator OpIter = Op.end();
   while (I != Psi.begin())
   {
      --I; --OpIter;
      E = contract_from_right(herm(*OpIter), *I, E, herm(*I));
   }
   return E[0];
}

MatrixOperator
inject_right(MatrixOperator const& m,
             LinearWavefunction const& Psi1,
             GenericMPO const& Op,
             LinearWavefunction const& Psi2)
{
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   StateComponent E(Op.Basis2(), m.Basis1(), m.Basis2());
   E[0] = m;
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   GenericMPO::const_iterator OpIter = Op.end();
   while (I1 != Psi1.begin())
   {
      --I1; --I2; --OpIter;
      E = contract_from_right(herm(*OpIter), *I1, E, herm(*I2));
   }
   return E[0];
}

MatrixOperator
inject_right_qshift(MatrixOperator const& m,
                    GenericMPO const& Op,
                    LinearWavefunction const& Psi,
                    QuantumNumber const& QShift)
{
   return delta_shift(inject_right(m, Op, Psi), adjoint(QShift));
}

MatrixOperator
inject_right(MatrixOperator const& m,
             LinearWavefunction const& Psi1,
             QuantumNumbers::QuantumNumber const& QShift,
             GenericMPO const& Op,
             LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   StateComponent E(Op.Basis2(), m.Basis1(), m.Basis2());
   E[0] = m;
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   GenericMPO::const_iterator OpIter = Op.end();
   while (OpIter != Op.begin())
   {
      if (I1 == Psi1.begin())
      {
         I1 = Psi1.end();
         I2 = Psi2.end();
         E = delta_shift(E, adjoint(QShift));
      }
      --I1; --I2; --OpIter;
      E = contract_from_right(herm(*OpIter), *I1, E, herm(*I2));
   }
   return delta_shift(E[0], adjoint(QShift));
}

StateComponent
contract_from_left_mask(OperatorComponent const& M,
                        HermitianProxy<StateComponent> const& A,
                        StateComponent const& E,
                        StateComponent const& B,
                        std::vector<int> const& Mask1,
                        std::vector<int> const& Mask2)
{
   StateComponent Result(M.Basis2(), A.base().Basis2(), B.Basis2());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M); I; ++I)
   {
      // skip over masked components
      if (!Mask1[I.index()])
         continue;

      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // skip over masked components
         if (!Mask2[J.index2()])
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

StateComponent
inject_left_mask(StateComponent const& In,
                 LinearWavefunction const& Psi1,
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

   StateComponent Result(In);

   while (OpIter != Op.end())
   {
      Result = contract_from_left_mask(*OpIter, herm(*I1), Result, *I2, *MaskIter, *(MaskIter+1));
      ++I1; ++I2; ++OpIter; ++MaskIter;
   }
   return Result;
}

StateComponent
inject_right_mask(StateComponent const& In,
                 LinearWavefunction const& Psi1,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());

   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   GenericMPO::const_iterator OpIter = Op.end();
   std::vector<std::vector<int> >::const_iterator MaskIter = Mask.end();

   StateComponent E;
   StateComponent Result(In);

   while (OpIter != Op.begin())
   {
      std::swap(E, Result);

      --I1; --I2; --OpIter; --MaskIter;

      Result = contract_from_right_mask(herm(*OpIter), *I1, E, herm(*I2), *(MaskIter-1), *MaskIter);

   }
   return Result;
}
