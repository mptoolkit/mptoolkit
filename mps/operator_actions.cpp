// -*- C++ -*- $Id$

#include "operator_actions.h"
#include "mpo/operator_component.h"

MatrixOperator apply_right(MatrixOperator const& m, 
			   MPOperator const& Op, 
			   LinearWavefunction const& Psi)
{
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   MatrixOperator Result = m;
   StateComponent E(Op.Basis1(), m.Basis1(), m.Basis2());
   E[0] = m;
   LinearWavefunction::const_iterator I = Psi.begin();
   MPOperator::const_iterator OpIter = Op.begin();
   while (I != Psi.end())
   {
      E = operator_prod(herm(*OpIter), herm(*I), E, *I);
      ++I; ++OpIter;
   }
   return E[0];
}


MatrixOperator apply_left(MatrixOperator const& m, 
			  MPOperator const& Op, 
			  LinearWavefunction const& Psi)
{
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   StateComponent E(Op.Basis2(), m.Basis1(), m.Basis2());
   E[0] = m;
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I = Psi.end();
   MPOperator::const_iterator OpIter = Op.end();
   while (I != Psi.begin())
   {
      --I; --OpIter;
      E = operator_prod(*OpIter, *I, E, herm(*I));
   }
   return E[0];
}

MatrixOperator 
apply_left_qshift(MatrixOperator const& m, 
		  MPOperator const& Op, 
		  LinearWavefunction const& Psi,
		  QuantumNumber const& QShift)
{
   return apply_left(delta_shift(m, adjoint(QShift)), Op, Psi);
}

StateComponent 
apply_right(StateComponent const& In, 
            LinearWavefunction const& Psi1, 
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   StateComponent Result = In;
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   MPOperator::const_iterator OpIter = Op.begin();

   while (OpIter != Op.end())
   {
      Result = operator_prod(herm(*OpIter), herm(*I1), Result, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return Result;
}

MatrixOperator
apply_right_qshift(MatrixOperator const& m, 
		   MPOperator const& Op, 
		   LinearWavefunction const& Psi,
		   QuantumNumber const& QShift)
{
   return delta_shift(apply_right(m, Op, Psi), QShift);
}
