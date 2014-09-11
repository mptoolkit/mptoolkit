// -*- C++ -*- $Id$

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

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   MatrixOperator Result = m;
   StateComponent E(Op.Basis1(), m.Basis1(), m.Basis2());
   E[0] = m;
   E.debug_check_structure();
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   GenericMPO::const_iterator OpIter = Op.begin();
   while (I1 != Psi1.end())
   {
      E = operator_prod(herm(*OpIter), herm(*I1), E, *I2);
      ++I1; ++I2; ++OpIter;
   }
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
      Result = operator_prod(herm(*OpIter), herm(*I1), Result, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return Result;
}

MatrixOperator
inject_right(MatrixOperator const& m, 
             GenericMPO const& Op, 
             LinearWavefunction const& Psi)
{
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
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
      E = operator_prod(*OpIter, *I, E, herm(*I));
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
      E = operator_prod(*OpIter, *I1, E, herm(*I2));
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

bool
is_local_basis_compatible(LinearWavefunction const& Psi, GenericMPO const& M)
{
   // If the MPO has a different support to Psi then we are not compatible
   if (Psi.size() != M.size())
      return false;

   // Also check the individual local basis at each step
   LinearWavefunction::const_iterator iPsi = Psi.begin();
   GenericMPO::const_iterator iM = M.begin();
   while (iPsi != Psi.end())
   {
      if (iPsi->LocalBasis() != iM->LocalBasis1()
          || iPsi->LocalBasis() != iM->LocalBasis2())
         return false;

      ++iPsi; ++iM;
   }
   return true;
}

void
local_basis_compatible_or_abort(LinearWavefunction const& Psi, GenericMPO const& M)
{
   CHECK_EQUAL(Psi.size(), M.size())("Wavefunction and MPO must have the same size unit cell!");

   int SiteNumber = 0;
   LinearWavefunction::const_iterator iPsi = Psi.begin();
   GenericMPO::const_iterator iM = M.begin();
   while (iPsi != Psi.end())
   {
      CHECK_EQUAL(iPsi->LocalBasis(), iM->LocalBasis1())
         ("Local basis for wavefunction does not match local basis 1 of the MPO at a site!")(SiteNumber);

      CHECK_EQUAL(iPsi->LocalBasis(), iM->LocalBasis2())
         ("Local basis for wavefunction does not match local basis 2 of the MPO at a site!")(SiteNumber);
      
      ++iPsi; ++iM; ++SiteNumber;
   }
}
