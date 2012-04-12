// -*- C++ -*- $Id$

#include "operator_actions.h"
#include "mpo/operator_component.h"

MatrixOperator transfer_from_left(MatrixOperator const& m, 
                                  LinearOperator const& Op, 
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
   LinearOperator::const_iterator OpIter = Op.begin();
   while (I != Psi.end())
   {
      E = operator_prod(herm(*OpIter), herm(*I), E, *I);
      ++I; ++OpIter;
   }
   return E[0];
}


MatrixOperator transfer_from_right(MatrixOperator const& m, 
                                   LinearOperator const& Op, 
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
   LinearOperator::const_iterator OpIter = Op.end();
   while (I != Psi.begin())
   {
      --I; --OpIter;
      E = operator_prod(*OpIter, *I, E, herm(*I));
   }
   return E[0];
}
