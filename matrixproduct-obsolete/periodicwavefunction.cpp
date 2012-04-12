// -*- C++ -*- $Id$

#include "periodicwavefunction.h"

PeriodicWavefunction::PeriodicWavefunction(int Size, MPStateComponent const& Data)
   : Size_(Size), Q_(MatrixOperator::make_identity(Data.Basis1())), Data_(Data)
{
   PRECONDITION_EQUAL(Data.Basis1(), Data.Basis2());
}

PeriodicWavefunction::PeriodicWavefunction(int Size, 
					   MatrixOperator const& Q, 
					   MPStateComponent const& Data)
   : Size_(Size), Q_(Q), Data_(Data)
{
   PRECONDITION_EQUAL(Data.Basis1(), Data.Basis2());
   PRECONDITION_EQUAL(Q.Basis1(), Q.Basis2());
   PRECONDITION_EQUAL(Q.Basis1(), Data.Basis1());
}

PeriodicWavefunction::PeriodicWavefunction(int Size, 
					   BasisList const& LocalBasis, 
					   VectorBasis const& Basis,
					   QuantumNumber const& Trans)
   : Size_(Size), Q_(Basis, Trans), Data_(LocalBasis, Basis, Basis)
{
}

PeriodicWavefunction::PeriodicWavefunction(int Size, 
					   BasisList const& LocalBasis, 
					   VectorBasis const& Basis)
   : Size_(Size), Q_(MatrixOperator::make_identity(Basis)), Data_(LocalBasis, Basis, Basis)
{
}

std::ostream& operator<<(std::ostream& out, 
			 PeriodicWavefunction const& psi)
{
   out << "Periodic operator, size = " << psi.size() << "\nTransforms as " << psi.TransformsAs()
       << "\nQ: " << psi.Q() << "\nMatrix: " << psi.Data() << '\n';
   return out;
}
