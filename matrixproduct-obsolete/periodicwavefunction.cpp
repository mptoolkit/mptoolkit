// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/periodicwavefunction.cpp
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
