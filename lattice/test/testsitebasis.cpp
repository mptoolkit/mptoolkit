// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/test/testsitebasis.cpp
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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

#include "siteoperator/sitebasis.h"
#include "quantumnumbers/su2.h"
#include "quantumnumbers/u1.h"

using QuantumNumbers::U1;
using QuantumNumbers::SU2;
using QuantumNumbers::QNConstructor;

int main()
{
   SymmetryList Symmetry("N:U(1),S:SU(2)");

   QNConstructor<U1,SU2> QN(Symmetry);

   SiteBasis HubbardBasis(Symmetry);
   HubbardBasis.push_back("Empty",  QN(0, 0));
   HubbardBasis.push_back("Single", QN(1, 0.5));
   HubbardBasis.push_back("Double", QN(2, 0));

   CHECK_EQUAL(HubbardBasis.size(), 3);

   CHECK_EQUAL(HubbardBasis[0].second, QN(0,0));
   CHECK_EQUAL(HubbardBasis[1].second, QN(1,0.5));
   CHECK_EQUAL(HubbardBasis[2].second, QN(2,0));

   CHECK_EQUAL(HubbardBasis[0].first, "Empty");
   CHECK_EQUAL(HubbardBasis[1].first, "Single");
   CHECK_EQUAL(HubbardBasis[2].first, "Double");
}
