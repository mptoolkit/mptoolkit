// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/test/testtriangularaddition.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/triangularoperator.h"
#include "models/spin-u1.h"

int main()
{
   SiteBlock Site = CreateU1SpinSite(0.5);

   MpOpTriangular SpSm = TriangularTwoSite(Site["Sp"], Site["Sm"]);
   MpOpTriangular SmSp = TriangularTwoSite(Site["Sm"], Site["Sp"]);
   MpOpTriangular SzSz = TriangularTwoSite(Site["Sz"], Site["Sz"]);

   MpOpTriangular SS = 0.5*(SpSm + SmSp) + SzSz;

   TRACE(SpSm.data())(SS.data());
}
