// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/test/testtriangularaddition.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
