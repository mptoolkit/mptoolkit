// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testvectorcomplex.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "linearalgebra/vector.h"
#include "linearalgebra/vectortransform.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/vectoroperations.h"
#include "linearalgebra/stdvector.h"
#include "linearalgebra/slice.h"

using namespace LinearAlgebra;

int main()
{
   typedef std::complex<double> complex;

#if 1
   typedef std::vector<double> vecdouble;
   typedef std::vector<complex> veccomplex;
#else
   typedef LinearAlgebra::Vector<double> vecdouble;
   typedef LinearAlgebra::Vector<complex> veccomplex;
#endif

   veccomplex v1(3, 0.7);
   CHECK_EQUAL(v1.size(), 3);
   CHECK_EQUAL(v1[0], 0.7);
   CHECK_EQUAL(v1[1], 0.7);
   CHECK_EQUAL(v1[2], 0.7);

   veccomplex v2(3);
   CHECK_EQUAL(v2.size(), 3);
   v2[0] = complex(0,1);
   v2[1] = complex(1,1);
   v2[2] = complex(7,-1);

   vecdouble v3(3);
   CHECK_EQUAL(v3.size(), 3);
   v3[0] = 4;
   v3[1] = 2;
   v3[2] = 8;

   // vector transform as an l-value
   imag(v2) = v3;

   CHECK_EQUAL(v2[0], complex(0,4));
   CHECK_EQUAL(v2[1], complex(1,2));
   CHECK_EQUAL(v2[2], complex(7,8));

   real(v2) = v3;
   CHECK_EQUAL(v2[0], complex(4,4));
   CHECK_EQUAL(v2[1], complex(2,2));
   CHECK_EQUAL(v2[2], complex(8,8));

   TRACE(typeid(real(v2)).name());
}
