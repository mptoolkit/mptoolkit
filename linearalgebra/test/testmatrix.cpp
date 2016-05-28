// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testmatrix.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/matrix.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
using tracer::typeid_name;

int main()
{
   Matrix<double> M(3,4,0);

#if 1
   M(1,2) = 5;
   M(0,1) = -10;

   TRACE(M);

   TRACE(transpose(M));

   TRACE(-M);

   Matrix<double> N(-M);
   TRACE(N);


   TRACE(typeid_name(-N));
   TRACE(typeid_name(- (-N)));

   TRACE(real(N));
   TRACE(transpose(N));

   TRACE(data(N));
   TRACE(data(transpose(N)));
   TRACE(typeid_name(transpose(N)));
   TRACE(herm(N));

#endif

   Matrix<std::complex<double> > C(M);

#if 1
#if 1
   C(0,0) = std::complex<double>(6,7);
   TRACE(C);
   TRACE(real(C));
   TRACE(transpose(C));
   TRACE(typeid_name(herm(C)));
   TRACE(herm(C));
   TRACE(real(herm(C)));
   TRACE(typeid_name(real(herm(C))));
   TRACE(typeid_name(imag(herm(C))));
   TRACE(typeid_name(-real(herm(C))));
   TRACE(typeid_name(-imag(herm(C))));

   TRACE(typeid_name(imag(transpose(C))));
   TRACE(typeid_name(transpose(imag(C))));
   TRACE(typeid_name(-transpose(imag(C))));

   TRACE(typeid_name(C+C));
   //   TRACE(typeid_name(C+transpose(C)));

   TRACE(C);
   TRACE(C+C);

#endif
   TRACE(real(C)+imag(C));

   real(C) = real(C) + imag(C) * 10;

   TRACE(C);

   TRACE(typeid_name(C+C*10.0 + transpose(0.5*transpose(C*4.0))));

#endif

   TRACE(C+C*10.0 + transpose(0.5*transpose(C*4.0)));

   TRACE(C*13.0);

   TRACE(norm_frob(C+C*10.0 + imag(transpose(0.5*transpose(C*4.0)))));

   TRACE(C * transpose(C));

   TRACE(transpose(C) * C);

   TRACE(real(C) * transpose(real(C)));

}
