// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testsparsematrix.cpp
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

#include "linearalgebra/sparsematrix.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
using tracer::typeid_name;

int main()
{
   //   typedef Vector<MapVector<std::complex<double> > > VV;
   //   typedef CompressedMatrix<VV> MatType;
   typedef SparseMatrix<std::complex<double> > MatType;
   MatType M(3,4);

      M(1,2) = 5;
      M(0,1) = -10;

   TRACE(M);

   TRACE(transpose(M));

   TRACE(get_element(M,1,2));

   TRACE(typeid(get_element(M,1,2)).name());

   get_element(M,1,2) = 5;

   M(2,1) = 7;

   TRACE(M);

   TRACE(-M);

   MatType N(-M);
   TRACE(N);


   TRACE(typeid_name(-N));
   TRACE(typeid_name(- (-N)));

   TRACE(real(N));
   TRACE(transpose(N));

   TRACE(typeid_name(transpose(N)));
   TRACE(herm(N));

   MatType C(M);
   get_element(C,1,1) = std::complex<double>(3,4);
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

   TRACE(transform(C, C*10.0, Addition<std::complex<double>, std::complex<double> >()));

   TRACE(C+C*10.0 + transpose(0.5*transpose(C*4.0)));
   TRACE(typeid_name(C+C*10.0 + transpose(0.5*transpose(C*4.0))));

   TRACE(C*13.0);
   
   TRACE(typeid_name(eval_expression(C+C*10.0 + transpose(0.5*transpose(C*4.0)))));
   TRACE(norm_frob(eval_expression(C+C*10.0 + imag(transpose(0.5*transpose(C*4.0))))));
   TRACE(norm_frob(C+C*10.0 + imag(transpose(0.5*transpose(C*4.0)))));
   TRACE(C+C*10.0 + imag(transpose(0.5*transpose(C*4.0))));
   TRACE(norm_frob(real(C+C*10.0 + imag(transpose(0.5*transpose(C*4.0))))));
   TRACE(norm_frob(imag(C+C*10.0 + imag(transpose(0.5*transpose(C*4.0))))));

   TRACE(C);

   TRACE(C * transpose(C));

   TRACE(C * herm(C));

   TRACE(transpose(C) * C);

   TRACE(herm(C) * C);

   TRACE(herm(C) * (imag(C) + real(C)));
}
