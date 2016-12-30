// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testfill.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(10,10);

   fill(M, 2.0);
   CHECK_CLOSE(norm_frob(M), 20);

   Matrix<std::complex<double> > X(10, 10);
   fill(real(X), 3.0);
   fill(imag(X), 4.0);
   CHECK_CLOSE(norm_frob(X), 50);

   Vector<double> V(100);
   fill(V, 2.0);
   CHECK_CLOSE(norm_frob(V), 20);

   Vector<std::complex<double> > W(100);
   fill(real(W), 3.0);
   fill(imag(W), 4.0);
   CHECK_CLOSE(norm_frob(W), 50);
}
