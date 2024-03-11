// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testmatrixvectorviewmem.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "linearalgebra/matrix.h"
#include "linearalgebra/matrix_utility.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M = random_matrix<double>(10, 15);

   TRACE(M);
   TRACE(vector_view(M));
   TRACE(vector_view(M)[0]);

   TRACE(vector_view(transpose(M))[0]);

   TRACE(vector_view(swap_sort_order(M))[0]);

   //typedef VectorView<SwapSortOrder<Matrix<double>&>::result_type>::result_type X;
   //TRACE(transform(vector_view(swap_sort_order(M)), NormFrob<X>()));

   Vector<Vector<double> > x = vector_view(swap_sort_order(M));

   TRACE(x);
   TRACE(transform(x, NormFrob<Vector<double> >()));
}
