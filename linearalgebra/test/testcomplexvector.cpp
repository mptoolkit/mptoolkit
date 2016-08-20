// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testcomplexvector.cpp
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

#include "linearalgebra/vector.h"
#include "linearalgebra/vectortransform.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/vectoroperations.h"
#include "linearalgebra/stdvector.h"
#include "linearalgebra/slice.h"
#include "linearalgebra/complexvector.h"
#include "linearalgebra/index.h"

using namespace LinearAlgebra;

int main()
{
   typedef std::complex<double> complex;

#if 0
   typedef std::vector<double> vecdouble;
   typedef std::vector<complex> veccomplex;
#else
   typedef LinearAlgebra::Vector<double> vecdouble;
   typedef LinearAlgebra::Vector<complex> veccomplex;
#endif

   vecdouble v1 = Slice(1, 50, 0);
   vecdouble v2 = Slice(1, 50, 0);

   ComplexVector<Vector<double> > v(v1, v2);

   TRACE(veccomplex(v));
   TRACE(v);

   TRACE(norm_1(v));

   TRACE(inner_prod(v,v));
   TRACE(inner_prod(v,conj(v)));
   TRACE(inner_prod(conj(v),conj(v)));

   veccomplex w(v);
   TRACE(w);
   TRACE(inner_prod(real(w),real(w)));
   TRACE(inner_prod(conj(w),conj(w)));
   TRACE(inner_prod(w,w));
   TRACE(inner_prod(conj(w),w));
   TRACE(inner_prod(w,conj(w)));


   TRACE(w[0]);
   TRACE(inner_prod(w[0],w[0]));

   TRACE(norm_2(w));

   TRACE(typeid(real(w)).name());
   real(w) *= 2;
   TRACE(real(w));
   TRACE(w);

   veccomplex v3 = v;

   std::vector<int> Ind(3);
   Ind[0] = 1;
   Ind[1] = 0;
   Ind[2] = 2;


   TRACE(v3);
   TRACE(index(v3,Ind));

   index(v3,Ind) = Ind;
   TRACE(v3);

   TRACE(LinearAlgebra::equal(v3,index(v3,Ind)));
   TRACE(LinearAlgebra::equal(index(v3,Ind), index(v3,Ind)));

   v1 *= 5;
   v2 *= -1;
   real(v) *= 100;
   TRACE(v)(v1)(v2);
   make_complex(v1,v2) = v;
   TRACE(v)(v1)(v2);

}
