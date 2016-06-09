// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testvectortransform.cpp
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

int main()
{
   LinearAlgebra::Vector<double> v1;   
   CHECK_EQUAL(v1.size(), 0);
   CHECK_EQUAL(size(v1), 0);

   LinearAlgebra::Vector<double> v2(3);
   CHECK_EQUAL(v2.size(), 3);
   CHECK_EQUAL(size(v2), 3);

   LinearAlgebra::Vector<double> v3(3, 0.7);
   CHECK_EQUAL(v3.size(), 3);
   CHECK_EQUAL(size(v3), 3);
   CHECK_EQUAL(v3[0], 0.7);
   CHECK_EQUAL(v3[1], 0.7);
   CHECK_EQUAL(v3[2], 0.7);

   LinearAlgebra::Vector<double> v4(-v3);
   CHECK_EQUAL(v4.size(), 3);
   CHECK_EQUAL(size(v4), 3);
   CHECK_EQUAL(v4[0], -0.7);
   CHECK_EQUAL(v4[1], -0.7);
   CHECK_EQUAL(v4[2], -0.7);

   // check that negating twice gets us back to the same type,
   // so the double negation is absorbed.
   CHECK(typeid(transform(transform(v3, LinearAlgebra::Negate<double>()), 
			  LinearAlgebra::Negate<double>())) 
	 == typeid(v3));

   // operator== for rhs proxy
   CHECK_EQUAL(v4, -v3);

   // sanity check
   CHECK(v4 != v3)(v4)(v3);
   
   // real, imag, conj, transpose, herm should all be identity transformations
   // for real vectors
   CHECK_EQUAL(real(v3), -v4);
   
   CHECK(typeid(real(v3)) == typeid(v3));
   CHECK(typeid(imag(v3)) == typeid(v3));
   CHECK(typeid(conj(v3)) == typeid(v3));
   CHECK(typeid(transpose(v3)) == typeid(v3));
   CHECK(typeid(herm(v3)) == typeid(v3));

   // in addition, they should be the same object
   CHECK(&real(v3) == &v3)(&real(v3))(&v3);
   CHECK(&imag(v3) == &v3)(&imag(v3))(&v3);
   CHECK(&conj(v3) == &v3)(&conj(v3))(&v3);
   CHECK(&transpose(v3) == &v3)(&transpose(v3))(&v3);
   CHECK(&herm(v3) == &v3)(&herm(v3))(&v3);

   std::cout << v3 << std::endl;

   std::cout << transform(v3, 
      LinearAlgebra::bind_second(LinearAlgebra::Multiplication<double, double>(), 3.0))
	     << std::endl;

   //   Vector<double> v4(v3);

   LinearAlgebra::Vector<double> v6(transform(v3, 
			       LinearAlgebra::bind_second(LinearAlgebra::
							  Multiplication<double, double>(), 3.0)));

   std::cout << v6 << std::endl;

   LinearAlgebra::Vector<double> v5;

   v5 = transform(v3, 
      LinearAlgebra::bind_second(LinearAlgebra::Multiplication<double, double>(), 3.0));

   std::cout << v5 << std::endl;

   std::cout << transform(v3, LinearAlgebra::Negate<double>()) << std::endl;

   std::cout << typeid(transform(transform(v3, LinearAlgebra::Negate<double>()), LinearAlgebra::Negate<double>())).name() << std::endl;


   std::cout << (-v4) << std::endl;

   std::cout << typeid(-v4).name() << std::endl;

   std::cout << transpose(v4) << std::endl;

   std::cout << (v4 * 7.5) << std::endl;

   std::cout << (v4 * scalar(v4) ) << std::endl;

   std::cout << LinearAlgebra::VectorScalarMultiplication<LinearAlgebra::Vector<double>, int>()(v4, 7) << std::endl;

   std::cout << norm_2(0.1 * v4) << std::endl;
   std::cout << norm_1(0.2 * v4) << std::endl;
   std::cout << norm_inf(0.1 * v4) << std::endl;


   std::cout << norm_2(v4 * scalar(v4)) << std::endl;

}
