// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/examples/example1.cpp
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
//
// example program for vectors
//
// compile with
// g++ -DNDEBUG -I. -I$BOOST_DIRECTORY common/poolallocator.cpp common/stackallocator.cpp example1.cpp -llapack -lblas -lg2c -o example1
//
// to see where BLAS level 1 calls are made, add
// -DBLAS1_TRACE_DETAILED
// to the compile line.
//

#include "linearalgebra/vector.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>

using namespace LinearAlgebra;
using namespace std;

int main()
{
   Vector<double> V(10);    // vector of length 10, uninitialized values
   Vector<double> U(10, 1.0); // vector length 10, initialized to 1.0

   fill(U,1.0);  // an alternative way to set all elements to 1.0

   TRACE(U);       // here, U = (1,1,1,1,1,1,1,1,1,1)

   Vector<int> FirstHalf  = Range(0, 5);	// range [0, 5), ie. elements 0,1,2,3,4
   Vector<int> SecondHalf = Range(5, 10);	// range [5, 10), ie. elements 5,6,7,8,9

   Vector<int> EvenElements = Slice(0,5,2);	// slice starting at element 0, length 5, stride 2
   Vector<int> OddElements = Slice(1,5,2);	// slice starting at element 1, length 5, stride 2

   TRACE(EvenElements)(U[EvenElements]);

   //U[EvenElements] = Vector<double>(5,10.);	//
   //U[EvenElements]*= 10;			// alternative
   U[EvenElements] = U[EvenElements]*10;	// alternative

   TRACE(U);					// here, U = (10,1,10,1,10,1,10,1,10,1)

   // A Range, Slice or Index is itself a model of VectorExpression<int>, therefore can be used like this:

   U[ OddElements] = FirstHalf;

   TRACE(U);				       // here, U = (10,0,10,1,10,2,10,3,10,4)

   VectorRef<double> VRef(V);			// a reference to V
   VectorRef<double> URef(U);			// a reference to U

   Vector<double> W(U);        // an ordinary vector, initially equal to U (but shares the representation!)

   fill( URef[SecondHalf], -1);			// sets the elements 5,6,7,8,9 of U (and URef) to -1
                                       		// but leaves W unaffected.  (This implies a copy is made).

   TRACE(U)(URef);				// here, U = (10,0,10,1,10,-1,-1,-1,-1,-1)
						 // URef is always equal to U

   TRACE(W);					// here, W = (10,0,10,1,10,2,10,3,10,4)

   VRef = U;					// equivalent to V = U    (or V = URef, or even VRef = URef)


   cout<<"\n\n############## VECTOR PRODUCTS #################################################\n";
   // parallel_prod(X,Y) calculates sum_i X_i * Y_i    (a "shallow" product)
   // NOTE: this is not a real inner product, if the type is complex etc
   TRACE( parallel_prod(U,W));

   // inner_prod calculates recursively: sum_i inner_prod(herm(X_i), Y_i)
   // herm() is hermitian conjugate
   TRACE( inner_prod(U,W));

   // real difference if we have a vector of vectors
   Vector<Vector<double> > VecOfVec(2);  // a vector of Vector<double>, length 2
   VecOfVec[0] = U;
   VecOfVec[1] = W;

   TRACE( VecOfVec );

   TRACE( inner_prod(VecOfVec, VecOfVec) );

   // this should be the same as inner_prod(U,U) + inner_prod(W,W) :
   TRACE( inner_prod(U,U) + inner_prod(W,W) );

   // Note that parallel_prod(VecOfVec, VecOfVec) won't compile: there is no operator* for Vector!
   // It would compile for Vector<Matrix<double> >, however, and return the sum of the products of the matrices


   std::vector<int> Indices;
   Indices.push_back(5);
   Indices.push_back(2);

   W[Indices] += 2.5 * U[Indices];

   TRACE(W);					// here, W = (10, 0, 35, 1, 10, -0.5, 10, 3, 10, 4)


   cout<<"\n\n############## max, min, abs, arg, real #####################################\n";
   Vector<complex<double> > X(size(U));
   real(X) = U;
   imag(X) = W;
   TRACE(U)(W)(X)(real(X))(imag(X))(abs(X))(arg(X));
   TRACE(max(U))(min(U));


   cout<<"\n\n############## NORMS ########################################################\n";
   TRACE(W)(VecOfVec);
   // calculate the Frobenius norm of W (= sqrt(sum of squares) )
   TRACE( norm_frob(W) )( norm_2(W) );
   TRACE( norm_frob(VecOfVec) )( norm_2(VecOfVec) );
   // calculate the square of the Frobenius norm of W
   TRACE( norm_frob_sq(W) );

   // calculate the inf-norm of W (= max(abs(W)) )
   TRACE( norm_inf(W) )(max(W));
   TRACE( norm_inf(VecOfVec) );



   /*// of course, we can combine operations in arbitrary ways
   cout << "scalar_prod(U.slice(EvenElements), W.slice(OddElements)) = "
	<< scalar_prod(U.slice(EvenElements), W.slice(OddElements)) << endl;
   // the previous line maps onto a single BLAS call
   */
}
