// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/test/test-arnoldi.cpp
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

#include "mp-algorithms/arnoldi.h"
#include "linearalgebra/matrix_utility.h"

using namespace LinearSolvers;
using namespace LinearAlgebra;

typedef std::complex<double> complex;

struct MultBy
{
   MultBy(Matrix<complex> const& M_) : M(M_) {}

   Vector<complex> operator()(Vector<complex> const& x) const
   {
      TRACE(size(x));
      Vector<complex> Result(size(x), 0.0);
      for (unsigned i = 0; i < size(x); ++i)
      {
         for (unsigned j = 0; j < size(x); ++j)
         {
            Result[i] += M(i,j) * x[j];
         }
      }
      return Result;
   }

   Matrix<complex> M;
};

void ShowVec(Matrix<complex> const& M, Vector<complex> xx)
{
   Vector<complex> Mxx = MultBy(M)(xx);
   for (int i = 0; i < size(xx); ++i)
   {
      TRACE(Mxx[i] / xx[i]);
   }
}

int main()
{
   Matrix<complex> M = LinearAlgebra::random_matrix<complex>(4,4);

   TRACE(M);

   Vector<complex> x(size1(M), 1.0);
   //   x[0] = 1.0;

   int Iter = 4;

   complex Theta = Arnoldi(x, MultBy(M), Iter);

   TRACE(Theta)(x)(Iter);

   TRACE((1.0/Theta)*MultBy(M)(x));

   Matrix<complex> L, R;
   Vector<complex> Eigen = Diagonalize(M, L, R);
   TRACE(Eigen);
   TRACE(R(0, all));

   ShowVec(M, R(0,all));
}
