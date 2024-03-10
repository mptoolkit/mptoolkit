// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/test/testcg.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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
#include "linearalgebra/vector.h"
#include "matrixproduct/conjugategradient.h"
#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/vector_utility.h"

typedef std::complex<double> complex;
typedef LinearAlgebra::Matrix<complex> Matrix;
typedef LinearAlgebra::Vector<complex> Vector;

struct DoMultiply
{
   DoMultiply(Matrix const& M) : M_(M) {}

   Vector operator()(Vector const& v) const
   {
      Vector Result(size1(M_), 0.0);
      for (std::size_t i = 0; i < size1(M_); ++i)
      {
         for (std::size_t j = 0; j < size2(M_); ++j)
         {
            Result[i] += M_(i,j) * v[j];
         }
      }
      return Result;
   }

   Matrix M_;
};

struct Precondition
{
   Precondition(Matrix const& M) : M_(M) {}

   Vector operator()(Vector const& v) const
   {
      Vector Res(v);
      for (std::size_t i = 0; i < size(Res); ++i)
      {
         Res[i] *= (1.0 / M_(i,i));
      }
      return Res;
   }
   Matrix M_;
};

int main()
{
   typedef std::complex<double> num;

   int const Sz = 100;

   Matrix M = LinearAlgebra::random_matrix<num>(Sz, Sz);

   // check CG for complex-hermitian
   Matrix P = M * herm(M) + LinearAlgebra::diagonal_matrix(Vector(Sz, 0.1));

   TRACE(EigenvaluesHermitian(P));

   Vector Rhs = LinearAlgebra::random_vector<num>(Sz);
   Vector v = Rhs;

   int MaxIter = 500;
   double Tol = 1E-14;

   ConjugateGradient(v, DoMultiply(P), Rhs, MaxIter, Tol, Precondition(P));
   TRACE(norm_2(Rhs - DoMultiply(P)(v)) / norm_2(Rhs));
   //   CHECK(LinearAlgebra::equal(DoMultiply(P)(v), Rhs, 1E-10))
   //      (DoMultiply(P)(v))(Rhs)(norm_2(Rhs - DoMultiply(P)(v)) / norm_2(Rhs));

   // check CG for complex symmetric
   P = M * trans(M)  + LinearAlgebra::diagonal_matrix(Vector(Sz, 1));

   MaxIter = 5000;
   Tol = 1E-14;
   ConjugateGradient(v, DoMultiply(P), Rhs, MaxIter, Tol,
                     LinearAlgebra::Identity<Vector>(),
                     LinearAlgebra::ParallelProd<Vector, Vector>(),
                     LinearAlgebra::ParallelProd<Vector, Vector>());

   TRACE(norm_2(Rhs - DoMultiply(P)(v)) / norm_2(Rhs))(MaxIter)(Tol);
   CHECK(LinearAlgebra::equal(DoMultiply(P)(v), Rhs, 1E-10))
      (DoMultiply(P)(v))(Rhs)(norm_2(Rhs - DoMultiply(P)(v)) / norm_2(Rhs));
}
