

// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectortransform.h
//
// Copyright (C) 2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2008-2009, Christian Mendl
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

#include "takagi.h"

// TakagiFactor.c
// computes the Takagi factorization of a complex symmetric matrix
// code adapted from the "Handbook" routines
// (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)

std::tuple<LinearAlgebra::Matrix<std::complex<double>>, LinearAlgebra::DiagonalMatrix<double>>
TakagiFactor(LinearAlgebra::Matrix<std::complex<double>> A)
{
   CHECK_EQUAL(A.size1(), A.size2());
   int n = A.size1();

   /* A matrix is considered diagonal if the sum of the squares
      of the off-diagonal elements is less than EPS.  SYM_EPS is
      half of EPS since only the upper triangle is counted for
      symmetric matrices.
      52 bits is the mantissa length for IEEE double precision. */

   double constexpr eps = 1e-15;
   double constexpr sym_eps = 0.5e-15;
   double constexpr dbl_eps = 1e-15;

   double constexpr MaxIter = 50;

   double Red = 0.4 / std::pow(n, 4);
   LinearAlgebra::Vector<std::complex<double>> ev1(n, 0.0), ev0(n, 0.0);

   LinearAlgebra::Matrix<std::complex<double>> U(n, n, 0.0);
   LinearAlgebra::Vector<double> d(n, 0.0); // the eigenvalues

   for (int i = 0; i < n; ++i)
   {
      ev1[i] = A(i,i);
   }

   int Iter = 0;
   while (Iter < MaxIter)
   {
      double thresh = 0;
      for (int q = 1; q < n; ++q)
      {
         for (int p = 0; p < q; ++p)
         {
            thresh += LinearAlgebra::norm_frob_sq(A(p,q));
         }
      }
      if(thresh <= eps)
         break;

      thresh = (Iter <= 4) ? thresh*Red : 0;

      for(int q = 1; q < n; ++q)
      {
         for(int p = 0; p < q; ++p )
         {
            std::complex<double> Apq = A(p,q);
            double off = LinearAlgebra::norm_frob_sq(Apq);
            double sqp = LinearAlgebra::norm_frob_sq(ev1[p]);
            double sqq = LinearAlgebra::norm_frob_sq(ev1[q]);

            if (Iter > 5 && off < sym_eps*(sqp + sqq))
               A(p,q) = 0;
            else if (off > thresh)
            {
               double t = 0.5*std::abs(sqp - sqq);
               std::complex<double> f = (t > dbl_eps) ?
                  std::copysign(1, sqp - sqq)*(ev1[q]*conj(Apq) + conj(ev1[p]*Apq))
               :
                  (sqp == 0) ? 1 : std::sqrt(ev1[q]/ev1[p]);

               t += std::sqrt(t*t + LinearAlgebra::norm_frob_sq(f));
               f /= t;

               ev0[p] += Apq*conj(f);
               ev1[p] = A(p,p) + ev0[p];
               ev0[q] -= Apq*f;
               ev1[q] = A(q,q) + ev0[q];

               t = LinearAlgebra::norm_frob_sq(f);
               double invc = std::sqrt(t + 1.0);
               f /= invc;
               t /= invc*(invc + 1);

               for (int j = 0; j < p; ++j)
               {
                  auto x = A(j,p);
                  auto y = A(j,q);
                  A(j,p) = x + (conj(f)*y - t*x);
                  A(j,q) = y - (f*x + t*y);
               }

               for(int j = p + 1; j < q; ++j)
               {
                  auto x = A(p,j);
                  auto y = A(j,q);
                  A(p,j) = x + (conj(f)*y - t*x);
                  A(j,q) = y - (f*x + t*y);
               }

               for(int j = q + 1; j < n; ++j)
               {
                  auto x = A(p,j);
                  auto y = A(q,j);
                  A(p,j) = x + (conj(f)*y - t*x);
                  A(q,j) = y - (f*x + t*y);
               }

               A(p,q) = 0;

               for(int j = 0; j < n; ++j )
               {
                  auto x = U(p,j);
                  auto y = U(q,j);
                  U(p,j) = x + (f*y - t*x);
                  U(q,j) = y - (conj(f)*x + t*y);
               }
            }
         }
      }
      for (int p = 0; p < n; ++p)
      {
         ev0[p] = 0;
         A(p,p) = ev1[p];
      }
      ++Iter;
   }

   if (Iter >= MaxIter)
   {
      std::cerr << "Failed convergence in TakagiFactor\n";
      std::abort();
   }

   // make the diagonal elements nonnegative

   for (int p = 0; p < n; ++p)
   {
      auto App = A(p,p);
      d[p] = LinearAlgebra::norm_frob_sq(App);
      if(d[p] > dbl_eps && d[p] != App.real())
      {
         auto f = std::sqrt(App/d[p]);
         for(int q = 0; q < n; ++q)
            U(p,q) *= f;
      }
   }

   // we could sort the eigenvalues if we wanted
   LinearAlgebra::DiagonalMatrix<double> D(n,n);
   D.diagonal() = d;
   return std::make_tuple(U, D);
}
