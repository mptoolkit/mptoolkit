// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/simdiag.cpp
//
// Copyright (C) 2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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


// Adapted from
// https://uk.mathworks.com/matlabcentral/fileexchange/46794-simdiag-m

/*
Copyright (c) 2009, Christian B. Mendl
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of  nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "simdiag.h"

double OffDiagonalNorm2(LinearAlgebra::Matrix<std::complex<double>> const& M)
{
   double Result = 0;
   for (int i = 0; i < M.size1(); ++i)
      for (int j = 0; j < M.size2(); ++j)
         if (i != j)
            Result += LinearAlgebra::norm_frob_sq(M(i,j));

   return Result;
}

void TimesR(LinearAlgebra::Matrix<std::complex<double>>& A, int j, int k, std::complex<double> c, std::complex<double> s)
{
   // A = A*R, R = R(j,k,c,s)
   for (int i = 0; i < A.size1(); ++i)
   {
      auto Aij = A(i,j);
      auto Aik = A(i,k);
      A(i,j) = c*Aij + s*Aik;
      A(i,k) = -conj(s)*Aij + conj(c)*Aik;
   }
}

void Rotate(LinearAlgebra::Matrix<std::complex<double>>& A, int j, int k, std::complex<double> c, std::complex<double> s)
{
   // A = R'*A*R, R = R(j,k,c,s)
   TimesR(A, j, k, c, s);
   for (int i = 0; i < A.size2(); ++i)
   {
      auto Aji = A(j,i);
      auto Aki = A(k,i);
      A(j,i) = conj(c)*Aji + conj(s)*Aki;
      A(k,i) = -s*Aji + c*Aki;
   }
}

double Target(std::complex<double> c, std::complex<double> s, LinearAlgebra::Matrix<std::complex<double>> const& V)
{
   return norm_frob_sq(s*c*conj(V(LinearAlgebra::all,0)) -c*c*conj(V(LinearAlgebra::all,1)) + s*s*conj(V(LinearAlgebra::all,2)))
      + norm_frob_sq(s*c*V(LinearAlgebra::all,0) + s*s*V(LinearAlgebra::all,1) -c*c*V(LinearAlgebra::all,2));
}

std::tuple<std::complex<double>, std::complex<double>>
ExactMin(std::complex<double> a0, std::complex<double> a21, std::complex<double> a12)
{
   // Exact minimizer of
   // |s c conj(a0) - c^2 conj(a21) + s^2 conj(a12)|^2 + |s c a0 + s^2 a21 - c^2 a12|^2;
   // Refer to
   //	H. H. Goldstine and L. P. Horwitz, A Procedure for the
   //	Diagonalization of Normal Matrices, J. ACM (1959)

   double u = a0.real();
   double v = a0.imag();

   auto tmp = 0.5*(a21+conj(a12));
   auto r = abs(tmp);
   auto beta = arg(tmp);

   tmp = 0.5*(a21-conj(a12));
   auto s = abs(tmp);
   auto gamma = arg(tmp);

   auto nu = beta - gamma;
   auto sin_nu = std::sin(nu);
   auto cos_nu = std::cos(nu);

   auto L = u*v - 4*r*s*sin_nu;
   auto M = u*u - v*v + 4*(r*r - s*s);

   auto A = L*(r*r - s*s)*sin_nu + M*r*s;
   auto B = L*(r*r + s*s)*cos_nu;
   auto C = L*(r*r - s*s) + M*r*s*sin_nu;

   auto t = r*s*cos_nu*std::sqrt(M*M + 4*L*L);
   auto phi = 0.5*(std::atan2(-A*C + B*t, B*C + A*t) - beta - gamma);

   auto r_cos_ba = r*cos(beta+phi);
   auto s_sin_ga = s*sin(gamma+phi);
   auto kappa = u*u + v*v - 4*(r_cos_ba*r_cos_ba + s_sin_ga*s_sin_ga);
   auto lambda = 4*(u*r_cos_ba + v*s_sin_ga);
   auto theta = -0.25*atan2(-lambda, kappa);

   return std::make_tuple(std::cos(theta), std::sin(theta)*std::complex<double>(std::cos(phi), std::sin(phi)));
}

std::tuple<std::complex<double>, std::complex<double>>
ApproxMin(LinearAlgebra::Matrix<std::complex<double>> const& V)
{
   std::complex<double> c, s;
   std::tie(c,s) = ExactMin(V(0,0), V(0,1), V(0,2));
   double m = Target(c,s,V);
   for (int j = 1; j < V.size1(); ++j)
   {
      std::complex<double> c1, s1;
      std::tie(c1,s1) = ExactMin(V(j,0), V(j,1), V(j,2));
      double x = Target(c1,s1,V);
      if (x < m)
      {
         m = x;
         c = c1;
         s = s1;
      }
   }
   return std::make_tuple(c,s);
}

LinearAlgebra::Matrix<std::complex<double>> eye(int N)
{
   LinearAlgebra::Matrix<std::complex<double>> Result(N,N, 0.0);
   for (int i = 0; i < N; ++i)
      Result(i,i) = 1.0;
   return Result;
}

std::tuple<LinearAlgebra::Matrix<std::complex<double>>, std::vector<LinearAlgebra::Vector<double>>>
SimultaneousDiagonalizeHermitian(std::vector<LinearAlgebra::Matrix<std::complex<double>>> M)
{
   CHECK(M.size() > 0);
   int N = M[0].size1();
   CHECK_EQUAL(M[0].size1(), M[0].size2());

   // The python code has a preprocessing step of diagonalizing all matrices in turn.  This should be skippable.
   // Or maybe diagonalize some random linear combination.

   double Tol = 1E-15;

   int MaxIter = 100;
   int Iter = 0;

   double Off2 = 0;
   double NScale = 0;
   for (auto const& x : M)
   {
      NScale += norm_frob(x);
      Off2 += OffDiagonalNorm2(x);
   }

   LinearAlgebra::Matrix<std::complex<double>> Q = eye(N);

   while (Off2 > Tol*NScale)
   {
      for (int j = 0; j < N; ++j)
      {
         for (int k = j+1; k < N; ++k)
         {
            LinearAlgebra::Matrix<std::complex<double>> V(M.size(), 3);
            for (int m = 0; m < M.size(); ++m)
            {
               V(m,0) = M[m](j,j) - M[m](k,k);
               V(m,1) = M[m](j,k);
               V(m,2) = M[m](k,j);
            }
            std::complex<double> c, s;
            std::tie(c,s) = ApproxMin(V);
            TimesR(Q, j, k, c, s);
            for (int m = 0; m < M.size(); ++m)
            {
               Rotate(M[m], j, k, c, s);
            }
         }
      }
      Off2 = 0;
      NScale = 0;
      for (auto const& x : M)
      {
         NScale += norm_frob(x);
         Off2 += OffDiagonalNorm2(x);
      }
      if (++Iter > MaxIter)
      {
         std::cerr << "Maximum number of iterations exceeded.  Current relative error: " << (Off2 / NScale) << '\n';
         TRACE(Q);
         std::abort();
      }
   }
   // Eigenvalues
   std::vector<LinearAlgebra::Vector<double>> Eval(M.size());
   for (int j = 0; j < M.size(); ++j)
   {
      Eval[j] = LinearAlgebra::Vector<double>(N);
      for (int i = 0; i < N; ++i)
      {
         Eval[j][i] = M[j](i,i).real();
      }
   }
   return std::make_tuple(Q, Eval);
}
