// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/gmres.h
//
// Copyright (C) 2004-2021 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
//
//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no onvergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

#if !defined(MPTOOLKIT_MP_ALGORITHMS_GMRES_H)
#define MPTOOLKIT_MP_ALGORITHMS_GMRES_H

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include "common/numerics.h"
#include <iostream>
#include <cmath>

using LinearAlgebra::range;
using LinearAlgebra::norm_2_sq;
using LinearAlgebra::norm_2;

double const DGKS_Threshold = 1.0 / std::sqrt(2.0); // 1.0; // between 0 and 1.

template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
void
Update(Vector1& x, int k, Matrix const& h, Vector2 const& s, Vector3 const& v)
{
   Vector2 y(s);

   // Backsolve:
   for (int i = k; i >= 0; i--)
   {
      y[i] /= h(i,i);
      for (int j = i - 1; j >= 0; j--)
         y[j] -= h(j,i) * y[i];
   }

   for (int j = 0; j <= k; j++)
      x += v[j] * y[j];
}

// Generate a plane rotation of a vector (dx, dy) that rotates it into (r, 0)
// Returns {cs, sn}, where the rotation matrix is
// R = ( cs         sn )
//     ( -conj(sn)  cs )
// following the scheme of arxiv:2211.04010
// This results in cs being real, so we can return a tuple of <real, complex>
template <typename Scalar>
auto GeneratePlaneRotation(Scalar dx, Scalar dy) {
   using RealType = decltype(std::real(dx));
   using std::conj;
   using std::sqrt;

   RealType const safmin = numerics::safmin<RealType>();

   RealType x2 = norm_2_sq(dx);
   RealType y2 = norm_2_sq(dy);
   RealType m2 = x2 + y2;

   // Check if m2/x2 could overflow
   if (m2*safmin <= x2)
   {
      RealType c = std::sqrt(x2/m2);
      auto r = dx / c;
      double const rtmin = sqrt(safmin);
      double const rtmax = 1.0 / rtmin;
      if (x2 > rtmin && m2 < rtmax)
      {
         auto s = conj(dy) * (dx / sqrt(x2 * m2));
         return std::make_tuple(c,s);
      }
      // else
      auto s = conj(dy) * (r / m2);
      return std::make_tuple(c,s);
   }
   // else m2/x2 might overflow
   RealType d = sqrt(x2*m2);
   RealType c = x2 / d;
   auto s = conj(dy) * (dx/d);
   return std::make_tuple(c,s);
}

// Apply the plane rotation to a vector (dx,dy)
// The rotation R is a tuple {cs, sn}, where the matrix is
// R = ( cs  sn          )
//     ( -conj(sn)   cs  )
template<typename Real, typename Scalar>
void ApplyPlaneRotation(std::tuple<Real, Scalar> const& R, Scalar &dx, Scalar &dy)
{
   using std::conj;
   Scalar temp = std::get<0>(R)*dx + std::get<1>(R)*dy;
   dy = -conj(std::get<1>(R))*dx + std::get<0>(R)*dy;
   dx = temp;
}

// GmRes algorithm
// Solve for x: MatVecMultiply(x) = b
// m = krylov subspace size (if this is reached, do a restart)
// using at most max_iter iterations, where we have done iter number of iterations so far
//
// On exit: tol = residual norm
// iter = number of iterations performed (added to the existing value)
// return value is 0

template <typename Vector, typename MultiplyFunc, typename Preconditioner>
int
GmRes(Vector &x, MultiplyFunc MatVecMultiply, double normb, Vector const& b,
      int m, int& iter, int max_iter, double& tol, Preconditioner Precondition, int Verbose = 0)
{
   //  typedef typename Vector::value_type value_type;
   typedef std::complex<double> value_type;
   using real_type = value_type::value_type;
   typedef LinearAlgebra::Vector<value_type> VecType;
   VecType s(m+1);
   std::vector<std::tuple<real_type, value_type>> R; // the array of plane rotations
   LinearAlgebra::Matrix<value_type> H(m+1, m+1, 0.0);

   Vector w = Precondition(MatVecMultiply(x));
   Vector r = Precondition(b) - w;
   double beta = norm_frob(r);

   // we want maximum number of iterations on this restart
   max_iter -= iter;

   // if the true norm of the right hand side is zero, then the residual calculation won't work.
   // The user needs to set eg normb=1 in that case.
   CHECK(normb > 0.0);

   double resid = norm_frob(r) / normb;

   int j = 1;
   while (j <= max_iter)
   {
      std::vector<Vector> v;  // the krylov vectors
      v.reserve(m+1);
      v.emplace_back((1.0 / beta) * r);
      zero_all(s);
      s[0] = beta;

      int i = 0;
      while (i < m && j <= max_iter && resid >= tol)
      {
         if (Verbose > 2)
            std::cerr << "GMRES: iteration " << (iter+i+1) << std::endl;

         double NormFrobSqH = 0; // for DGKS correction

         w = Precondition(MatVecMultiply(v[i]));
         for (int k = 0; k <= i; k++)
         {
            H(k, i) = inner_prod(v[k], w);
            w -= H(k, i) * v[k];
            NormFrobSqH += LinearAlgebra::norm_frob_sq(H(k,i));
         }

         // Apply DGKS correction, if necessary
         double NormFrobSqF = norm_frob_sq(w);
         if (NormFrobSqF < DGKS_Threshold * DGKS_Threshold * NormFrobSqH)
         {
            if (Verbose > 2)
               std::cerr << "GMRES: DGKS correction\n";
           //DEBUG_TRACE("DGKS correction in GMRES")(NormFrobSqF / (DGKS_Threshold * DGKS_Threshold * NormFrobSqH));
           for (int k = 0; k <= i; k++)
           {
              value_type z = inner_prod(v[k], w);
              H(k, i) += z;
              w -= z * v[k];
           }
         }

         // Continue with our normal schedule...
         H(i+1, i) = norm_frob(w);
         v.emplace_back(w * (1.0 / H(i+1, i)));

         for (int k = 0; k < i; k++)
           ApplyPlaneRotation(R[k], H(k,i), H(k+1,i));

         R.push_back(GeneratePlaneRotation(H(i,i), H(i+1,i)));
         ApplyPlaneRotation(R[i], H(i,i), H(i+1,i));
         ApplyPlaneRotation(R[i], s[i], s[i+1]);

         resid = norm_2(s[i+1]) / normb;

         if (Verbose > 2)

           std::cerr << "GMRES: resid=" << resid << '\n';

         ++i;
         ++j;
      }
      Update(x, i-1, H, s, v);
      r = Precondition(b - MatVecMultiply(x));
      beta = norm_frob(r);

      //DEBUG_TRACE(r)(beta);

      // use the old value of resid here, to avoid cases
      // where the recalculation no longer satisfies the convergence criteria
      if (resid < tol)
      {
         double UpdatedResid = beta / normb;
         double ExpectedResid = tol;
         tol = UpdatedResid;
         iter += j;
         //DEBUG_TRACE(beta)(normb);
         if (Verbose > 0)
           std::cerr << "GMRES: finished, iter=" << iter << ", approx resid=" << resid
                     << ", actual resid=" << UpdatedResid << std::endl;
         // If the exact residual is not close to the required residual, then indicate an error.
         // What we should do here is go back and do some more iterations.
         if (UpdatedResid > ExpectedResid*10)
            return 1;
         return 0;
      }
      else
      {
         if (Verbose > 1)
            std::cerr << "GMRES: restarting, iter=" << iter << ", resid=" << resid << '\n';
      }
      resid = beta / normb;
      //DEBUG_TRACE(resid)(norm_frob(Precondition(b - MatVecMultiply(x))) / normb)
      //   (norm_frob(b - MatVecMultiply(x)) / norm_frob(b));
   }

   tol = resid;
   iter += j;
   return 1;
}

// Old style GMRES, for compatability
template <typename Vector, typename MultiplyFunc, typename Preconditioner>
int
GmRes(Vector &x, MultiplyFunc MatVecMultiply, Vector const& b,
      int m, int& max_iter, double& tol, Preconditioner P, int Verbose = 0)
{
   int MaxIter = max_iter;
   int iter = 0;
   int Ret = GmRes(x, MatVecMultiply, norm_frob(P(b)), b, m, iter, MaxIter, tol, P, Verbose);
   max_iter = iter;
   return Ret;
}

template <typename Vector, typename MultiplyFunc, typename Preconditioner>
int
GmResRefineInitial(Vector &x, MultiplyFunc MatVecMultiply, Vector const& b,
      int m, int& max_iter, double& tol, Preconditioner P, int Verbose = 0)
{
   // GMRES with iterative refinement, starting from an initial good guess
   // the algorithm:
   // 1. let xRefine = 0, bRefine = A(xRefine) = 0
   // 2. solve for x: A(x) = b - bRefine
   // 3. set xRefine = xRefine + x
   // 4. set bRefine = a(xRefine)
   // 5. go back to step 2 until converged
   // The fixed point is x=0, at which point b = bRefine, and the solution is xRefine

   double OriginalTol = tol;
   int Iter = 0;
   double normb = norm_frob(b);

   Vector xRefine = x;
   x *= 0.0;
   int Ret = 1;
   while (Ret == 1 && Iter < max_iter)
   {
      Vector MyB = b - MatVecMultiply(xRefine);
      tol = OriginalTol;
      Ret = GmRes(x, MatVecMultiply, normb, MyB, m, Iter, Iter+m, tol, P, Verbose);
      xRefine = xRefine + x;
      x *= 0.0;
   }

   x = xRefine;
   max_iter = Iter;
   return Ret;
}

template <typename Vector, typename MultiplyFunc, typename Preconditioner>
int
GmResRefine(Vector &x, MultiplyFunc MatVecMultiply, Vector const& b,
      int m, int& max_iter, double& tol, Preconditioner P, int Verbose = 0)
{
   // Assume the initial vector isn't good, so do one round of GMRES before iterative refinement
   // We can force this simply by setting the initial vector x to zero
   double OriginalTol = tol;
   int Iter = 0;
   double normb = norm_frob(b);

   x *= 0.0;
   Vector xRefine = x;
   int Ret = 1;
   while (Ret == 1 && Iter < max_iter)
   {
      Vector MyB = b - MatVecMultiply(xRefine);
      tol = OriginalTol;
      Ret = GmRes(x, MatVecMultiply, normb, MyB, m, Iter, Iter+m, tol, P, Verbose);
      xRefine = xRefine + x;
      x *= 0.0;
   }

   x = xRefine;
   max_iter = Iter;
   return Ret;
}


template <typename Vector, typename MultiplyFunc, typename Preconditioner>
int
GmResRefineOrtho(Vector &x, Vector const& OrthoLeft, Vector const& OrthoRight, MultiplyFunc MatVecMultiply, Vector const& b,
      int m, int& max_iter, double& tol, Preconditioner P, int Verbose = 0)
{
   // Assume the initial vector isn't good, so do one round of GMRES before iterative refinement
   // We can force this simply by setting the initial vector x to zero
   double OriginalTol = tol;
   int Iter = 0;
   double normb = norm_frob(b);

   x *= 0.0;
   Vector xRefine = x;
   int Ret = 1;
   while (Ret == 1 && Iter < max_iter)
   {
      Vector MyB = b - MatVecMultiply(xRefine);
      tol = OriginalTol;
      Ret = GmRes(x, MatVecMultiply, normb, MyB, m, Iter, Iter+m, tol, P, Verbose);
      x -= std::conj(inner_prod(x, OrthoLeft)) * OrthoRight;
      xRefine = xRefine + x;
      x *= 0.0;
   }

   x = xRefine;
   max_iter = Iter;
   return Ret;
}

#endif
