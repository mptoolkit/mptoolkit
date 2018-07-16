// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/gmres.h
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

#include "common/proccontrol.h"
#include <iostream>
#include <cmath>
#include "blas/functors.h"

using blas::norm_frob;
using blas::norm_frob_sq;

double const DGKS_Threshold = 1.0 / std::sqrt(2.0); // 1.0; // between 0 and 1.

template <typename Matrix, typename Vector1, typename Vector2>
void
Update(Vector1& x, int k, Matrix const& h, Vector2 const& s, Vector1 v[])
{
   Vector2 y(copy(s));

   // Backsolve:
   for (int i = k; i >= 0; i--) {
      y[i] /= h(i,i);
      for (int j = i - 1; j >= 0; j--)
	 y[j] -= h(j,i) * y[i];
   }

   for (int j = 0; j <= k; j++)
      x += v[j] * y[j];
}

template<typename Real>
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
   if (dy == 0.0)
   {
      cs = 1.0;
      sn = 0.0;
   }
   else if (norm_frob_sq(dy) > norm_frob_sq(dx))
   {
      Real temp = dx / dy;
      sn = 1.0 / std::sqrt( 1.0 + norm_frob_sq(temp) );
      cs = temp * sn;
   }
   else
   {
      Real temp = dy / dx;
      cs = 1.0 / std::sqrt( 1.0 + norm_frob_sq(temp) );
      sn = temp * cs;
   }

#if 0
   // debug: This is the matrix that we want to be Unitary
   blas::Matrix<std::complex<double> > M(2,2);
   M(0,0) = cs;
   M(1,1) = cs;
   M(0,1) = conj(sn) * (cs / conj(cs));
   M(1,0) = -sn;

   TRACE(M*herm(M));
#endif
}

template<typename Real>
void ApplyPlaneRotation(Real &dx, Real &dy, Real cs, Real sn)
{
   Real temp  = conj(cs) * dx + conj(sn) * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}

// GmRes algorithm
// Solve for x: MatVecMultiply(x) = b
// m = krylov subspace size (if this is reached, do a restart)
// using at most max_iter iterations
//
// On exit: tol = residual norm
// max_iter = number of iterations performed

template <typename Vector, typename MultiplyFunc, typename PrecFunc>
int
GmRes(Vector &x, MultiplyFunc MatVecMultiply, double normb, Vector const& b,
      int m, int& max_iter, double& tol, PrecFunc Precondition, int Verbose = 0)
{
  //  typedef typename Vector::value_type value_type;
  typedef std::complex<double> value_type;
  typedef blas::Vector<value_type> VecType;
  VecType s(m+1), cs(m+1), sn(m+1);
  blas::Matrix<value_type> H(m+1, m+1, 0.0);

  Vector w = Precondition(MatVecMultiply(x));
  Vector r = Precondition(copy(b)) - w; // - MatVecMultiply(x));
  double beta = norm_frob(r);

  // Initial guess for x:
  // let r = MatVecMultiply(x);
  // minimize choose prefactor to minimize norm_frob_sq(b - a*r)
  // This means
  // norm_frob_sq(b) + |a|^2 norm_frob_sq(r) - <b|ar> - <ar|b>
  // let <b|r> = c
  // then norm_frob_sq(b) + |a|^2 norm_frob_sq(r) - 2*real(ac)
  // So choose a = k*conj(c), for k real
  // Then minimize:
  // f(k) = norm_frob_sq(b) + k^2 |c|^2 norm_frob_sq(r) - 2*k*|c|^2
  // minimize with respect to k:
  // df/dk = 0 =>
  // 2k |c|^2 norm_frob_sq(r) - 2|c|^2 = 0
  // -> k = 1/norm_frob_sq(r)
  // so scale r => r*conj(<b|r>)/norm_frob(r)
  // r -> <r|b>/norm_frob(r)
  // In practice, this seems to have negligble effect versus simply scaling
  // the guess vector by the norm of b, so it is currently disabled.
  // Also I don't understand the complex conjugation here.

#if 0
  TRACE(norm_frob(r));
  r = MatVecMultiply(x);
  r = Precondition(b - (inner_prod(b,r)/norm_frob_sq(r))*r);
  beta = norm_frob(r);
  TRACE(norm_frob(r));
#endif

  DEBUG_TRACE(normb);
  if (normb == 0.0)
    normb = 1;
  double resid = norm_frob(r) / normb;
  DEBUG_TRACE(norm_frob(b))(norm_frob(MatVecMultiply(x)))(beta);
  DEBUG_TRACE(resid)(norm_frob(b - MatVecMultiply(x)) / norm_frob(b));

  DEBUG_TRACE(r)(norm_frob(w))(norm_frob(x));

  if ((resid = norm_frob(r) / normb) <= tol || norm_frob(w) / norm_frob(x) <= tol)
  {
#if 0
     if (norm_frob(w) / norm_frob(x) <= tol && Verbose > 0)
     {
        // This means that the matrix is effectively zero
        std::cerr << "GMRES: early exit at norm_frob(w)/norm_frob(b) < tol\n";
     }
     tol = resid;
     max_iter = 0;
     return 0;
#endif
  }

  Vector* v = new Vector[m+1];

  int j = 1;
  while (j <= max_iter)
  {
     v[0] = (1.0 / beta) * r;
     clear(s);
     s[0] = beta;

     int i = 0;
     while (i < m && j <= max_iter && resid >= tol)
     {
        if (Verbose > 2)
           std::cerr << "GMRES: iteration " << i << std::endl;

        double NormFrobSqH = 0; // for DGKS correction

        w = Precondition(MatVecMultiply(v[i]));
        for (int k = 0; k <= i; k++)
        {
           H(k, i) = inner_prod(v[k], w);
           w -= H(k, i) * v[k];
	   using blas::norm_frob_sq;
           NormFrobSqH += norm_frob_sq(H(k,i));
        }

        // Apply DGKS correction, if necessary
        double NormFrobSqF = norm_frob_sq(w);
        if (NormFrobSqF < DGKS_Threshold * DGKS_Threshold * NormFrobSqH)
        {
           DEBUG_TRACE("DGKS correction in GMRES")(NormFrobSqF / (DGKS_Threshold * DGKS_Threshold * NormFrobSqH));
           for (int k = 0; k <= i; k++)
           {
              value_type z = inner_prod(v[k], w);
              H(k, i) += z;
              w -= z * v[k];
           }
        }

        // Continue with our normal schedule...
        H(i+1, i) = norm_frob(w);
        v[i+1] = w * (1.0 / H(i+1, i));

        for (int k = 0; k < i; k++)
           ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);

        GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

        resid = norm_frob(s[i+1]) / normb;

        if (Verbose > 2)
           std::cerr << "GMRES: resid=" << resid << '\n';

        ++i;
        ++j;

#if 0
        // debugging check
        {
           Vector X2 = x;
           Update(X2, i-1, H, s, v);
           Vector R = Precondition(b - MatVecMultiply(X2));
           TRACE(i)(norm_frob(s[i]))(norm_frob(R));
	   X2 = R;
           for (int k = 0; k <= i; k++)
           {
	      X2 -= inner_prod(v[k], X2) * v[k];
           }
	   TRACE(norm_frob(X2));

	   
        }
#endif

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
        tol = UpdatedResid;
        max_iter = j;
        delete [] v;
        DEBUG_TRACE(beta)(normb);
        if (Verbose > 0)
           std::cerr << "GMRES: finished, iter=" << (j-1) << ", approx resid=" << resid
                     << ", actual resid=" << UpdatedResid << std::endl;
        return 0;
     }
     else
     {
        if (Verbose > 1)
           std::cerr << "GMRES: restarting, iter=" << (j-1) << ", resid=" << resid << '\n';
     }
     resid = beta / normb;
     DEBUG_TRACE(resid)(norm_frob(Precondition(b - MatVecMultiply(x))) / normb)
        (norm_frob(b - MatVecMultiply(x)) / norm_frob(b));
  }

  tol = resid;
  delete [] v;
  return 1;
}

template <typename Vector, typename MultiplyFunc, typename PrecFunc>
int
GmRes(Vector &x, MultiplyFunc MatVecMultiply, Vector const& b,
      int m, int& max_iter, double& tol, PrecFunc Precondition, int Verbose = 0)
{
   return GmRes(x, MatVecMultiply, norm_frob(b), b, m, max_iter, tol, Precondition, Verbose);
}

#endif
