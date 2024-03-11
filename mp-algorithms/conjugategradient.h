// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/conjugategradient.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
//
//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// CG follows the algorithm described on p. 15 in the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

using LinearAlgebra::range;

template <typename Vector, typename MultiplyFunctor,
          typename PreFunctor, typename InnerProdLRFunctor,
          typename InnerProdRLFunctor>
void
ConjugateGradient(Vector &x, MultiplyFunctor MatVecMultiply, Vector const& b,
                  int& max_iter, double& tol,
                  PreFunctor Precondition,
                  InnerProdLRFunctor InnerProd_LR,
                  InnerProdRLFunctor InnerProd_RL)
{
  double resid;
  Vector p, z, q;
  typedef std::complex<double> value_type;
  value_type alpha, beta, rho, rho_1;

  Vector Tx; Vector Px;

  double normb = norm_frob(b);
  Vector r = b - MatVecMultiply(x);

#if defined(CHECK_KRYLOV)
  std::vector<Vector> Krylov;
  Krylov.push_back(r);
#endif

  if (normb == 0.0)
     normb = 1;

  if ((resid = norm_frob(r) / normb) <= tol)
  {
     tol = resid;
     max_iter = 0;
     return;
  }

  TRACE(resid);

  for (int i = 1; i <= max_iter; i++)
  {
     z = Precondition(r);
     rho = InnerProd_RL(r, z);

     //TRACE(rho);

     //    rho = conj(rho);

     //     TRACE(rho)(norm_frob(z))(norm_frob(r))(InnerProd(z, conj(r)));

     if (i == 1)
        p = z;
     else
     {
        beta = rho / rho_1;
        //        TRACE(beta);

        //      Px = q;
        p = z + beta * p;

        //      TRACE(inner_prod(p,Px));
     }

     q = MatVecMultiply(p);
     alpha = rho / InnerProd_LR(p, q);
     TRACE(norm_frob(r))(InnerProd_RL(r,r))(rho)(alpha)(norm_frob(p))(norm_frob(q));

     //     TRACE(norm_frob_sq(p))(norm_frob_sq(q))(InnerProd(p,q))(InnerProd(p,conj(q)));

     //     TRACE(norm_frob(r - alpha*q))(norm_frob(r + alpha*q));
     //     TRACE(norm_frob(r + conj(alpha)*q))(norm_frob(r - conj(alpha)*q));
     //     TRACE(norm_frob(r + conj(alpha*q)))(norm_frob(r - conj(alpha*q)));
     //     TRACE(norm_frob(r + alpha*conj(q)))(norm_frob(r - alpha*conj(q)));

     //     Tx = r;

     x += alpha * p;
     r -= alpha * q;

     //TRACE(inner_prod(Tx, r));

#if defined(CHECK_KRYLOV)
     for (std::size_t n = 0; n < Krylov.size(); ++n)
     {
        TRACE(i)(n)(inner_prod(r, Krylov[n]));
     }
     Krylov.push_back(r);
#endif

     if ((resid = norm_frob(r) / normb) <= tol)
     {
        tol = resid;
        max_iter = i;
      return;
     }

     TRACE(resid);
     rho_1 = rho;
  }

  tol = resid;
  // at this point, we did not converge
  return;
}

template <typename Vector, typename MultiplyFunctor, typename Precondition>
inline
void
ConjugateGradient(Vector &x, MultiplyFunctor const& MatVecMultiply, Vector const& b,
                  int& max_iter, double& tol, Precondition P)
{
   return ConjugateGradient(x, MatVecMultiply, b, max_iter, tol,
                            P,
                            LinearAlgebra::InnerProd<Vector, Vector>(),
                            LinearAlgebra::InnerProd<Vector, Vector>());
}

template <typename Vector, typename MultiplyFunctor>
inline
void
ConjugateGradient(Vector &x, MultiplyFunctor const& MatVecMultiply, Vector const& b,
                  int& max_iter, double& tol)
{
   return ConjugateGradient(x, MatVecMultiply, b, max_iter, tol,
                            LinearAlgebra::Identity<Vector>(),
                            LinearAlgebra::InnerProd<Vector, Vector>(),
                            LinearAlgebra::InnerProd<Vector, Vector>());
}
