// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ddmrg_functional.h
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
// Calculates the minimum value of the functional
// F = <x|A|x> + 2 Re <b|x>
// where A is some operator, and |b> is a supplied vector.
//
// If we choose A = (w+E-H)^2 + \eta^2
// and |b> = \eta * |lv>
// then this corresponds to dynamical DMRG, with |x> as the imaginary part of the correction vector.
//
// If we choose A = (w+E-H)^2 + \eta^2,
// |b> = -(w+E-H - i\eta) |lv>  (it is safe to project this onto the local Hilbert space!)
// then this corresponds to the residual norm shifted by a constant,
// ||r||^2 = F + <b|b>
// where |r> = |lv> - (w+E-H+i\eta)|x>
// and |x> is the correction vector.
//

#if !defined(DDMRG_FUNCTIONAL_H_FJH43Y34Y498HP)
#define DDMRG_FUNCTIONAL_H_FJH43Y34Y498HP

#include "linearalgebra/eigen.h"

double const OrthoThreshold = 1E-14;  // threshold for happy breakdown

template <typename Vector, typename MultiplyFunc,
          //          typename InnerProdFunc,
          typename PrecFunc>
double
FunctionalMinimize(Vector& x, MultiplyFunc MatVecMultiply, Vector const& b, int& m,
                   //                   InnerProdFunc InnerProd,
                   PrecFunc Precondition)
{
   using LinearAlgebra::norm_2_sq;

  typedef std::complex<double> value_type;
  typedef double real_type;

  LinearAlgebra::Matrix<value_type> H(m+1, m+1, 0.0);
  LinearAlgebra::Matrix<value_type> bR(m+1, 1);

  LinearAlgebra::Vector<Vector> Krylov(m+1);

  double normx = norm_frob(x);
  // Generate our Krylov subspace
  Krylov[0] = x;
  Krylov[0] *= 1.0 / normx;
  Vector w = MatVecMultiply(Krylov[0]);
  H(0,0) = inner_prod(Krylov[0], w);
  bR(0,0) = inner_prod(Krylov[0], b);
  w -= H(0,0) * Krylov[0];
  double wNorm = norm_frob(w);

  std::complex<double> f = (H(0,0) * normx + 2.0 * bR(0,0).real()) * normx;
  DEBUG_TRACE(f);

  int i = 1;
  for ( ; i <= m && wNorm >= OrthoThreshold; ++i)
  {
     w *= 1.0 / wNorm;
     Krylov[i] = w;
     w = MatVecMultiply(Krylov[i]);
     for (int Hack = 0; Hack < 2; ++Hack)
     {
        for (int k = 0; k < i; ++k)
        {
           value_type Overlap = inner_prod(Krylov[k], w);
           w -= Overlap * Krylov[k];
           H(k,i) += Overlap;
           H(i,k) = conj(H(k,i));
 }
     }
     H(i,i) = inner_prod(Krylov[i], w);
     w -= H(i,i) * Krylov[i];
     wNorm = norm_frob(w);
     //     bR(i,0) = InnerProd(Krylov[i], b);
     bR(i,0) = inner_prod(Krylov[i], b);

     {
        // Diagonalize our operator in the krylov subspace
        LinearAlgebra::Matrix<value_type> U = H(LinearAlgebra::range(0,i+1), LinearAlgebra::range(0,i+1));
        LinearAlgebra::Vector<double> Eigen = LinearAlgebra::DiagonalizeHermitian(U);

        // convert bR to our diagonal basis
        LinearAlgebra::Matrix<value_type> bRDiag = U * bR;

        // construct the output vector in the Krylov space
        double Out = 0.0;
        for (int j = 0; j <= i; ++j)
        {
           Out -= norm_2_sq(bRDiag(j,0)) / Eigen[j];
        }
        TRACE(Out);
     }
  }

  m = i-1;

  // Diagonalize our operator in the krylov subspace
  LinearAlgebra::Matrix<value_type> U = H;
  LinearAlgebra::Vector<double> Eigen = LinearAlgebra::DiagonalizeHermitian(U);

  // convert bR to our diagonal basis
  LinearAlgebra::Matrix<value_type> bRDiag = U * bR;
  DEBUG_TRACE(bRDiag);

  // construct the output vector in the Krylov space
  double Out = 0.0;
  LinearAlgebra::Matrix<value_type> Res(m+1, 1);
  for (int i = 0; i <= m; ++i)
  {
     Out -= norm_2_sq(bRDiag(i,0)) / Eigen[i];
     Res(i,0) = -bRDiag(i,0) / Eigen[i];
  }
  Res = herm(U) * Res;  // convert back to our Krylov basis

  // construct the real output vector
  x = Res(0,0) * Krylov[0];
  for (int i = 1; i <= m; ++i)
  {
     x += Res(i,0) * Krylov[i];
  }

  return Out;
}

// In this version, it is fully preconditioned.  We ignore H*|k_{n-1}> and instead
// use the Precondition( k_{n-1} ).  This is not normal preconditioning, which would
// instead act as Precondition( H*|k_{n-1}> ).
template <typename Vector, typename MultiplyFunc,
          //          typename InnerProdFunc,
          typename PrecFunc>
double
FunctionalMinimizeWithPre(Vector& x, MultiplyFunc MatVecMultiply, Vector const& b, int& m,
                          //                   InnerProdFunc InnerProd,
                          PrecFunc Precondition,
                          double Delta = -1, bool TwoStep = false, int MinIter = 0, double DebugShift = 0.0,
                          Vector const& DebugR = Vector())
{
   using LinearAlgebra::norm_2_sq;

  typedef std::complex<double> value_type;
  typedef double real_type;

  LinearAlgebra::Matrix<value_type> H(m+1, m+1, 0.0);
  LinearAlgebra::Matrix<value_type> bR(m+1, 1);

  LinearAlgebra::Vector<Vector> Krylov(m+1);

  // residual
  //double normR = norm_frob_sq(DebugR);
  //double r = norm_frob_sq(Precondition(x) - DebugR) / normR;
  //TRACE(r);

  double normx = norm_frob(x);
  // Generate our Krylov subspace
  Krylov[0] = x;
  Krylov[0] *= 1.0 / normx;
  bR(0,0) = inner_prod(Krylov[0], b);
  Vector Hw = MatVecMultiply(Krylov[0]);
  H(0,0) = inner_prod(Krylov[0], Hw);

  //  Vector w = Precondition(Krylov[0]);
  //  w -= inner_prod(Krylov[0], w);
  double wNorm = 1.0; // norm_frob(w);

  double LastF = (H(0,0).real() * normx + 2.0 * bR(0,0).real()) * normx;
  DEBUG_TRACE(LastF);

  // it seems that Delta is often small on the first iteration, force at least 2 iters.
  LastF = 0;

  bool H2Next = false;

  int i = 1;
  for ( ; i <= m && wNorm >= OrthoThreshold; ++i)
  {
     // Add another Krylov vector
     Vector w; // = Precondition(x) - DebugR;
#if 1
     if (TwoStep && H2Next)
     {
        H2Next = false;
        w = MatVecMultiply(Krylov[i-2]);
     }
     else
     {
        H2Next = true;
        w = Precondition(Krylov[i-1]);
     }
#endif

     // orthogonalize twice, for stability
     for (int Hack = 0; Hack < 2; ++Hack)
     {
        for (int k = 0; k < i; ++k)
        {
           w -= inner_prod(Krylov[k], w) * Krylov[k];
        }
     }
     wNorm = norm_frob(w);
     w *= 1.0 / wNorm;
     Krylov[i] = w;
     bR(i,0) = inner_prod(Krylov[i], b);

     // H matrix elements
     Hw = MatVecMultiply(Krylov[i]);
     for (int k = 0; k < i; ++k)
     {
        H(k,i) = inner_prod(Krylov[k], Hw);
        H(i,k) = conj(H(k,i));
     }
     H(i,i) = inner_prod(Krylov[i], Hw);


     // stopping test
     {
        // Diagonalize our operator in the krylov subspace
        LinearAlgebra::Matrix<value_type> U = H(LinearAlgebra::range(0,i+1), LinearAlgebra::range(0,i+1));
        LinearAlgebra::Vector<double> Eigen = LinearAlgebra::DiagonalizeHermitian(U);

        // convert bR to our diagonal basis
        LinearAlgebra::Matrix<value_type> bRDiag =
           U * bR(LinearAlgebra::range(0,i+1), LinearAlgebra::range(0,1));

        // construct the value of the functional
        double Out = 0.0;
        LinearAlgebra::Matrix<value_type> Res(i+1, 1);
        for (int j = 0; j <= i; ++j)
        {
           Out -= norm_2_sq(bRDiag(j,0)) / Eigen[j];
           Res(j,0) = -bRDiag(j,0) / Eigen[j];
        }
        Res = herm(U) * Res;  // convert back to our Krylov basis

#if 0
        // construct the real output vector
        x = Res(0,0) * Krylov[0];
        for (int j = 1; j <= i; ++j)
        {
           x += Res(j,0) * Krylov[j];
        }

        // residual
        double r = norm_frob_sq(Precondition(x) - DebugR) / normR;
        TRACE(r);
#endif

        // if we've converged, force an early return.  This is slightly dubious
        if (fabs(LastF - Out) < Delta && i >= MinIter)
           wNorm = 0.0;

        DEBUG_TRACE(Out)(fabs(LastF-Out));
        //        TRACE(Out+DebugShift);

        LastF = Out;
     }
  }

  m = i-1;

  // Diagonalize our operator in the krylov subspace
  LinearAlgebra::Matrix<value_type> U = H(LinearAlgebra::range(0,m+1), LinearAlgebra::range(0,m+1));
  LinearAlgebra::Vector<double> Eigen = LinearAlgebra::DiagonalizeHermitian(U);

  // convert bR to our diagonal basis
  LinearAlgebra::Matrix<value_type> bRDiag = U * bR(LinearAlgebra::range(0,m+1), LinearAlgebra::range(0,1));
  //  TRACE(bRDiag);

  // construct the output vector in the Krylov space
  double Out = 0.0;
  LinearAlgebra::Matrix<value_type> Res(m+1, 1);
  for (int i = 0; i <= m; ++i)
  {
     Out -= norm_2_sq(bRDiag(i,0)) / Eigen[i];
     Res(i,0) = -bRDiag(i,0) / Eigen[i];
  }
  Res = herm(U) * Res;  // convert back to our Krylov basis

  // construct the real output vector
  x = Res(0,0) * Krylov[0];
  for (int i = 1; i <= m; ++i)
  {
     x += Res(i,0) * Krylov[i];
  }

  DEBUG_TRACE(Out);
  //TRACE(Out+DebugShift);

  return Out;
}

// In this version, it is fully preconditioned.  We ignore H*|k_{n-1}> and instead
// use the Precondition( k_{n-1} ).  This is not normal preconditioning, which would
// instead act as Precondition( H*|k_{n-1}> ).
template <typename Vector, typename MultiplyFunc,
          //          typename InnerProdFunc,
          typename PrecFunc>
double
FunctionalMinimizeWithResidVector(Vector& x, MultiplyFunc MatVecMultiply, Vector const& b, int& m,
                          //                   InnerProdFunc InnerProd,
                                  PrecFunc Precondition, Vector const& RightSide, double DebugShift,
                                  double Delta = -1, bool TwoStep = false, int MinIter = 0)
{
   using LinearAlgebra::norm_2_sq;

  typedef std::complex<double> value_type;
  typedef double real_type;

  LinearAlgebra::Matrix<value_type> H(m+1, m+1, 0.0);
  LinearAlgebra::Matrix<value_type> bR(m+1, 1);

  LinearAlgebra::Vector<Vector> Krylov(m+1);

  // residual
  double normR = norm_frob_sq(RightSide);
  double r = norm_frob_sq(Precondition(x) - RightSide) / normR;
  //  TRACE(r);

  double normx = norm_frob(x);
  // Generate our Krylov subspace
  Krylov[0] = x;
  Krylov[0] *= 1.0 / normx;
  bR(0,0) = inner_prod(Krylov[0], b);
  Vector Hw = MatVecMultiply(Krylov[0]);
  H(0,0) = inner_prod(Krylov[0], Hw);

  //  Vector w = Precondition(Krylov[0]);
  //  w -= inner_prod(Krylov[0], w);
  double wNorm = 1.0; // norm_frob(w);

  double LastF = (H(0,0).real() * normx + 2.0 * bR(0,0).real()) * normx;
  DEBUG_TRACE(LastF);

  // it seems that Delta is often small on the first iteration, force at least 2 iters.
  LastF = 0;

  int i = 1;
  for ( ; i <= m && wNorm >= OrthoThreshold; ++i)
  {
     // Add another Krylov vector, this time using the best guess residual
     Vector w = Precondition(x) - RightSide;
     r = norm_frob_sq(w) / normR;
     //TRACE(r);

     // orthogonalize twice, for stability
     for (int Hack = 0; Hack < 2; ++Hack)
     {
        for (int k = 0; k < i; ++k)
        {
           w -= inner_prod(Krylov[k], w) * Krylov[k];
        }
     }
     wNorm = norm_frob(w);

     if (wNorm < OrthoThreshold)
        break;

     w *= 1.0 / wNorm;
     Krylov[i] = w;
     bR(i,0) = inner_prod(Krylov[i], b);

     // H matrix elements
     Hw = MatVecMultiply(Krylov[i]);
     for (int k = 0; k < i; ++k)
     {
        H(k,i) = inner_prod(Krylov[k], Hw);
        H(i,k) = conj(H(k,i));
     }
     H(i,i) = inner_prod(Krylov[i], Hw);


     // stopping test
     {
        // Diagonalize our operator in the krylov subspace
        LinearAlgebra::Matrix<value_type> U = H(LinearAlgebra::range(0,i+1), LinearAlgebra::range(0,i+1));
        LinearAlgebra::Vector<double> Eigen = LinearAlgebra::DiagonalizeHermitian(U);

        // convert bR to our diagonal basis
        LinearAlgebra::Matrix<value_type> bRDiag =
           U * bR(LinearAlgebra::range(0,i+1), LinearAlgebra::range(0,1));

        // construct the value of the functional
        double Out = 0.0;
        LinearAlgebra::Matrix<value_type> Res(i+1, 1);
        for (int j = 0; j <= i; ++j)
        {
           Out -= norm_2_sq(bRDiag(j,0)) / Eigen[j];
           Res(j,0) = -bRDiag(j,0) / Eigen[j];
        }
        Res = herm(U) * Res;  // convert back to our Krylov basis

        // construct the real output vector
        x = Res(0,0) * Krylov[0];
        for (int j = 1; j <= i; ++j)
        {
           x += Res(j,0) * Krylov[j];
        }

        // residual
        //        double r = norm_frob_sq(Precondition(x) - DebugR) / normR;
        //        TRACE(r);

        // if we've converged, force an early return.  This is slightly dubious
        if (fabs(LastF - Out) < Delta && i >= MinIter)
           wNorm = 0.0;

        DEBUG_TRACE(Out)(fabs(LastF-Out));
        //TRACE(Out+DebugShift);

        LastF = Out;
     }
  }

  m = i-1;

  // Diagonalize our operator in the krylov subspace
  LinearAlgebra::Matrix<value_type> U = H(LinearAlgebra::range(0,m+1), LinearAlgebra::range(0,m+1));
  LinearAlgebra::Vector<double> Eigen = LinearAlgebra::DiagonalizeHermitian(U);

  // convert bR to our diagonal basis
  LinearAlgebra::Matrix<value_type> bRDiag = U * bR(LinearAlgebra::range(0,m+1), LinearAlgebra::range(0,1));
  //  TRACE(bRDiag);

  // construct the output vector in the Krylov space
  double Out = 0.0;
  LinearAlgebra::Matrix<value_type> Res(m+1, 1);
  for (int i = 0; i <= m; ++i)
  {
     Out -= norm_2_sq(bRDiag(i,0)) / Eigen[i];
     Res(i,0) = -bRDiag(i,0) / Eigen[i];
  }
  Res = herm(U) * Res;  // convert back to our Krylov basis

  // construct the real output vector
  x = Res(0,0) * Krylov[0];
  for (int i = 1; i <= m; ++i)
  {
     x += Res(i,0) * Krylov[i];
  }

  DEBUG_TRACE(Out);
  //TRACE(Out+DebugShift);

  return Out;
}

#endif
