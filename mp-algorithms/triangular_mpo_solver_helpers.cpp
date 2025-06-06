// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/triangular_mpo_solver_helpers.cpp
//
// Copyright (C) 2009-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "triangular_mpo_solver_helpers.h"

// maximum number of iterations used in GMRES.  In some cases of slow convergence
// it is necessary to increase this.
int const MaxIter = getenv_or_default("MP_GMRES_MAXITER", 50000);

// A note on orthogonalizing vectors in non-orthogonal Hilbert spaces.
// Suppose <l| and |r> are left and right eigenvectors, satisfying <l|r>=1.
// If we want to orthogonalize some bra <x| against <l|, we want do do:
// <x'| = <x| - conj(c)*<l|
// where c = <x|r>
// Proof: <x'|r> = <x|r> - c*<l|r> = 0
//
// In the other direction, if we want to orthogonalize |y> against |r>, then
// |y'> = |y> - c*|r>
// where c = <l|y>

// binomial function.  This can be implemented with integer arithmetic.
// gmp contains a version with arbitrary size integers, but we're never going
// to overflow a long (or even an int!) with an MPO contraction.
long Binomial(int n, int k)
{
   if (k > n/2)
      k = n-k;     // take advantage of symmetry
   long r = 1;
   long n_k = n-k;
   for (long i = 1; i <= k; ++i)
   {
      DEBUG_CHECK((r*(n_k-i))%i == 0);
      r = (r * (n_k+i)) / i;
   }
   return r;
}

// Given a matrix polynomial C, extract the components proportional to the
// Rho into a c-number polynomial
ComplexPolyType
FindParallelParts(MatrixPolyType& C,
                  MatrixOperator const& Identity,
                  MatrixOperator const& Rho, int MaxDegree)
{
   ComplexPolyType CParallel;
   // iterate over components of the polynomial
   for (MatrixPolyType::iterator I = C.begin(); I != C.end(); ++I)
   {
      // Orthogonalize against the identity
      std::complex<double> Overlap = inner_prod(I->second, Rho);
      DEBUG_TRACE(Overlap);
      I->second -= std::conj(Overlap)*Identity;

      // **comparison** we always want to add here to get the degree of the polynomial correct.
      // This is the important one
      //      if (norm_frob(Overlap) > 1E-16)

      // Only keep parts if they are degree one less (or more) than the maximum degree.
      if ((MaxDegree == 0 || I->first <= MaxDegree))
      {
         CParallel[I->first] = Overlap;
         //DEBUG_TRACE("Adding component")(I->first)(Overlap);
      }
      else if (norm_frob(Overlap) > 1E-10)
      {
         std::cerr << "Warning: ignoring component at degree " << I->first << " magnitude " << norm_frob(Overlap) << '\n';
      }
   }
   return CParallel;
}

// Update the components of the E polynomial, given CParallel
ComplexPolyType
UpdateParallelParts(ComplexPolyType const& CParallel, std::complex<double> MomentumFactor)
{
   ComplexPolyType EParallel;
   for (int m = CParallel.degree(); m >= 0; --m)
   {
      if (CParallel.has_term(m))
      {
         EParallel[m+1] = std::conj(MomentumFactor) * CParallel[m]; // include momentum
      }

      for (int k = m+2; k <= CParallel.degree()+1; ++k)
      {
         if (EParallel.has_term(k))
         {
            EParallel[m+1] -= double(Binomial(k,m)) * EParallel[k];
         }
      }
      if (EParallel.has_term(m+1))
      {
         EParallel[m+1] *= 1.0 / (1.0 + m);
         //DEBUG_TRACE("Component at ")(m+1)(EParallel[m+1]);
      }
   }
   return EParallel;
}

// Extract the parallel parts of C, and return the corresponding E
ComplexPolyType
DecomposeParallelParts(MatrixPolyType& C,
                       MatrixOperator const& Identity,
                       MatrixOperator const& Rho, double UnityEpsilon, int Degree)
{
   ComplexPolyType CParallel = FindParallelParts(C, Identity, Rho, Degree);
   return UpdateParallelParts(CParallel);
}

// Extract the parallel parts of C, and return the corresponding E, with momentum.
// The Identity and Rho matrices can be any left/right eigenpair.
KComplexPolyType
DecomposeParallelPartsWithMomentum(KMatrixPolyType& C,
                                   std::complex<double> Factor,
                                   MatrixOperator const& Identity,
                                   MatrixOperator const& Rho, double UnityEpsilon, int Degree)
{
   KComplexPolyType EParallel;
   // diagonal element is the identity, up to a unitary factor
   DEBUG_TRACE("Unit diagonal element")(Factor);

   // decompose C into components parallel and perpendicular to the identity
   // The only part we have to care about is a component with the same momentum as our unit operator
   for (KMatrixPolyType::iterator Ck = C.begin(); Ck != C.end(); ++Ck) // sum over momenta
   {
      std::complex<double> K = Ck->first;
      ComplexPolyType CParallel = FindParallelParts(Ck->second, Identity, Rho, Degree);

      // Is this the same momentum as our unit operator?
      if (norm_frob(K - Factor) < UnityEpsilon*10)
      {
         DEBUG_TRACE("Component at equal momenta")(K);
         // same momenta, these components diverge
         DEBUG_TRACE(CParallel);
         for (int m = CParallel.degree(); m >= 0; --m)
         {
            EParallel[Factor] = UpdateParallelParts(CParallel, Factor);
         }
      }
      else // different momentum
      {
         for (int m = CParallel.degree(); m >= 0; --m)
         {
            // different momenta, we get a contribution both at K and Factor
            if (CParallel.has_term(m))
            {
               std::complex<double> Term = CParallel[m] / (K - Factor);
               EParallel[K][m] += Term;
               EParallel[Factor][m] -= Term;
            }
         }
      }
   }
   return EParallel;
}

MatrixPolyType
DecomposePerpendicularPartsLeft(MatrixPolyType const& C, std::complex<double> K,
                           BasicFiniteMPO const& Diag,
                           MatrixOperator const& Identity,
                           MatrixOperator const& Rho,
                           LinearWavefunction const& Psi1,
                           LinearWavefunction const& Psi2,
                           QuantumNumber const& QShift,
                           double TCond,
                           bool HasEigenvalue1,
                           double Tol,
                           int Verbose)
{
   // Identity and Rho are only used if HasEigenvalue1 is true
   // Components perpendicular to the identity satisfy equation (24)
   MatrixPolyType E;

   for (int m = C.degree(); m >= 0; --m)
   {
      if (Verbose > 0)
         std::cerr << "Degree " << m << std::endl;

      //DEBUG_TRACE("degree")(m);
      MatrixOperator Rhs = std::conj(K) * C[m];

      for (int k = m+1; k <= E.degree(); ++k)
      {
         // avoid accessing E[k] if it doesn't exist, to avoid adding a null term
         if (E.has_term(k))
            Rhs -= double(Binomial(k,m)) * E[k];
      }

      // orthogonalize Rhs against the identity again, which is a kind of
      // DGKS correction
      if (HasEigenvalue1 && Rhs.TransformsAs() == Rho.TransformsAs())
      {
         //DEBUG_TRACE(inner_prod(Rhs, Rho))("should be small");
         Rhs -= std::conj(inner_prod(Rhs, Rho)) * Identity;
         Rhs -= std::conj(inner_prod(Rhs, Rho)) * Identity;
         //DEBUG_TRACE(inner_prod(Rhs, Rho))("should be zero");
         //DEBUG_TRACE(inner_prod(Rhs, Identity));
      }

      double RhsNorm2 = norm_frob_sq(Rhs);
      RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
      //DEBUG_TRACE(RhsNorm2);

      if (Verbose > 2)
      {
         std::cerr << "Orthogonality of perpendicular part, RhsNorm2=" << RhsNorm2 << '\n';
      }

      if (RhsNorm2 > 0.0)
      {
         // Initial guess vector -- scale it by the norm of Rhs, improves the stability
         E[m] = std::sqrt(RhsNorm2) *
            MakeRandomMatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
         if (HasEigenvalue1 && Rhs.TransformsAs() == Rho.TransformsAs())
         {
            // Orthogonalize the initial guess -- this is important for the numerical stability
            E[m] -= std::conj(inner_prod(E[m], Rho)) * Identity;
            //DEBUG_TRACE("should be zero")(inner_prod(E[m], Rho));

            LinearSolveOrtho(E[m], Rho, Identity, OneMinusTransferLeft_Ortho(Psi1, QShift, K*Diag, Psi2,
                             Identity, Rho, HasEigenvalue1),
                             TCond, Rhs, Tol, Verbose);

            // do another orthogonalization -- this should be unncessary but for the paranoid...
            std::complex<double> z = inner_prod(E[m], Rho);
            E[m] -= std::conj(z) * Identity;
            z = inner_prod(E[m], Rho);
            if (Verbose > 2)
            {
               std::cerr << "z=" << norm_frob_sq(z) << '\n';
            }
            if (LinearAlgebra::norm_frob_sq(z) > Tol * RhsNorm2)
            {
               WARNING("Possible numerical instability in triangular MPO solver")(z)(RhsNorm2);
            };
            E[m] -= std::conj(z) * Identity;
            //DEBUG_TRACE(inner_prod(E[m], Rho))("should be zero");
         }
         else
         {
            LinearSolve(E[m], OneMinusTransferLeft_Ortho(Psi1, QShift, K*Diag, Psi2,
                        Identity, Rho, HasEigenvalue1),
                        Rhs, Tol, Verbose);
         }
      }
      else
      {
         E[m] = MatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
      }
   }
   return E;
}

KMatrixPolyType
DecomposePerpendicularPartsLeft(KMatrixPolyType const& C,
                            BasicFiniteMPO const& Diag,
                            MatrixOperator const& Identity,
                            MatrixOperator const& Rho,
                            LinearWavefunction const& Psi1,
                            LinearWavefunction const& Psi2,
                            QuantumNumber const& QShift,
                            double TCond,
                            bool HasEigenvalue1,
                            double Tol,
                            int Verbose)
{
   // Identity and Rho are only used if HasEigenvalue1 is true
   // Components perpendicular to the identity satisfy equation (24)
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      std::complex<double> K = I->first;  // the momentum (complex phase)
      E[K] = DecomposePerpendicularPartsLeft(I->second, I->first, Diag, Identity, Rho,
                  Psi1, Psi2, QShift, TCond, HasEigenvalue1, Tol, Verbose);
   }
   return E;
}

// On entry, Rho ~ Identity, InitMatrixLeft ~ Rho
MatrixPolyType
DecomposePerpendicularPartsRight(MatrixPolyType const& C, std::complex<double> K,
                                 BasicFiniteMPO const& Diag,
                                 MatrixOperator const& Identity,
                                 MatrixOperator const& Rho,
                                 LinearWavefunction const& Psi1,
                                 LinearWavefunction const& Psi2,
                                 QuantumNumber const& QShift,
                                 double TCond,
                                 bool HasEigenvalue1,
                                 double Tol,
                                 int Verbose)
{
   // Identity and Rho are only used if HasEigenvalue1 is true
   // Components perpendicular to the identity satisfy equation (24)
   MatrixPolyType F;

   for (int m = C.degree(); m >= 0; --m)
   {
      if (Verbose > 0)
         std::cerr << "Degree " << m << std::endl;

      //DEBUG_TRACE("degree")(m);
      MatrixOperator Rhs = std::conj(K) * C[m];
      for (int k = m+1; k <= F.degree(); ++k)
      {
         // avoid accessing E[k] if it doesn't exist, to avoid adding a null term
         if (F.has_term(k))
            Rhs -= double(Binomial(k,m)) * F[k];
      }

      // orthogonalize Rhs against the identity again, which is a kind of
      // DGKS correction
      if (HasEigenvalue1 && Rhs.TransformsAs() == Identity.TransformsAs())
      {
         Rhs -= std::conj(inner_prod(Rhs, Rho)) * Identity;
         Rhs -= std::conj(inner_prod(Rhs, Rho)) * Identity;
      }

      double RhsNorm2 = norm_frob_sq(Rhs);
      RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
      //DEBUG_TRACE(RhsNorm2);

      if (RhsNorm2 > 0.0)
      {
         // Initial guess vector -- scale it by the norm of Rhs, improves the stability
         F[m] = std::sqrt(RhsNorm2) *
            MakeRandomMatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
         if (HasEigenvalue1 && Rhs.TransformsAs() == Rho.TransformsAs())
         {
            // Orthogonalize the initial guess -- this is important for the numerical stability
            F[m] -= std::conj(inner_prod(F[m], Rho)) * Identity;

            LinearSolveOrtho(F[m], Rho, Identity, OneMinusTransferRight_Ortho(Psi1, QShift, K*Diag, Psi2,
                             Identity, Rho, HasEigenvalue1),
                             TCond, Rhs, Tol, Verbose);

            //DEBUG_TRACE(m)(norm_frob(F[m]))(inner_prod(F[m], Rho));

            // do another orthogonalization -- this should be unncessary but for the paranoid...
            std::complex<double> z = inner_prod(F[m], Rho);
            F[m] -= std::conj(z) * Identity;
            z = inner_prod(F[m], Rho);
            if (Verbose > 2)
            {
               std::cerr << "z=" << norm_frob_sq(z) << '\n';
            }
            if (LinearAlgebra::norm_frob_sq(z) > Tol * RhsNorm2)
            {
               WARNING("Possible numerical instability in triangular MPO solver")(z);
            };
            F[m] -= std::conj(z) * Identity;
         }
         else
         {
            LinearSolve(F[m], OneMinusTransferRight_Ortho(Psi1, QShift, K*Diag, Psi2,
                        Identity, Rho, HasEigenvalue1),
                        Rhs, Tol, Verbose);
         }
      }
      else
      {
         F[m] = MatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
      }
   }
   return F;
}

KMatrixPolyType
DecomposePerpendicularPartsRight(KMatrixPolyType const& C,
                            BasicFiniteMPO const& Diag,
                            MatrixOperator const& Identity,
                            MatrixOperator const& Rho,
                            LinearWavefunction const& Psi1,
                            LinearWavefunction const& Psi2,
                            QuantumNumber const& QShift,
                            double TCond,
                            bool HasEigenvalue1,
                            double Tol,
                            int Verbose)
{
   // Identity and Rho are only used if HasEigenvalue1 is true
   // Components perpendicular to the identity satisfy equation (24)
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      std::complex<double> K = I->first;  // the momentum (complex phase)
      E[K] = DecomposePerpendicularPartsRight(I->second, I->first, Diag, Identity, Rho,
                  Psi1, Psi2, QShift, TCond, HasEigenvalue1, Tol, Verbose);
   }
   return E;
}

// Solve the components for the case where the diagonal operator is zero
MatrixPolyType
SolveZeroDiagonal(MatrixPolyType const& C, std::complex<double> K)
{
   // See equation (20)
   // E(n) = C(n-1)
   // which expands to, for the degree d component of the polynomial with momentum k,
   // E^k_d = e^{-ik} C^k_d - \sum_{j=d+1}^p (j d) E^k_j
   MatrixPolyType E;
   int MaxDegree = C.degree();
   for (int i = MaxDegree; i >= 0; --i)
   {
      if (C.has_term(i))
      {
         E[i] = std::conj(K) * C[i];
      }
      for (int j = i+1; j <= MaxDegree; ++j)
      {
         DEBUG_CHECK(!E[j].is_null());
         E[i] -= double(Binomial(j,i)) * E[j];
      }
   }
   return E;
}

KMatrixPolyType
SolveZeroDiagonal(KMatrixPolyType const& C)
{
   // See equation (20)
   // E(n) = C(n-1)
   // which expands to, for the degree d component of the polynomial with momentum k,
   // E^k_d = e^{-ik} C^k_d - \sum_{j=d+1}^p (j d) E^k_j
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      E[I->first] = SolveZeroDiagonal(I->second, I->first);
   }
   return E;
}
