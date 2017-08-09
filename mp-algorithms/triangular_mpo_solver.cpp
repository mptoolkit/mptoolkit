// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/triangular_mpo_solver.cpp
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "triangular_mpo_solver.h"
#include "wavefunction/operator_actions.h"
#include "mp-algorithms/arnoldi.h"
#include "mp-algorithms/gmres.h"
#include "common/environment.h"

// maximum number of iterations used in GMRES.  In some cases of slow convergence
// it is necessary to increase this.
int const MaxIter = getenv_or_default("MP_GMRES_MAXITER", 10000);

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

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator const& Rhs, double Tol = 1E-14, int Verbose = 0)
{
   MatrixOperator Guess = MakeRandomMatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
   LinearSolve(Guess, F, Rhs, Tol, Verbose);
   return Guess;
}

//#define USE_ITERATIVE_REFINEMENT
// iterative refinement doesn't appear to add anything useful here, but increasing the
// krylov length has a useful effect to avoid stagnation
template <typename Func>
void
LinearSolve(MatrixOperator& x, Func F, MatrixOperator const& Rhs, double Tol = 1E-14, int Verbose = 0)
{
   int m = 30;     // krylov subspace size
   int iter = 0;   // total number of iterations performed

   double normb = norm_frob(Rhs);

   int IterThisRound = m*50; // if it takes more than 50 rounds to converge, then we've proably stagnated
   double tol = Tol;
   int Ret = GmRes(x, F, normb, Rhs, m, IterThisRound, tol, LinearAlgebra::Identity<MatrixOperator>(), Verbose);
   iter += IterThisRound;

   while (Ret != 0 && iter < MaxIter)
   {
      // Attempt to avoid stagnation by increasing the number of iterations
      m=m+10;
      if (Verbose > 1)
      {
	 std::cerr << "Refinement step, increasing m to " << m << '\n';
      }

      //      TRACE("Refinement step")(iter);
      // iterative refinement step
#if defined(USE_ITERATIVE_REFINEMENT)
      MatrixOperator R = Rhs- F(x);
      MatrixOperator xRefine = R;
#else
      MatrixOperator R = Rhs;
      MatrixOperator xRefine = x;
#endif
      IterThisRound = m*50;
      double tol = Tol;
      Ret = GmRes(xRefine, F, normb, R, m, IterThisRound, tol, LinearAlgebra::Identity<MatrixOperator>(), Verbose);

      iter += IterThisRound;
#if defined(USE_ITERATIVE_REFINEMENT)
      x += xRefine;
#else
      x = xRefine;
#endif

      if (Verbose > 1)
      {
	 double Resid = norm_frob(F(x) - Rhs) / normb;
	 std::cerr << "Residual after refinement step = " << Resid << '\n';
      }
   }

   if (Ret != 0)
   {
      // failed
      PANIC("Linear solver failed to converge after max_iter iterations")(MaxIter);
   }
}

// Calculate the (complex) eigenvalue that is closest to 1.0
// using Arnoldi.
template <typename T>
std::complex<double>
FindClosestUnitEigenvalue(MatrixOperator& M, T Func, double tol, int Verbose)
{
   int Iterations = 20;
   double Tol = tol;
   std::complex<double> EtaL;
   EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, Verbose);
   while (Iterations == 20)
   {
      Tol = Tol;
      EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, Verbose);
   }
   return EtaL;
}

KComplexPolyType
DecomposeParallelParts(KMatrixPolyType& C, std::complex<double> Factor,
                       MatrixOperator const& UnitMatrixLeft,
                       MatrixOperator const& UnitMatrixRight, double UnityEpsilon)
{
   KComplexPolyType EParallel;
   // diagonal element is the identity, up to a unitary factor
   DEBUG_TRACE("Unit diagonal element")(Factor);

   // decompose C into components parallel and perpendicular to the identity
   // The only part we have to care about is a component with the same momentum as our unit operator
   for (KMatrixPolyType::iterator Ck = C.begin(); Ck != C.end(); ++Ck) // sum over momenta
   {
      ComplexPolyType CParallel;
      std::complex<double> K = Ck->first;  // the momentum (complex phase)
      for (MatrixPolyType::iterator I = Ck->second.begin(); I != Ck->second.end(); ++I)
      {
         if (I->second.is_null())
         {
            TRACE(Ck->second);
            PANIC("oops");
         }

         //if (!is_scalar(I->second.TransformsAs()))
         //continue;

         std::complex<double> Overlap = inner_prod(I->second, UnitMatrixRight);
         DEBUG_TRACE(Overlap);
         I->second -= conj(Overlap)*UnitMatrixLeft;
         DEBUG_TRACE(inner_prod(I->second, UnitMatrixRight))("should be zero");
         DEBUG_TRACE(inner_prod(UnitMatrixLeft, UnitMatrixRight));

         // **comparison** we always want to add here to get the degree of the polynomial correct.
         // This is the important one
         //      if (norm_frob(Overlap) > 1E-16)

            CParallel[I->first] = Overlap;
      }

      // Is this the same momentum as our unit operator?
      if (norm_frob(K - Factor) < UnityEpsilon*10)
      {
         DEBUG_TRACE("Component at equal momenta")(K);
         // same momenta, these components diverge
         for (int m = CParallel.degree(); m >= 0; --m)
         {
            if (CParallel.has_term(m))
            {
               EParallel[Factor][m+1] = conj(Factor) * CParallel[m]; // include momentum
               for (int k = m+2; k <= CParallel.degree()+1; ++k)
               {
                  if (EParallel[Factor].has_term(k))
                  {
                     EParallel[Factor][m+1] -= double(Binomial(k,m))
                        * EParallel[Factor][k];
                  }
               }
               EParallel[Factor][m+1] *= 1.0 / (1.0 + m);
            }
         }
      }
      else
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

KMatrixPolyType
DecomposePerpendicularParts(KMatrixPolyType& C,
                            BasicFiniteMPO const& Diag,
                            MatrixOperator const& UnitMatrixLeft,
                            MatrixOperator const& UnitMatrixRight,
                            LinearWavefunction const& Psi,
                            QuantumNumber const& QShift,
                            bool HasEigenvalue1,
                            double Tol,
                            int Verbose)
{
   // Components perpendicular to the identity satisfy equation (24)
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      std::complex<double> K = I->first;  // the momentum (complex phase)

      if (Verbose > 0)
         std::cerr << "Momentum " << K << std::endl;

      DEBUG_TRACE("Momentum")(K)(HasEigenvalue1);
      for (int m = I->second.degree(); m >= 0; --m)
      {
         if (Verbose > 0)
            std::cerr << "Degree " << m << std::endl;

         DEBUG_TRACE("degree")(m);
         MatrixOperator Rhs = conj(K) * I->second[m];
         for (int k = m+1; k <= I->second.degree(); ++k)
         {
            // avoid accessing E[K][k] if it doesn't exist, to avoid adding a null term
            if (E.has_element(K) && E[K].has_term(k))
               Rhs -= double(Binomial(k,m)) * E[K][k];
         }

         // orthogonalize Rhs against the identity again, which is a kind of
         // DGKS correction
         if (HasEigenvalue1 && Rhs.TransformsAs() == UnitMatrixRight.TransformsAs())
         {
            DEBUG_TRACE(inner_prod(Rhs, UnitMatrixRight))("should be small");
            Rhs -= conj(inner_prod(Rhs, UnitMatrixRight)) * UnitMatrixLeft;
            DEBUG_TRACE(inner_prod(Rhs, UnitMatrixRight))("should be zero");
            DEBUG_TRACE(inner_prod(Rhs, UnitMatrixLeft));
         }


         double RhsNorm2 = norm_frob_sq(Rhs);
         RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
         DEBUG_TRACE(RhsNorm2);
         //      DEBUG_TRACE(HasEigenvalue1)(UnitMatrixLeft)(UnitMatrixRight)(K)(Diag)(Rhs);

         //if (RhsNorm2 > 1E-22)

         {
            //TRACE(norm_frob_sq(Rhs))(inner_prod(Rhs, UnitMatrixRight));
            // Initial guess vector -- scale it by the norm of Rhs, improves the stability
            E[K][m] = std::sqrt(RhsNorm2) *
               MakeRandomMatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
            // Orthogonalize the initial guess -- this is important for the numerical stability
            if (HasEigenvalue1 && Rhs.TransformsAs() == UnitMatrixRight.TransformsAs())
            {
               E[K][m] -= conj(inner_prod(E[K][m], UnitMatrixRight)) * UnitMatrixLeft;
            }

            LinearSolve(E[K][m], OneMinusTransferLeft_Ortho(K*Diag, Psi, QShift,
                                                            UnitMatrixLeft, UnitMatrixRight, HasEigenvalue1),
                        Rhs, Tol, Verbose);

            // do another orthogonalization -- this should be unncessary but for the paranoid...
            if (HasEigenvalue1 && E[K][m].TransformsAs() == UnitMatrixRight.TransformsAs())
            {
               std::complex<double> z = inner_prod(E[K][m], UnitMatrixRight);
               DEBUG_TRACE(z);
               if (LinearAlgebra::norm_frob_sq(z) > 1E-10)
               {
                  WARNING("Possible numerical instability in triangular MPO solver")(z);
               };
               E[K][m] -= conj(z) * UnitMatrixLeft;
               DEBUG_TRACE(inner_prod(E[K][m], UnitMatrixRight))("should be zero");
            }
         }
      }
   }
   return E;
}

// Solve the components for the case where the diagonal operator is zero
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
      std::complex<double> K = I->first;  // the momentum (complex phase)
      int MaxDegree = I->second.degree();
      for (int i = MaxDegree; i >= 0; --i)
      {
         if (I->second.has_term(i))
         {
            DEBUG_CHECK(!I->second[i].is_null());
            E[K][i] = conj(K) * I->second[i];
         }
         for (int j = i+1; j <= MaxDegree; ++j)
         {
            DEBUG_CHECK(!E[K][j].is_null());
            E[K][i] -= double(Binomial(j,i)) * E.lookup_or_default(K)[j];
         }
      }
   }
   return E;
}


// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator

void
SolveMPO_Left(std::vector<KMatrixPolyType>& EMatK,
              LinearWavefunction const& Psi, QuantumNumber const& QShift,
              BasicTriangularMPO const& Op, MatrixOperator const& LeftIdentity,
              MatrixOperator const& RightIdentity, bool NeedFinalMatrix,
              double Tol,
              double UnityEpsilon, int Verbose)
{
   CHECK_EQUAL(RightIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(RightIdentity.Basis2(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis2(), Psi.Basis1());

   int Dim = Op.Basis1().size();       // dimension of the MPO
   EMatK.reserve(Dim);

   if (Verbose > 0)
      std::cerr << "SolveMPO_Left: dimension is " << Dim << std::endl;

   if (EMatK.empty())
   {
      // Make sure the (0,0) part is identity
      DEBUG_TRACE(UnityEpsilon);
      OperatorClassification CheckIdent = classify(Op(0,0), UnityEpsilon);
      if (!CheckIdent.is_identity())
      {
         std::cerr << "SolveMPO_Left: fatal: (0,0) component of the MPO must be the identity operator.\n";
         // the (0,0) component isn't the identity operator, which is a fatal error.
         // Show some diagnosics and quit.
         if (!CheckIdent.is_product())
         {
            std::cerr << "SolveMPO_Left: fatal: component is not a product operator!\n";
         }
         else
         {
            if (CheckIdent.is_unitary())
               std::cerr << "SolveMPO_Left: fatal: component is a unitary operator, but not the identity.\n";
            else if (CheckIdent.is_prop_identity())
            {
               std::cerr << "SolveMPO_Left: fatal: component is proportional "
                         << "to the identity, with prefactor " << CheckIdent.factor() << '\n';
            }
            else
               std::cerr << "SolveMPO_Left: fatal: component has unknown classification.\n";
         }
         PANIC("Fatal")(CheckIdent)(Op(0,0));
      }

      // Initialize the first E matrix.  These are operators acting in the Basis1()
      EMatK.push_back(KMatrixPolyType());
      EMatK[0][1.0] = MatrixPolyType(LeftIdentity);
   }

   int StartCol = EMatK.size();

   // fill out the remainder of EMat with zero
   while (int(EMatK.size()) < Dim)
      EMatK.push_back(KMatrixPolyType());

   // solve recursively column 1 onwards
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving column " << Col << " of " << (Dim-1) << '\n';
      }

      // Generate the next C matrices, C(n) = sum_{j<Col} Op(j,Col) E_j(n)
      KMatrixPolyType C;

      std::vector<std::vector<int> > Mask = mask_column(Op, Col);
      C = inject_left_mask(EMatK, Psi, QShift, Op.data(), Psi, Mask)[Col];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;
         EMatK[Col] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

         // Non-zero diagonal element.
         // In this case we have to consider a possible component with eigenvalue magnitude 1.
         // We call this component (if it exists) the unit matrix.

         KComplexPolyType EParallel;  // components parallel to the identity at some momentum, may be zero

         // Multiplication factor and the left and right eigenvalues.
         std::complex<double> Factor = Classification.factor();
         MatrixOperator UnitMatrixLeft = LeftIdentity;
         MatrixOperator UnitMatrixRight = RightIdentity; //delta_shift(RightIdentity, QShift);

         // If the diagonal operator is unitary, it might have an eigenvalue of magnitude 1
         // Or if we are not in the scalar sector, then there might be an additional eigenvalue 1
         // due to symmetry.  We don't attempt to handle the case where we have more than one
         // eigenvalue 1 in the same symmetry sector.
#if defined(NDEBUG)
         if (Classification.is_unitary() && (!is_scalar(Diag.Basis2()[0]) || !Classification.is_complex_identity()))
#else
         // if we're debugging, then find the eigenvector anyway, even if its the identity
         if (Classification.is_unitary())
#endif
         {
            if (Verbose > 0)
            {
               if (Classification.is_complex_identity())
               {
               std::cerr << "Identity in non-scalar sector at " << Col
                         << ", finding left eigenvalue..." << std::endl;
               }
               else
               {
                  std::cerr << "Non-trivial unitary at column " << Col
                            << ", finding left eigenvalue..." << std::endl;
               }
            }

            // Find the largest eigenvalue
            // We need initial guess vectors in the correct symmetry sector
            UnitMatrixLeft = MakeRandomMatrixOperator(LeftIdentity.Basis1(), LeftIdentity.Basis2(), Diag.Basis2()[0]);

            std::complex<double> EtaL = FindClosestUnitEigenvalue(UnitMatrixLeft,
                                                                  InjectLeftQShift(Diag, Psi, QShift),
                                                                  Tol, Verbose);
            EtaL = conj(EtaL); // left eigenvalue, so conjugate (see comment at operator_actions.h)
            if (Verbose > 0)
               std::cerr << "Eigenvalue of unitary operator is " << EtaL << std::endl;
            Factor = EtaL;

            if (std::abs(norm_frob(EtaL) - 1.0) < UnityEpsilon)
            {
               if (Verbose > 0)
                  std::cerr << "Found an eigenvalue 1, so we need the right eigenvector..." << std::endl;

               UnitMatrixRight = UnitMatrixLeft;
               // we have an eigenvalue of magnitude 1.  Find the right eigenvalue too
               std::complex<double> EtaR = FindClosestUnitEigenvalue(UnitMatrixRight,
                                                                     InjectRightQShift(Diag, Psi,
                                                                                       QShift),
                                                                     Tol, Verbose);
               if (Verbose > 0)
                  std::cerr << "Right eigenvalue is " << EtaR << std::endl;

               CHECK(norm_frob(EtaL-EtaR) < UnityEpsilon)("Left and right eigenvalues do not agree!")(EtaL)(EtaR);
               // we already determined that the norm is sufficiently close to 1, but
               // fine-tune normalization
               Factor = EtaL / norm_frob(EtaL);

               // normalize the left/right eigenvector pair
               //              UnitMatrixLeft *= 1.0 / (inner_prod(UnitMatrixRight, UnitMatrixLeft));
               UnitMatrixRight *= 1.0 / (inner_prod(UnitMatrixLeft, UnitMatrixRight));
               DEBUG_TRACE(inner_prod(UnitMatrixLeft, UnitMatrixRight));
            }
            else if (std::abs(norm_frob(EtaL) - 1.0) < UnityEpsilon*100)
            {
               // If we get here, then we have an eigenvalue that is close to 1, but
               // misses the test
               std::cerr << "SolveMPO_Left: warning: Eigenvalue misses tolerance but is near 1"
                         << ", epsilon=" << std::abs(norm_frob(EtaL)-1.0) << '\n';
            }
         }
         else if (Verbose && Classification.is_identity())
         {
            std::cerr << "Diagonal component is the identity\n";
         }
         else if (Verbose && Classification.is_complex_identity())
         {
            std::cerr << "Diagonal component is proportional to the identity\n";
         }

         // If we have an eigenvalue equal to 1, then decompose C into parallel and perpendicular parts
         bool HasEigenvalue1 = false;
         if (std::abs(norm_frob(Factor) - 1.0) < UnityEpsilon)
         {
            HasEigenvalue1 = true;
            //DEBUG_TRACE(UnitMatrixLeft)(UnitMatrixRight);
            if (Verbose > 0)
               std::cerr << "Decomposing parts parallel to the unit matrix\n";
            EParallel = DecomposeParallelParts(C, Factor, UnitMatrixLeft, UnitMatrixRight, UnityEpsilon);
         }
         else
         {
            if (Verbose > 0)
            {
               std::cerr << "Diagonal component is not unitary, assuming spectral radius < 1\n";
               DEBUG_TRACE(Diag)(Classification);
            }
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1

         KMatrixPolyType E;
         // if we are on the last column and we don't need the matrix elements, then we can
         // skip this operation
         if (Col < Dim-1 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix\n";
            E = DecomposePerpendicularParts(C, Diag, UnitMatrixLeft, UnitMatrixRight,
                                            Psi, QShift, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
         {
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last column.\n";
         }

         // Reinsert the components parallel to the unit matrix (if any)
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
            for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
            {
               // Conj here because this comes from an overlap(x, RightUnitMatrix)
               E[I->first][J->first] += conj(J->second) * UnitMatrixLeft;
               DEBUG_TRACE(inner_prod(E[I->first][J->first], UnitMatrixRight));
            }
         }

         // Finally, set the E matrix element at this column
         //DEBUG_TRACE(E[1.0]);
         EMatK[Col] = E;
      }

   }
}

//
// SolveSimpleMPO_Left
//

struct SubProductLeftProject
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductLeftProject(LinearWavefunction const& Psi_, QuantumNumber const& QShift_,
                         MatrixOperator const& Proj_, MatrixOperator const& Ident_)
      : Psi(Psi_), QShift(QShift_), Proj(Proj_),
        Ident(Ident_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In; //delta_shift(In, QShift);
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
       {
          Result = operator_prod(herm(*I), Result, *I);
       }
      Result = In - delta_shift(Result, QShift);
      //Result = 0.5 * (Result + adjoint(Result));
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator Proj;
   MatrixOperator Ident;
};

std::complex<double>
SolveSimpleMPO_Left(StateComponent& E, LinearWavefunction const& Psi,
                    QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                    MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (E.is_null())
      E = Initial_E(Op, Psi.Basis1());
   CHECK_EQUAL(E.Basis1(), Psi.Basis1());
   CHECK_EQUAL(E.Basis2(), Psi.Basis1());
   CHECK_EQUAL(E.LocalBasis(), Op.Basis1());
   CHECK_EQUAL(Rho.Basis1(), Psi.Basis1());
   CHECK_EQUAL(Rho.Basis2(), Psi.Basis1());

   MatrixOperator Ident = E[0];

   // The UnityEpsilon is just a paranoid check here, as we don't support
   // eigenvalue 1 on the diagonal
   double UnityEpsilon = 1E-12;

   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Left: fatal: MPO(0,0) must be the identity operator.\n";
      PANIC("Fatal");
   }

   int Dim = Op.Basis1().size();       // dimension of the MPO
   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Left: dimension is " << Dim << std::endl;

   // Column 0 (E[0]) is the Identity, we don't need to solve for it
   for (int Col = 1; Col < Dim-1; ++Col)
   {
      std::vector<std::vector<int> > Mask = mask_column(Op, Col);
      MatrixOperator C = inject_left_mask(E, Psi, Op.data(), Psi, Mask)[Col];
      C = delta_shift(C, QShift);

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;
         E[Col] = C;
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

         // non-zero diagonal element.  The only case that we support here is
         // an operator with spectral radius strictly < 1
         if (Classification.is_unitary())
         {
            std::cerr << "SolveSimpleMPO_Left: Unitary operator on the diagonal is not supported!\n";
            PANIC("Fatal: unitary")(Col);
         }

         // Initial guess for linear solver
         E[Col] = C;

	 if (Col == 2)
	    TRACE(C);

	 LinearSolve(E[Col], OneMinusTransferLeft_Ortho(Diag, Psi, QShift, Ident, Rho, false), C, Tol, Verbose);
         //LinearSolve(E[Col], OneMinusTransferLeft(Diag, Psi, QShift), C, Ident, Rho, Tol, Verbose);
      }
   }
   // Final column, must be identity
   int const Col = Dim-1;
   if (!classify(Op(Col,Col), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Left: fatal: MPO(d,d) must be the identity operator.\n";
      PANIC("Fatal");
   }

   std::vector<std::vector<int> > Mask = mask_column(Op, Col);
   MatrixOperator C = inject_left_mask(E, Psi, Op.data(), Psi, Mask)[Col];
   C = delta_shift(C, QShift);

   // The component in the direction of the identity is proportional to the energy
   std::complex<double> Energy = inner_prod(Rho, C);

   // orthogonalize
   C -= Energy * Ident;

   // solve for the first component
   SubProductLeftProject ProdL(Psi, QShift, Rho, Ident);
   E[Col] = C;
   LinearSolve(E[Col], ProdL, C, Tol, Verbose);

   // Make it Hermitian
   //   E.back() = 0.5 * (E.back() + adjoint(E.back()));

   // Stability fix: remove overall constant
   if (Verbose > 0)
      std::cerr << "Overall constant: " << inner_prod(E.back(), E.front()) << '\n';
   E.back() -= inner_prod(E.front(), E.back()) * E.front();

   // remove the spurious constant term from the energy
   if (Verbose > 0)
      std::cerr << "Spurious constant: " << inner_prod(E.back(), Rho) << '\n';
   E.back() -= inner_prod(Rho, E.back()) * E.front();


#if 0
   // residual
   MatrixOperator R = E.back();
   for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
      R = operator_prod(herm(*I), R, *I);
   }
   R = delta_shift(R, QShift);
   R += C;

   DEBUG_TRACE("Residual norm")(norm_frob(E.back() - R));

   E.back() = R;

   // Make it Hermitian
   E.back() = 0.5 * (E.back() + adjoint(E.back()));
#endif


#if !defined(NDEBUG)
   {
      std::vector<KMatrixPolyType> CheckEMat;
      SolveMPO_Left(CheckEMat, Psi, QShift, Op, Ident, Rho, true);
      ComplexPolyType EValues = ExtractOverlap(CheckEMat.back()[1.0], Rho);
      TRACE(EValues);
      MatrixOperator HCheck = CheckEMat.back()[1.0][0];
      TRACE("H matrix elements check")(inner_prod(HCheck - E.back(), Rho));
      TRACE(inner_prod(HCheck, Rho));
   }
#endif

   return Energy;
}

std::complex<double>
SolveSimpleMPO_Left(StateComponent& E, InfiniteWavefunctionLeft const& Psi,
                    BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(PsiLinear, Lambda) = get_left_canonical(Psi);
   MatrixOperator Rho = Lambda*Lambda;

   return SolveSimpleMPO_Left(E, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}

//
// SolveSimpleMPO_Right
//

struct SubProductRightProject
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductRightProject(LinearWavefunction const& Psi_, QuantumNumber const& QShift_,
                          MatrixOperator const& Proj_, MatrixOperator const& Ident_)
      : Psi(Psi_), QShift(QShift_), Proj(Proj_), Ident(Ident_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In;
      LinearWavefunction::const_iterator I = Psi.end();
      while (I != Psi.begin())
      {
         --I;
         Result = operator_prod(*I, Result, herm(*I));
      }
      Result = delta_shift(Result, adjoint(QShift));
      Result = In - Result;
      //Result = 0.5 * (Result + adjoint(Result));
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator const& Proj;
   MatrixOperator const& Ident;
};

std::complex<double>
SolveSimpleMPO_Right(StateComponent& F, LinearWavefunction const& Psi,
                    QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                    MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (F.is_null())
      F = Initial_F(Op, Psi.Basis2());
   CHECK_EQUAL(F.Basis1(), Psi.Basis2());
   CHECK_EQUAL(F.Basis2(), Psi.Basis2());
   CHECK_EQUAL(F.LocalBasis(), Op.Basis1());
   CHECK_EQUAL(Rho.Basis1(), Psi.Basis2());
   CHECK_EQUAL(Rho.Basis2(), Psi.Basis2());

   MatrixOperator Ident = F.back();

   // The UnityEpsilon is just a paranoid check here, as we don't support
   // eigenvalue 1 on the diagonal
   double UnityEpsilon = 1E-12;

   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Right: fatal: MPO(0,0) must be the identity operator.\n";
      PANIC("Fatal");
   }

   int Dim = Op.Basis1().size();       // dimension of the MPO
   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Right: dimension is " << Dim << std::endl;

   // Row Dim-1 (F[Dim-1]) is the Identity, we don't need to solve for it
   for (int Row = Dim-2; Row >= 1; --Row)
   {
      std::vector<std::vector<int> > Mask = mask_row(Op, Row);
      MatrixOperator C = inject_right_mask(F, Psi, Op.data(), Psi, Mask)[Row];
      C.delta_shift(adjoint(QShift));

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Row)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at row " << Row << std::endl;
         F[Row] = C;
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at row " << Row << std::endl;

         // non-zero diagonal element.  The only case that we support here is
         // an operator with spectral radius strictly < 1
         if (Classification.is_unitary())
         {
            std::cerr << "SolveSimpleMPO_Right: Unitary operator on the diagonal is not supported!\n";
            PANIC("Fatal: unitary")(Row);
         }

         // Initial guess for linear solver
         F[Row] = C;

	 LinearSolve(F[Row], OneMinusTransferRight(Diag, Psi, QShift), C, Tol, Verbose);
      }
   }
   // Final row, must be identity
   int const Row = 0;
   if (!classify(Op(Row,Row), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Right: fatal: MPO(d,d) must be the identity operator.\n";
      PANIC("Fatal");
   }

   std::vector<std::vector<int> > Mask = mask_row(Op, Row);
   MatrixOperator C = inject_right_mask(F, Psi, Op.data(), Psi, Mask)[Row];
   C.delta_shift(adjoint(QShift));

   // The component in the direction of the identity is proportional to the energy
   std::complex<double> Energy = inner_prod(Rho, C);
   // orthogonalize
   C -= Energy * Ident;

   // solve for the first component
   SubProductRightProject ProdR(Psi, QShift, Rho, Ident);
   F[Row] = C;
   LinearSolve(F[Row], ProdR, C, Tol, Verbose);

   // Make it Hermitian
   //   F.front() = 0.5 * (F.front() + adjoint(F.front()));

   // stability fix
   if (Verbose > 0)
      std::cerr << "Overall constant " << inner_prod(F.front(), F.back()) << '\n';
   F.front() -= inner_prod(F.back(), F.front()) * F.back();

   // remove the spurious constant term from the energy
   if (Verbose > 0)
      std::cerr << "Spurius constant " << inner_prod(F.front(), Rho) << '\n';
   F.front() -= inner_prod(Rho, F.front()) * F.back();

   // Everything here is in the Hermitian representation, so the actual energy is
   // the conjugate
   return conj(Energy);
}

std::complex<double>
SolveSimpleMPO_Right(StateComponent& F, InfiniteWavefunctionRight const& Psi,
                     BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(Lambda, PsiLinear) = get_right_canonical(Psi);
   MatrixOperator Rho = Lambda*Lambda;
   return SolveSimpleMPO_Right(F, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}
