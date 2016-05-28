// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iexpectation-tri.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "mpo/triangular_operator.h"
#include "mps/infinitewavefunction.h"
#include "mps/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"
#include "common/polynomial.h"
#include "tensor/tensor_eigen.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"
#include "models/hubbard-u1u1.h"
#include "models/bosehubbard-spinless.h"
#include "models/bosehubbard-spinless-u1.h"

long Binomial(int n, int k)
{
   if (k > n/2)
      k = n-k;     // take advantage of symmetry
   double r = 1.0;
   for (int i = 1; i <= k; ++i)
   {
      r *= double(n-k+i) / double(i);
   }
   return long(r+0.5); // round to nearest
}

// Functor to evaluate (1-T)(Psi) where T is the generalized transfer matrix
struct OneMinusTransfer
{
   OneMinusTransfer(MPOperator const& T, LinearWavefunction const& Psi, QuantumNumber const& QShift, 
                    MatrixOperator const& Rho,
                    MatrixOperator const& Identity, bool Orthogonalize)
      : T_(T), Psi_(Psi), 
	QShift_(QShift), Rho_(Rho), 
	Identity_(Identity), Orthogonalize_(Orthogonalize) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-delta_shift(apply_right(x, T_, Psi_), QShift_);
      //      r = delta_shift(r, QShift_);
      if (Orthogonalize_ && is_scalar(r.TransformsAs()))
          r -= inner_prod(r, Rho_) * Identity_; // orthogonalize to the identity
      return r;
   }

   MPOperator const& T_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
   MatrixOperator const& Rho_;
   MatrixOperator const& Identity_;
   bool Orthogonalize_;
};

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator Rhs)
{
   MatrixOperator Guess = Rhs;
   int m = 30;
   int max_iter = 10000;
   double tol = 1e-15;
   GmRes(Guess, F, Rhs, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());
   return Guess;
}

template <typename T>
struct EigenPair
{
   T LeftVector;
   T RightVector;
   std::complex<double> Eigenvalue;
};

// Comparitor for complex numbers.  This is so that we can put them in a map,
// the choice of comparison operation is arbitrary
struct CompareComplex
{
   typedef std::complex<double> first_argument_type;
   typedef std::complex<double> second_argument_type;
   typedef bool result_type;
   bool operator()(std::complex<double> const& x, std::complex<double> const& y) const
   {
      return (x.real() < y.real()) || (x.real() == y.real() && x.imag() < y.imag());
   }
};



// polynomial with matrix coefficients
typedef Polynomial<MatrixOperator> MatrixPolyType;

// polynomial with complex coefficients
typedef Polynomial<std::complex<double> > ComplexPolyType;

// Momentum-dependent complex polynomial
typedef std::map<std::complex<double>, ComplexPolyType, CompareComplex> KComplexPolyType;

// momentum-dependent matrix polynomial,
// this represents an E matrix
typedef std::map<std::complex<double>, MatrixPolyType, CompareComplex> KMatrixPolyType;

MatrixPolyType
delta_shift(MatrixPolyType const& In, QuantumNumber const& QShift)
{
   MatrixPolyType Result(In);
   for (MatrixPolyType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

#if 0
// define here.
// Or better as vector<Polynomial<MatrixOperator> > ?
// Probably even better, to do finite momentum at the same time.
std::vector<MatrixPolyType>
apply_right(std::vector<MatrixPolyType> const& In, 
            LinearWavefunction const& Psi1, 
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());

   std::vector<MatrixPolyType> Result(Op.Basis2().size());
   int MaxDegree = 0;
   for (unsigned i = 0; i < In.size(); ++i)
      MaxDegree = std::max(In[i].degree(), MaxDegree);

   for (int Degree = 0; Degree <= MaxDegree; ++Degree)
   {
      StateComponent E(Op.Basis1(), Psi1.Basis1(), Psi2.Basis1());
      for (unsigned k = 0; k < E.size(); ++k)
      {
         E[k] = In[k][Degree];
      }

      E = apply_right(E, Psi1, Op, Psi2);

      CHECK_EQUAL(E.size(), Result.size());
      for (unsigned i = 0; i < Result.size(); ++i)
      {
         Result[i][Degree] += E[i];
      }
   }
   return Result;
}

MatrixPolyType
apply_right(MatrixPolyType const& In, 
            LinearWavefunction const& Psi1, 
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   std::vector<MatrixPolyType> Vec(1, In);
   Vec = apply_right(Vec, Psi1, Op, Psi2);
   CHECK_EQUAL(Vec.size(), 1);
   return Vec[0];
}

// Calculates result' = C[column] = sum_{j > Column} E[j] * T_Op(j, Column)
// Assumes that E[j] is defined, for j > Column
MatrixPolyType
MultiplyLeft(std::vector<MatrixPolyType> const& E, 
             TriangularOperator const& Op, 
             LinearWavefunction const& Psi, 
             QuantumNumber const& QShift, int Column)
{
   CHECK_EQUAL(Op.size(), Psi.size());
   MatrixPolyType Result;

   // replace this stuff with the apply_right implementation, and extract_column() in TriangularOperator

   MPOperator OpCol = extract_column(Op, Column);
   std::vector<MatrixPolyType> C = apply_right(E, Psi, OpCol, Psi);
   return delta_shift(C[Column], QShift);
}
#endif


KMatrixPolyType
delta_shift(KMatrixPolyType const& In, QuantumNumber const& QShift)
{
   KMatrixPolyType Result(In);
   for (KMatrixPolyType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

std::vector<KMatrixPolyType>
delta_shift(std::vector<KMatrixPolyType> const& In, QuantumNumber const& QShift)
{
   std::vector<KMatrixPolyType> Result(In);
   for (std::vector<KMatrixPolyType>::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = delta_shift(*I, QShift);
   }
   return Result;
}

void
add_triple_prod(KMatrixPolyType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                KMatrixPolyType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp)
{
   // loop over momenta
   for (KMatrixPolyType::const_iterator K = E.begin(); K != E.end(); ++K)
   {
      // loop over degrees of the polynomial
      for (MatrixPolyType::const_iterator D = K->second.begin(); D != K->second.end(); ++D)
      {
         Result[K->first][D->first] += Factor * triple_prod(x, D->second, y, qxy, qEp);
      }
   }
}

std::vector<KMatrixPolyType>
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              std::vector<KMatrixPolyType> const& E, 
              StateComponent const& B)
{
   //   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().Basis1());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1().size(), E.size());
   
   std::vector<KMatrixPolyType> Result(M.base().Basis2().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index2()], herm(*S), 
                                  herm(A.base()[S.index1()]), 
                                  E[J.index1()], 
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.base().Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<KMatrixPolyType>
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              std::vector<KMatrixPolyType> const& E, 
              StateComponent const& B,
              std::vector<int> const& OutMask,
              std::vector<int> const& InMask)
{
   //   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1().size(), E.size());
   
   std::vector<KMatrixPolyType> Result(M.base().Basis2().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // skip over masked components
      if (!InMask[I.index()])
         continue;

      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // skip over masked components
         if (!OutMask[J.index2()])
            continue;

         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index2()], herm(*S), 
                                  herm(A.base()[S.index1()]), 
                                  E[J.index1()], 
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.base().Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<KMatrixPolyType>
apply_right(std::vector<KMatrixPolyType> const& In, 
            LinearWavefunction const& Psi1, 
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   PRECONDITION_EQUAL(Op.Basis1().size(), In.size());

   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   MPOperator::const_iterator OpIter = Op.begin();

   std::vector<KMatrixPolyType> E;
   std::vector<KMatrixPolyType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      //      std::vector<KMatrixPolyType>(OpIter->Basis2().size()).swap(Result);

      Result = operator_prod(herm(*OpIter), herm(*I1), E, *I2);

      ++I1; ++I2; ++OpIter;
   }
   return Result;
}

std::vector<KMatrixPolyType>
apply_right_mask(std::vector<KMatrixPolyType> const& In, 
                 LinearWavefunction const& Psi1, 
                 QuantumNumber const& QShift,
                 MPOperator const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   PRECONDITION_EQUAL(Op.Basis1().size(), In.size());

   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   MPOperator::const_iterator OpIter = Op.begin();
   std::vector<std::vector<int> >::const_iterator MaskIter = Mask.begin();

   std::vector<KMatrixPolyType> E;
   std::vector<KMatrixPolyType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      //      std::vector<KMatrixPolyType>(OpIter->Basis2().size()).swap(Result);

      Result = operator_prod(herm(*OpIter), herm(*I1), E, *I2, *(MaskIter+1), *MaskIter);

      ++I1; ++I2; ++OpIter; ++MaskIter;
   }
   return delta_shift(Result, QShift);
}

template <typename T>
std::complex<double>
FindLargestEigenvalue(MatrixOperator& M, T Func, bool Verbose)
{
   int Iterations = 20;
   double Tol = 1E-14;
   std::complex<double> EtaL;
   EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, Verbose);
   while (Iterations == 20)
   {
      Tol = 1E-14;
      EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, Verbose);
   }
   return EtaL;
}

KComplexPolyType
DecomposeParallelParts(KMatrixPolyType& C, std::complex<double> Factor, 
		       MatrixOperator const& UnitMatrixLeft, 
		       MatrixOperator const& UnitMatrixRight)
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
	 if (!is_scalar(I->second.TransformsAs()))
	    continue;

	 std::complex<double> Overlap = inner_prod(I->second, UnitMatrixRight);
	 I->second -= Overlap*UnitMatrixLeft;
	 if (norm_frob(Overlap) > 1E-16)
	    CParallel[I->first] = Overlap;
      }

      // Is this the same momentum as our unit operator?
      if (norm_frob(K - Factor) < 1E-12)
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
			    MPOperator const& Diag, 
			    MatrixOperator const& UnitMatrixLeft, 
			    MatrixOperator const& UnitMatrixRight, 
			    LinearWavefunction const& Psi,
			    QuantumNumber const& QShift,
			    bool Magnitude1)
{
   // Components perpendicular to the identity satisfy equation (24)
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      std::complex<double> K = I->first;  // the momentum (complex phase)
      DEBUG_TRACE("Momentum")(K);
      for (int m = I->second.degree(); m >= 0; --m)
      {
	 DEBUG_TRACE("degree")(m);
	 MatrixOperator Rhs = conj(K) * I->second[m];
	 for (int k = m+1; k <= I->second.degree(); ++k)
	 {
	    Rhs -= double(Binomial(k,m)) * E[K][k];
	 }
            
	 // orthogonalize Rhs against the identity again, which is a kind of
	 // DGKS correction
	 if (is_scalar(Rhs.TransformsAs()))
	    Rhs -= inner_prod(UnitMatrixRight, Rhs) * UnitMatrixLeft;

	 double RhsNorm2 = norm_frob_sq(Rhs);
	 RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
	 //TRACE(RhsNorm2);
	 if (RhsNorm2 > 1E-22)
	 {
	    E[K][m] = LinearSolve(OneMinusTransfer(K*Diag, Psi, QShift, UnitMatrixRight, UnitMatrixLeft, Magnitude1), 
				  Rhs);
	    // do another orthogonalization
	    if (Magnitude1 && is_scalar(E[K][m].TransformsAs()))
	       E[K][m] -= inner_prod(UnitMatrixRight, E[K][m]) * UnitMatrixLeft;
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
	    E[K][i] = conj(K) * I->second[i];
	 for (int j = i+1; j <= MaxDegree; ++j)
	 {
	    E[K][i] -= double(Binomial(j,i)) * E[K][j];
	 }
      }
   }
   return E;
}

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator
KMatrixPolyType
SolveMPO_Left(LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularOperator const& Op, MatrixOperator const& Rho,
              MatrixOperator const& Identity, bool Verbose = false)
{
   int Dim = Op.Basis1().size();       // dimension of the MPO
   std::vector<KMatrixPolyType> EMatK(Dim);  // the vector of E matrices

   // Initialize the first E matrix
   EMatK[Dim-1][1.0] = MatrixPolyType(Identity);

   // solve recursively
   int Col = Dim-2;
   while (Col >= 0)
   {
      if (Verbose)
      {
         std::cerr << "Solving column " << Col << '\n';
      }

      // Generate the next C matrices, C(n) = sum_{j>Col} Op(j,Col) E_j(n)
      KMatrixPolyType C;

#if 0
      MPOperator M = extract_lower_column(Op, Col);
      C = apply_right(EMatK, Psi, M, Psi)[0];
#else
      std::vector<std::vector<int> > Mask;
      mask_lower_column(Op, Col, Mask);
      C = apply_right_mask(EMatK, Psi, QShift, Op.data(), Psi, Mask)[Col];
#endif

      // Now do the classification, based on the properties of the diagonal operator
      MPOperator Diag = Op(Col, Col);
      MPOperatorClassification Classification = classify(Diag);
      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
	 if (Verbose)
	    std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;
         EMatK[Col] = SolveZeroDiagonal(C);
      }
      else
      {
	 if (Verbose)
	    std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

	 // Non-zero diagonal element.
	 // In this case we have to consider a possible component with eigenvalue magnitude 1.
	 // We call this component (if it exists) the unit matrix.

         KComplexPolyType EParallel;  // components parallel to the identity at some momentum, may be zero

	 // Multiplication factor and the left and right eigenvalues.
	 std::complex<double> Factor = Classification.factor();
	 MatrixOperator UnitMatrixLeft = Identity;
	 MatrixOperator UnitMatrixRight = Rho;

	 // If the diagonal operator is unitary, it might have an eigenvalue of magnitude 1
	 if (Classification.is_unitary() && !Classification.is_complex_identity())
	 {
	    if (Verbose)
	       std::cerr << "Non-trivial unitary at column " << Col << std::endl;
	    
	    // Find the largest eigenvalue
	    // We need initial guess vectors in the correct symmetry sector
	    UnitMatrixLeft = MakeRandomMatrixOperator(Identity.Basis1(), Identity.Basis2(), Diag.Basis2()[0]);

	    std::complex<double> EtaL = FindLargestEigenvalue(UnitMatrixLeft, 
							      ApplyLeftQShift(Diag, Psi, QShift), Verbose);

	    if (Verbose)
	       std::cerr << "Eigenvalue of unitary operator is " << EtaL << std::endl;

	    if (std::abs(norm_frob(EtaL) - 1.0) < 1E-12)
	    {

	       UnitMatrixRight = UnitMatrixLeft;
	       // we have an eigenvalue of magnitude 1.  Find the right eigenvalue too
	       std::complex<double> EtaR = FindLargestEigenvalue(UnitMatrixRight, 
								 ApplyRightQShift(Diag, Psi, QShift), Verbose);
	       CHECK(norm_frob(EtaL-EtaR) < 1E-12)("Left and right eigenvalues do not agree!")(EtaL)(EtaR);
	       Factor = EtaL / norm_frob(EtaL);
	    }
	 }

	 // If we have an eigenvalue with magnitude 1, then decompose C into parallel and perpendicular parts
	 bool Magnitude1 = false;
         if (std::abs(norm_frob(Factor) - 1.0) < 1E-12)
         {
	    Magnitude1 = true;
	    EParallel = DecomposeParallelParts(C, Factor, UnitMatrixLeft, UnitMatrixRight);
	 }

         // Now the remaining components, which is anything that is not proportional
	 // to an eigenvector of magnitude 1
      
         KMatrixPolyType E = DecomposePerpendicularParts(C, Diag, UnitMatrixLeft, UnitMatrixRight, 
							 Psi, QShift, Magnitude1);

         // Reinsert the components parallel to the unit matrix (if any)
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
	    for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
	    {
	       E[I->first][J->first] += J->second * Identity;
	    }
         }

	 // Finally, set the E matrix element at this column
	 DEBUG_TRACE(E[1.0]);
         EMatK[Col] = E;
      }

      --Col;
   }

   // The return value is column 0 of our E matrices
   return EMatK[0];
}

void remove_redundant(OperatorComponent& Op);

Polynomial<std::complex<double> >
ExtractOverlap(Polynomial<MatrixOperator> const& E, MatrixOperator const& Rho)
{
   Polynomial<std::complex<double> > Overlap;
   for (Polynomial<MatrixOperator>::const_iterator I = E.begin(); I != E.end(); ++I)
   {
      Overlap[I->first] = inner_prod(I->second, Rho);
   }
   return Overlap;
}

int main(int argc, char** argv)
{
   if (argc < 2 || argc > 3)
   {
      std::cout << "usage: mp-iexpectation-tri <psi> [divisions]\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(argv[1], CacheSize, true);
   InfiniteWavefunction Psi = *PsiPtr;

   int Divisions = 1000;
   if (argc > 2)
      Divisions = boost::lexical_cast<int>(argv[2]);

   // Set up the operators
   int UnitCellSize = PsiPtr->size();
   SiteBlock Site = CreateU1U1HubbardSite();
   std::vector<BasisList> Sites(UnitCellSize, Site["I"].Basis().Basis());
   std::vector<SimpleOperator> FermionString(UnitCellSize, SimpleOperator(Site["P"]));
   

   // Shift Psi to the symmetric orthogonalization
   MatrixOperator Identity = MatrixOperator::make_identity(Psi.Psi.Basis1());
   MatrixOperator Rho = scalar_prod(Psi.C_right, herm(Psi.C_right));
   LinearWavefunction Phi = Psi.Psi; // no need to bugger around with C_old,C_right
 
   MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
   MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);

   bool Verbose = false;

   Phi.set_front(prod(LambdaInvSqrt, Phi.get_front()));
   Phi.set_back(prod(Phi.get_back(), delta_shift(LambdaSqrt, adjoint(Psi.QShift))));
   Rho = Psi.C_old;
   Identity = Rho;

   //   TRACE(norm_frob(operator_prod(herm(Phi), Identity, Phi) - Identity));
   //   TRACE(norm_frob(operator_prod(Phi, Rho, herm(Phi)) - Rho));

   std::cout.precision(14);

   // print the title line
   std::cout << "#k ";
   for (int i = 0; i < UnitCellSize; ++i)
   {
      for (int j = 0; j < UnitCellSize; ++j)
      {
	 std::cout << "#" << i << "," << j << "-real ";
	 std::cout << "#" << i << "," << j << "-imag ";
      }
   }
   std::cout << std::endl;

   // Loop over momenta
   double Delta = 1.0 / Divisions;
   for (double k = 0; k < 1+Delta/2; k += Delta)
   {
      std::cout << k << ' ';
      for (int i = 0; i < UnitCellSize; ++i)
      {
	 for (int j = 0; j < UnitCellSize; ++j)
	 {
	    TriangularOperator Op1 = OnePointStringOperator(Sites, FermionString, i, Site["CHup"], k*math_const::pi);
	    TriangularOperator Op2 = OnePointStringOperator(Sites, FermionString, j, Site["Cup"], -k*math_const::pi);

	    TriangularOperator Op = Op1 * Op2;

	    std::map<std::complex<double>, Polynomial<MatrixOperator>, CompareComplex> 
	       E = SolveMPO_Left(Phi, Psi.QShift, Op, Rho, Identity, Verbose);
	    Polynomial<std::complex<double> > aNorm = ExtractOverlap(E[1.0], Rho);

	    if (Verbose)
	    {
	       for (int i = 0; i <= aNorm.degree(); ++i)
	       {
		  std::cout << i << ' ' << aNorm[i].real() << ' ' << aNorm[i].imag() << '\n';
	       }
	    }
	    else
	       std::cout << aNorm[1].real() << ' ' << aNorm[1].imag() << ' ';
	 }
      }
      std::cout << std::endl;
   }

   pheap::Shutdown();
}
