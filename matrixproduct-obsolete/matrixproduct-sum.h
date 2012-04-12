// -*- C++ -*- $Id$

/*
  sum.h

  Matrix product sum function, using an improved algorithm.

  In the usual formula for the sum, we construct the matrices
  R^s = B^s \otimes B^s \oplus C^s \oplus ...
  This requires O(N^3 m^3) time and O(N^2 m^2) space to construct the sum of
  N matrix product states, each of m states, and an additional sweep to
  orthogonalize the basis.

  Alternatively, we can construct the final density matrix directly:
  psi = psi_A + psi_B + psi_C, construct the density matrix and then project onto
  the target Hilbert space:
  rho_L = P_A psi_A psi_A^\dagger P_A^\dagger
        + P_A psi_A P_{AB} psi^B\dagger P_B^\dagger
        + ....

  This takes O(N^2 m^2) space, but only O(N^2 m^3) time, so we win O(N) speed
  (plus we only need one sweep).

  The truncation error is measured in absolute units.
*/

#if !defined(MATRIXPRODUCT_SUM_H_hdsjhiu438y89y8hf98ry8y437fh)
#define MATRIXPRODUCT_SUM_H_hdsjhiu438y89y8hf98ry8y437fh

#include "centerwavefunction.h"
#include "handlestack.h"

inline
CenterWavefunction
mpwavefunction_sum(std::vector<CenterWavefunction> X, StatesInfo const& SInfo, bool ShowStates = false)
{
   int const n = X.size();

   typedef MPStateComponent Component;
   typedef Component::OperatorType OperatorType;

   QuantumNumber Ident(X[0].GetSymmetryList());  // the scalar quantum number

   // The lower-triangular matrix of Transform operators.
   typedef LinearAlgebra::Matrix<OperatorType> InternalTransType;
   HandleStack<InternalTransType> LeftInternal;
   InternalTransType Top(n,n);
   OperatorType LeftVac = OperatorType::make_identity(X[0].LeftVacuumBasis());
   for (int i = 0; i < n; ++i)
   {
      for (int j = 0; j < i; ++j)
      {
         Top(i,j) = LeftVac;
      }
   }
   LeftInternal.push(Top);

   // We start from the right hand side here, build the internal transform stack at the same time
   for (int i = 0; i < n; ++i)
   {
      CHECK_EQUAL(X[i].LeftSize(), 1);
   }
   while (X[0].RightSize() > 1)
   {
      for (int i = 0; i < n; ++i)
      {
         for (int j = 0; j < i; ++j)
         {
            Top(i,j) = operator_prod(herm(X[i].Left()), Top(i,j), X[j].Left());
         }
      }
      LeftInternal.push(Top);
      for (int i = 0; i < n; ++i)
      {
         X[i].RotateRight();
      }
   }
   for (int i = 0; i < n; ++i)
   {
      for (int j = 0; j < i; ++j)
      {
         Top(i,j) = operator_prod(herm(X[i].Left()), Top(i,j), X[j].Left());
      }
   }
   LeftInternal.push(Top);
   
   CenterWavefunction Result;
   Result.PushRight(Component::ConstructFullBasis1(X[0].Right().SiteBasis(), 
						   X[0].RightVacuumBasis()));

   LinearAlgebra::Vector<OperatorType> RightMap(n);        // mapping of state i to the result

   OperatorType RightVac = OperatorType::make_identity(X[0].RightVacuumBasis());
   for (int i = 0; i < n; ++i)
   {
      RightMap[i] = RightVac;
   }
   for (int i = 0; i < n; ++i)
   {
      RightMap[i] = operator_prod(Result.Right(), RightMap[i], herm(X[i].Right()));
   }

   while (X[0].LeftSize() > 1)
   {
      // construct right density matrix, a short cut is to construct only half the off-diagonal
      // parts then add the hermitian conjugate
      OperatorType Rho;

      for (int i = 1; i < n; ++i)
      {
	 // (j,i)
	 for (int j = 0; j < i; ++j)
	 {
	    Rho += triple_prod(RightMap[i], 
			       triple_prod(herm(X[i].Center()), LeftInternal.top()(i,j), X[j].Center()), 
			       herm(RightMap[j]));
	 }
      }
      Rho += adjoint(Rho);

      // diagonal parts
      for (int i = 0; i < n; ++i)
      {
	 Rho += triple_prod(RightMap[i], 
			    scalar_prod(herm(X[i].Center()), X[i].Center()), 
			    herm(RightMap[i]));
      }

      // Construct the truncator for the result
      DensityMatrix<OperatorType> DM(Rho);
      TruncationInfo Info;
      MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
                                                                                              SInfo,
                                                                                              Info));
      if (ShowStates)
         std::cout << "forming sum, states=" << Info.KeptStates() << ", trunc=" << Info.TruncationError() << '\n';
      //DM.DensityMatrixReport(std::cout);
      DEBUG_TRACE(Info.KeptStates());
      Result.Right() = prod(U, Result.Right());

      // truncate the transformers, and rotate the original states
      for (int i = 0; i < n; ++i)
      {
	 RightMap[i] = prod(U, RightMap[i], Ident);
	 X[i].RotateLeft();
      }
      LeftInternal.pop();

      // next right operator
      Result.PushRight(Component::ConstructFullBasis1(X[0].Right().SiteBasis(), Result.Right().Basis1()));

      for (int i = 0; i < n; ++i)
      {
	 RightMap[i] = operator_prod(Result.Right(), RightMap[i], herm(X[i].Right()));
      }
   }

   // Now approach from the left hand side
   LinearAlgebra::Vector<OperatorType> LeftMap(n);        // mapping of state i to the result

   Result.PushLeft(Component::ConstructFullBasis2(X[0].LeftVacuumBasis(),
						  X[0].Left().SiteBasis()));

   //   OperatorType LeftVac = OperatorType::make_identity(X[0].LeftVacuumBasis());
   for (int i = 0; i < n; ++i)
   {
      LeftMap[i] = scalar_prod(herm(Result.Left()), X[i].Left());
   }

   // Construct the center matrix
   Result.Center() = triple_prod(LeftMap[0], X[0].Center(), herm(RightMap[0]));
   for (int i = 1; i < n; ++i)
   {
      Result.Center() += triple_prod(LeftMap[i], X[i].Center(), herm(RightMap[i]));
   }

   return Result;
}

#endif
