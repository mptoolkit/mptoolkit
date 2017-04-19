// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-idmrg-s3e.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// iDMRG with single-site subspace expansion
//
// A note on initialization.
// The wavefunction is stored in the left-canonical form, of A-matrices.
// We also have the lambdaR matrix on the right hand side.
// For determining the left Hamiltonian matrix elements, lambdaR isn't used
// since the wavefunction is defined by the infinite string AAAAA....
// Thus, lambdaR is only a short-cut for (1) determining the density matrix in the left canonical basis,
// and (2) speeding up finding the right-canonical form.
//
// To find the right canonical form, we insert lambdaR and othogonalize matrices to right-orthogonalized,
// ending up with AAAAA lambdaR = U lambdaL BBBBB
// The U is the matrix left over from the SVD.  In principle, this commutes with lambdaL, but there is
// no guarantee that U.Basis1() == U.Basis2() (they may be in a different order, for example),
// and similarly in principle lambdaL == lambdaR, but there may be numerical differences, ordering of
// singular values etc.
//
// So, to find the right Hamiltonian matrix elements, we shift the B matrix to the right hand side.
// Now, we have the identity
// U lambdaL BBBBB = AAAAA lambdaR
// so upon shifting U to the other side, we are changing to the basis
// lambdaL (U^\dagger U) BBBBB U = U^\dagger AAAAA (U U^\dagger) lambdaR U
//
// Hence we can also interpret U as the unitary that maps the basis of lambdaR into the basis of lambdaL.
// Having obtained the matrix elements of the right Hamiltonian in this basis (BBBBB U), we could either shift
// back to the original basis, or shift the left Hamiltonian to this basis.
// To shift the right Hamiltonian back to the original basis, we need to convert back to the
// BBBBB basis, which is effected by U BBBBB U^\dagger, so we need to act on the right Hamiltonian with
// this operator, H -> U * H * U^\dagger
//
// Alternatively, we could shift the left Hamiltonian to the new basis.  To do this, we write
// U^\dagger AAAAA U U^\dagger lambdaR U = lambdaL BBBBB U
//
// And we see that we transform AAAAA -> U^\dagger AAAAA U
// lambdaR -> U^\dagger lambdaR U
// and the Hamiltonian matrix elements will change by H -> U^\dagger H U

#include "mpo/triangular_mpo.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include "common/prog_options.h"
#include <iostream>
#include "common/environment.h"
#include "common/statistics.h"
#include "common/prog_opt_accum.h"
#include "mp-algorithms/gmres.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mp-algorithms/stateslist.h"
#include "wavefunction/operator_actions.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

#if !defined(NDEBUG)
#include "mp-algorithms/triangular_mpo_solver.h"
#endif

#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "mp-algorithms/eigensolver.h"

namespace prog_opt = boost::program_options;

using statistics::moving_average;
using statistics::moving_exponential;

double const iTol = 1E-7;

MatrixOperator GlobalU;

bool EarlyTermination = false;  // we set this to true if we get a checkpoint


LinearAlgebra::DiagonalMatrix<double>
InvertDiagonal(LinearAlgebra::DiagonalMatrix<double> const& D, double Tol = 1E-15)
{
   LinearAlgebra::DiagonalMatrix<double> Result(size1(D), size2(D));
   for (unsigned i = 0; i < size1(D); ++i)
   {
      Result.diagonal()[i] = norm_frob(D.diagonal()[i]) < Tol ? 0.0 : 1.0 / D.diagonal()[i];
   }
   return Result;
}

double const Alpha = 100;

// function to calculate
// (D1 * U) * Inverse(D2)
// where D1 and D2 are diagonal, and U is unitary.
// U may be non-square, in which case it is only row or column unitary.
// This depends on whether we are increasing or reducing the number of states.
#if 0
MatrixOperator
Solve_DU_DInv(MatrixOperator const& DU, RealDiagonalOperator const& D)
{
   int Dim1 = DU.Basis1().total_dimension();
   int Dim2 = DU.Basis2().total_dimension();
   MatrixOperator Result = DU * InvertDiagonal(D,1E-8);
   TRACE(norm_frob_sq(Result)/(Dim1*Dim2));
   return Result;
}
#else
MatrixOperator
Solve_DU_DInv(MatrixOperator const& DU, RealDiagonalOperator const& D)
{
   DU.check_structure();

   MatrixOperator Result(DU.Basis1(), D.Basis1());
   for (MatrixOperator::const_iterator I = iterate(DU); I; ++I)
   {
      for (MatrixOperator::const_inner_iterator J = iterate(I); J; ++J)
      {
         LinearAlgebra::Matrix<std::complex<double> > Component = (*J) * InvertDiagonal(D(J.index2(),J.index2()),iTol);
         Result(J.index1(), J.index2()) = Component;
         //TRACE(norm_frob_sq(Component) / (size1(Component)*size2(Component)));
         //TRACE(Component);
      }
   }
   return Result;
}
#endif

MatrixOperator
Solve_DInv_UD(RealDiagonalOperator const& D, MatrixOperator const& UD)
{
   int Dim1 = UD.Basis1().total_dimension();
   int Dim2 = UD.Basis2().total_dimension();
   MatrixOperator Result = InvertDiagonal(D,iTol) * UD;
   TRACE(norm_frob_sq(Result)/sqrt(Dim1*Dim2));
   return Result;
}


// solves D*U*InvertDiagonal(E)
#if 0
MatrixOperator
Solve_D_U_DInv(RealDiagonalOperator const& D, MatrixOperator const& U, RealDiagonalOperator const& E)
{
   return D*U*InvertDiagonal(E,iTol);
}

MatrixOperator
Solve_DInv_U_D(RealDiagonalOperator const& D, MatrixOperator const& U, RealDiagonalOperator const& E)
{
   return InvertDiagonal(D,iTol)*U*E;
}

#else
MatrixOperator
Solve_D_U_DInv(RealDiagonalOperator const& D, MatrixOperator const& U, RealDiagonalOperator const& E)
{
   MatrixOperator Result(D.Basis1(), E.Basis1());
   for (MatrixOperator::const_iterator I = iterate(U); I; ++I)
   {
      for (MatrixOperator::const_inner_iterator J = iterate(I); J; ++J)
      {
         int Dim1 = D.Basis1().dim(J.index1());
         int Dim2 = E.Basis1().dim(J.index2());

         LinearAlgebra::DiagonalMatrix<double>
            DD = D(J.index1(), J.index1()) + LinearAlgebra::DiagonalMatrix<double>(Dim1, Dim1, Alpha);
         LinearAlgebra::DiagonalMatrix<double>
            EE = E(J.index2(), J.index2()) + LinearAlgebra::DiagonalMatrix<double>(Dim2, Dim2, Alpha);

         LinearAlgebra::Matrix<std::complex<double>> Component =
            DD * (*J) * InvertDiagonal(EE, iTol);
         Result(J.index1(), J.index2()) = Component;
         //TRACE(norm_frob_sq(Component) / (size1(Component)*size2(Component)));
         //TRACE(Component);
      }
   }
   return Result;
}

MatrixOperator
Solve_DInv_U_D(RealDiagonalOperator const& D, MatrixOperator const& U, RealDiagonalOperator const& E)
{
   MatrixOperator Result(D.Basis1(), E.Basis1());
   for (MatrixOperator::const_iterator I = iterate(U); I; ++I)
   {
      for (MatrixOperator::const_inner_iterator J = iterate(I); J; ++J)
      {
         int Dim1 = D.Basis1().dim(J.index1());
         int Dim2 = E.Basis1().dim(J.index2());

         LinearAlgebra::DiagonalMatrix<double>
            DD = D(J.index1(), J.index1()) + LinearAlgebra::DiagonalMatrix<double>(Dim1, Dim1, Alpha);
         LinearAlgebra::DiagonalMatrix<double>
            EE = E(J.index2(), J.index2()) + LinearAlgebra::DiagonalMatrix<double>(Dim2, Dim2, Alpha);

         LinearAlgebra::Matrix<std::complex<double>> Component =
            InvertDiagonal(DD, iTol) * (*J) * EE;
         Result(J.index1(), J.index2()) = Component;
         //TRACE(norm_frob_sq(Component) / (size1(Component)*size2(Component)));
         //TRACE(Component);
      }
   }
   return Result;
}



#endif




struct ProductLeft
{
   typedef StateComponent result_type;
   typedef StateComponent argument_type;

   ProductLeft(LinearWavefunction const& Psi_, TriangularMPO const& Op_,
               QuantumNumber const& QShift_)
      : Psi(Psi_), Op(Op_), QShift(QShift_)
   {
   }

   StateComponent operator()(StateComponent const& In) const
   {
      StateComponent Guess = delta_shift(In, QShift);
      TriangularMPO::const_iterator HI = Op.begin();
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I, ++HI)
      {
         Guess = contract_from_left(*HI, herm(*I), Guess, *I);
      }
      DEBUG_CHECK(HI == Op.end());
      return Guess; //delta_shift(Guess, adjoint(QShift));
   }

   LinearWavefunction const& Psi;
   TriangularMPO const& Op;
   QuantumNumber QShift;
};

struct ProductRight
{
   typedef StateComponent result_type;
   typedef StateComponent argument_type;

   ProductRight(LinearWavefunction const& Psi_, TriangularMPO const& Op_,
                QuantumNumber const& QShift_)
      : Psi(Psi_), Op(Op_), QShift(QShift_)
   {
   }

   StateComponent operator()(StateComponent const& In) const
   {
      StateComponent Guess = In;
      LinearWavefunction::const_iterator I = Psi.end();
      TriangularMPO::const_iterator HI = Op.end();
      while (I != Psi.begin())
      {
         --I;
         --HI;
         Guess = contract_from_right(herm(*HI), *I, Guess, herm(*I));
         //Guess = operator_prod(conj(*HI), *I, Guess, herm(*I));
      }
      DEBUG_CHECK(HI == Op.begin());
      return delta_shift(Guess, adjoint(QShift));
   }

   LinearWavefunction const& Psi;
   TriangularMPO const& Op;
   QuantumNumber QShift;
};

// Does Iter iterations of the MPO solver, without assuming that Psi
// is perfectly orthogonal - that is, Rho and Identity are approximations.
// The number of GMRES iterations is fixed by Iter, which is typically rather small.
// On input, E is an approximation for the final solution,

std::complex<double>
PartialSolveSimpleMPO_Left(StateComponent& E, StateComponent const& OldE,
                           LinearWavefunction const& Psi,
                           QuantumNumber const& QShift, TriangularMPO const& Op,
                           MatrixOperator const& Rho, int Iter)
{
}

struct MPSMultiply
{
   MPSMultiply(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_)
      : E(E_), H(H_), F(F_)
   {
   }

   StateComponent operator()(StateComponent const& Psi) const
   {
      StateComponent R = operator_prod_inner(H, E, Psi, herm(F));
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
};

struct MixInfo
{
   double MixFactor;
   double RandomMixFactor;
};

//#define SSC

// Apply subspace expansion / truncation on the left (C.Basis1()).
// Returns a matrix Lambda (diagonal) and a unitary
// Postcondition: U' Lambda' C' = C (up to truncation!)
std::pair<MatrixOperator, RealDiagonalOperator>
SubspaceExpandBasis1(StateComponent& C, OperatorComponent const& H, StateComponent const& RightHam,
                     MixInfo const& Mix, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& LeftHam)
{
   // truncate - FIXME: this is the s3e step
#if defined(SSC)
   MatrixOperator Lambda;
   SimpleStateComponent CX;
   std::tie(Lambda, CX) = ExpandBasis1_(C);
#else
   MatrixOperator Lambda = ExpandBasis1(C);
#endif

   MatrixOperator Rho = scalar_prod(herm(Lambda), Lambda);
   if (Mix.MixFactor > 0)
   {
#if defined(SSC)
      StateComponent RH = contract_from_right(herm(H), CX, RightHam, herm(CX));
#else
      StateComponent RH = contract_from_right(herm(H), C, RightHam, herm(C));
#endif
      MatrixOperator RhoMix;
      MatrixOperator RhoL = scalar_prod(Lambda, herm(Lambda));

      // Skip the identity and the Hamiltonian
      for (unsigned i = 1; i < RH.size()-1; ++i)
      {
         double Prefactor = trace(triple_prod(herm(LeftHam[i]), RhoL, LeftHam[i])).real();
         if (Prefactor == 0)
            Prefactor = 1;
         RhoMix += Prefactor * triple_prod(herm(RH[i]), Rho, RH[i]);
      }
      // check for a zero mixing term - can happen if there are no interactions that span
      // the current bond
      double RhoTrace = trace(RhoMix).real();
      if (RhoTrace != 0)
         Rho += (Mix.MixFactor / RhoTrace) * RhoMix;
   }
   if (Mix.RandomMixFactor > 0)
   {
      MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
      RhoMix = herm(RhoMix) * RhoMix;
      Rho += (Mix.RandomMixFactor / trace(RhoMix)) * RhoMix;
   }

   //TRACE(Rho);

   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);
   Lambda = Lambda * herm(U);

   //TRACE(Lambda);

#if defined(SSC)
   C = U*CX; //prod(U, CX);
#else
   C = prod(U, C);
#endif

   MatrixOperator Vh;
   RealDiagonalOperator D;
   SingularValueDecompositionKeepBasis2(Lambda, U, D, Vh);

   //TRACE(U)(D)(Vh);

   C = prod(Vh, C);
   return std::make_pair(U, D);
}

// Apply subspace expansion / truncation on the right (C.Basis2()).
// Returns Lambda matrix (diagonal) and a unitary matrix
// Postcondition: C' Lambda' U' = C (up to truncation!)
std::pair<RealDiagonalOperator, MatrixOperator>
SubspaceExpandBasis2(StateComponent& C, OperatorComponent const& H, StateComponent const& LeftHam,
                     MixInfo const& Mix, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& RightHam)
{
   // truncate - FIXME: this is the s3e step
   MatrixOperator Lambda = ExpandBasis2(C);

   MatrixOperator Rho = scalar_prod(Lambda, herm(Lambda));
   if (Mix.MixFactor > 0)
   {
      StateComponent LH = contract_from_left(H, herm(C), LeftHam, C);
      MatrixOperator RhoMix;

      MatrixOperator RhoR = scalar_prod(herm(Lambda), Lambda);

      for (unsigned i = 1; i < LH.size()-1; ++i)
      {
         double Prefactor = trace(triple_prod(herm(RightHam[i]), RhoR, RightHam[i])).real();
         if (Prefactor == 0)
            Prefactor = 1;
         RhoMix += Prefactor * triple_prod(LH[i], Rho, herm(LH[i]));
      }
      double RhoTrace = trace(RhoMix).real();
      if (RhoTrace != 0)
         Rho += (Mix.MixFactor / RhoTrace) * RhoMix;
   }
   if (Mix.RandomMixFactor > 0)
   {
      MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
      RhoMix = herm(RhoMix) * RhoMix;
      Rho += (Mix.RandomMixFactor / trace(RhoMix)) * RhoMix;
   }
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);

   Lambda = U * Lambda;
   C = prod(C, herm(U));

   MatrixOperator Vh;
   RealDiagonalOperator D;
   SingularValueDecompositionKeepBasis1(Lambda, U, D, Vh);

   C = prod(C, U);

   return std::make_pair(D, Vh);
}

class iDMRG
{
   public:
      // Construct an iDMRG object.  It is assumed that Psi_ is in left-canonical form, with
      // LambdaR being the lambda matrix on the right edge.
      iDMRG(LinearWavefunction const& Psi_, RealDiagonalOperator const& LambdaR,
            MatrixOperator const& UR,
            QuantumNumber const& QShift_, TriangularMPO const& Hamiltonian_,
            StateComponent const& LeftHam, StateComponent const& RightHam,
            std::complex<double> InitialEnergy = 0.0, int Verbose = 0);

      void SetMixInfo(MixInfo const& m);

      void SetInitialFidelity(double f)
      {
         Solver_.SetInitialFidelity(Psi.size(), f);
      }

      LocalEigensolver& Solver() { return Solver_; }

      void SweepRight(StatesInfo const& SInfo, double HMix = 0);
      void SweepLeft(StatesInfo const& SInfo, double HMix = 0, bool NoUpdate = false);

      // call after a SweepRight() to make Psi an infinite wavefunction
      void Finish(StatesInfo const& SInfo);

      // returns the wavefunction, which is in center-regular form
      // (all MPS orthogonalized except for one site, and Psi.Basis1() == Psi.Basis2())
      LinearWavefunction const& Wavefunction() const { return Psi; }

      // Above functions are implemented in terms of:

      // construct initial hamiltonian, called by the constructor
      void Initialize(RealDiagonalOperator const& LambdaR, MatrixOperator const& UR,
                      std::complex<double> InitialEnergy);

      void Solve();

      // at the left boundary of the unit cell, construct and save the RightHamiltonian for later use
      void SaveRightBlock(StatesInfo const& States);

      // at the right boundary of the unit cell, construct and save the LeftHamiltonian for later use
      void SaveLeftBlock(StatesInfo const& States);

      // at the left boundary of the unit cell, update the LeftHamiltonian and C to SaveLeftHamiltonian
      // At the end of this step, the wavefunction is in a 'regular' form with
      // C.Basis1() == Psi.Basis2()
      void UpdateLeftBlock(double HMix = 0);

      void UpdateRightBlock(double HMix = 0);

      void TruncateAndShiftLeft(StatesInfo const& States);
      void TruncateAndShiftRight(StatesInfo const& States);

      void ShowInfo(char c);

      void CheckConsistency() const;

      //void SubspaceExpandLeft();
      //void SubspaceExpandRight();


      //   private:
      TriangularMPO Hamiltonian;
      LinearWavefunction Psi;
      QuantumNumber QShift;

      std::deque<StateComponent> LeftHamiltonian;
      std::deque<StateComponent> RightHamiltonian;

      LinearWavefunction::iterator C;
      TriangularMPO::const_iterator H;

      int Verbose;

      // iterators pointing to the edges of the unit cell.
      // FirstSite = Psi.begin()
      // LastSite = Psi.end() - 1
      LinearWavefunction::const_iterator FirstSite, LastSite;

      LocalEigensolver Solver_;

      StateComponent SaveLeftHamiltonian;
      RealDiagonalOperator SaveLambda2;
      MatrixOperator SaveU2;

      StateComponent SaveRightHamiltonian;
      MatrixOperator SaveU1;
      RealDiagonalOperator SaveLambda1;

      MixInfo    MixingInfo;
      TruncationInfo Info;

      int SweepNumber;
};

iDMRG::iDMRG(LinearWavefunction const& Psi_, RealDiagonalOperator const& LambdaR, MatrixOperator const& UR,
             QuantumNumber const& QShift_, TriangularMPO const& Hamiltonian_,
             StateComponent const& LeftHam, StateComponent const& RightHam,
             std::complex<double> InitialEnergy, int Verbose_)
   : Hamiltonian(Hamiltonian_), Psi(Psi_), QShift(QShift_),
     LeftHamiltonian(1, LeftHam), RightHamiltonian(1, RightHam), Verbose(Verbose_), SweepNumber(1)
{
   this->Initialize(LambdaR, UR, InitialEnergy);
}

void
iDMRG::Initialize(RealDiagonalOperator const& LambdaR, MatrixOperator const& UR, std::complex<double> InitialEnergy)
{
   // the starting point is at the right hand side, so fill out the LeftHamiltonian
   // fill out the LeftHamiltonian
   C = Psi.begin();
   H = Hamiltonian.begin();

   FirstSite = Psi.begin();
   LastSite = Psi.end();
   --LastSite;

   CHECK_EQUAL(LeftHamiltonian.back().Basis2(), C->Basis1());

   // Generate the left Hamiltonian matrices at each site, up to (but not including) the last site.
   while (C != LastSite)
   {
      LeftHamiltonian.push_back(contract_from_left(*H, herm(*C), LeftHamiltonian.back(), *C));
      ++C;
      ++H;
   }

   // make *C the 'centre MPS' by incorporating LambdaR into it
   *C = prod(*C, LambdaR);

   // The InitialEnergy is the energy per unit cell.  At this point,
   // the Hamiltonian we have represents the sum of interactions across
   // one unit cell, but the boundary interactions are counted at both boundaries.
   // So this gives an energy that is (UnitCellSize+1) times the energy per site.
   // To compensate, we subtract a term off the RightHamiltonian
   // so that the initial energy is correct.

   if (InitialEnergy != 0.0)
   {
      if (Verbose > 0)
         std::cerr << "Correcting initial energy to " << InitialEnergy << '\n';

      StateComponent R = operator_prod_inner(*H, LeftHamiltonian.back(), *C, herm(RightHamiltonian.front()));
      std::complex<double> E = inner_prod(*C, R) / inner_prod(*C, *C);

      if (Verbose > 0)
      {
         std::cerr << "Raw initial energy was " << E << '\n';
         if (Verbose > 1)
            std::cerr << "Initial wavefunction norm is " << norm_frob(*C) << '\n';
         StateComponent Resid = R - E*(*C);
         std::cerr << "Initial residual norm is " << norm_frob(Resid) << '\n';
      }

      RightHamiltonian.front().front() += (InitialEnergy - E) * RightHamiltonian.front().back();
   }

   // We are going to start by sweeping left, so we need to initialize the
   // left blocks such that UpdateLeftBlock() will function
   SaveLeftHamiltonian = LeftHamiltonian.front();
   SaveRightHamiltonian = RightHamiltonian.front();
   SaveLambda2 = LambdaR; //MatrixOperator::make_identity(Psi.Basis1());
   SaveU2 = UR;

   SaveU1 = adjoint(UR);
   SaveLambda1 = delta_shift(LambdaR, QShift); // MatrixOperator::make_identity(Psi.Basis2());

   this->CheckConsistency();
}

//
// we need to store SaveLamda1 as a combination of a RealDiagonalOperator and
// a unitary MatrixOperator.
//

void
iDMRG::UpdateLeftBlock(double HMix)
{
   if (Verbose > 2)
   {
      std::cerr << "Updating left Hamiltonian.  Old basis has "
                << LeftHamiltonian.back().Basis1().total_dimension()
                << " states, new basis has states = "
                << SaveLeftHamiltonian.Basis1().total_dimension() << "\n";
   }

   StateComponent E = LeftHamiltonian.back();
   LeftHamiltonian = std::deque<StateComponent>(1, delta_shift(SaveLeftHamiltonian, QShift));

   MatrixOperator T = Solve_D_U_DInv(delta_shift(SaveLambda2, QShift), delta_shift(SaveU2, QShift), SaveLambda1);

   T = T * herm(SaveU1);

   (*C) = prod(T, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   //HMix = 0;
   if (HMix != 0)
   {
      // Adjust basis of E
      MatrixOperator V = delta_shift(SaveU2, QShift) * herm(SaveU1);
      E = triple_prod(V, E, herm(V));
      LeftHamiltonian.back() = HMix * E + (1.0 - HMix) * LeftHamiltonian.back();
   }

   // Subtract off the energy
   LeftHamiltonian.back().back() -= Solver_.LastEnergy() * LeftHamiltonian.back().front();

   this->CheckConsistency();
}

void
iDMRG::UpdateRightBlock(double HMix)
{
   if (Verbose > 2)
   {
      std::cerr << "Updating right Hamiltonian.  Old basis has "
                << RightHamiltonian.front().Basis1().total_dimension()
                << " states, new basis has states = "
                << SaveRightHamiltonian.Basis1().total_dimension() << "\n";
   }

   StateComponent F = RightHamiltonian.front();
   RightHamiltonian = std::deque<StateComponent>(1, delta_shift(SaveRightHamiltonian, adjoint(QShift)));

   //TRACE(SaveLambda2)(SaveU1)(SaveLambda1);
   MatrixOperator T = Solve_DInv_U_D(SaveLambda2, delta_shift(SaveU1, adjoint(QShift)),
                                     delta_shift(SaveLambda1, adjoint(QShift)));

   T = herm(SaveU2) * T;
   //TRACE(T);

   (*C) = prod(*C, T);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   //HMix = 0;
   if (HMix != 0)
   {
      // Adjust basis of F
      MatrixOperator V = herm(SaveU2) * delta_shift(SaveU1, adjoint(QShift));
      V = T;
      F = triple_prod(herm(V), F, V);

      for (unsigned n = 0; n < F.size(); ++n)
      {
         TRACE(norm_frob(F[n] - RightHamiltonian.front()[n]));
      }

      RightHamiltonian.front() = HMix * F + (1.0 - HMix) * RightHamiltonian.front();
   }

   // Subtract off the energy
   RightHamiltonian.front().front() -= Solver_.LastEnergy() * RightHamiltonian.front().back();

   this->CheckConsistency();
}

void
iDMRG::SaveLeftBlock(StatesInfo const& States)
{
   // When we save the block, we need to end up with
   // C.Basis2() == SaveLeftHamiltonian.Basis()
   CHECK(C == LastSite);
   StateComponent L = *C;
   std::tie(SaveLambda2, SaveU2) = SubspaceExpandBasis2(L, *H, LeftHamiltonian.back(),
                                                          MixingInfo, States, Info, RightHamiltonian.front());

   if (Verbose > 1)
   {
      std::cerr << "Saving left block for idmrg, states=" << Info.KeptStates()
                << " L.Basis2() is " << L.Basis2().total_dimension() << '\n';
   }
   SaveLeftHamiltonian = contract_from_left(*H, herm(L), LeftHamiltonian.back(), L);
   this->CheckConsistency();
}

void
iDMRG::SaveRightBlock(StatesInfo const& States)
{
   CHECK(C == FirstSite);
   StateComponent R = *C;
   std::tie(SaveU1, SaveLambda1) = SubspaceExpandBasis1(R, *H, RightHamiltonian.front(),
                                                          MixingInfo, States, Info, LeftHamiltonian.back());

   //   TRACE(SaveU1);

   if (Verbose > 1)
   {
      std::cerr << "Saving right block for idmrg, states=" << Info.KeptStates() << '\n';
   }
   SaveRightHamiltonian = contract_from_right(herm(*H), R, RightHamiltonian.front(), herm(R));
   this->CheckConsistency();
}

void
iDMRG::TruncateAndShiftLeft(StatesInfo const& States)
{
   this->CheckConsistency();
   // Truncate right
   MatrixOperator U;
   RealDiagonalOperator Lambda;
   std::tie(U, Lambda) = SubspaceExpandBasis1(*C, *H, RightHamiltonian.front(), MixingInfo, States, Info,
                                                LeftHamiltonian.back());

   if (Verbose > 1)
   {
      std::cerr << "Truncating left basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks
   RightHamiltonian.push_front(contract_from_right(herm(*H), *C, RightHamiltonian.front(), herm(*C)));
   LeftHamiltonian.pop_back();

   // next site
   --H;
   --C;

   *C = prod(*C, U*Lambda);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   this->CheckConsistency();
}

void
iDMRG::CheckConsistency() const
{
   CHECK_EQUAL(LeftHamiltonian.back().Basis2(), C->Basis1());
   CHECK_EQUAL(RightHamiltonian.front().Basis1(), C->Basis2());
   CHECK_EQUAL(LeftHamiltonian.back().LocalBasis(), H->Basis1());
   CHECK_EQUAL(RightHamiltonian.front().LocalBasis(), H->Basis2());
   CHECK_EQUAL(H->LocalBasis2(), C->LocalBasis());
}

void
iDMRG::TruncateAndShiftRight(StatesInfo const& States)
{
   // Truncate right
   RealDiagonalOperator Lambda;
   MatrixOperator U;
   std::tie(Lambda, U) = SubspaceExpandBasis2(*C, *H, LeftHamiltonian.back(), MixingInfo, States, Info,
                                                RightHamiltonian.front());
   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks
   LeftHamiltonian.push_back(contract_from_left(*H, herm(*C), LeftHamiltonian.back(), *C));
   RightHamiltonian.pop_front();

   // next site
   ++H;
   ++C;

   *C = prod(Lambda*U, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   this->CheckConsistency();
}

void
iDMRG::SweepRight(StatesInfo const& States, double HMix)
{
   this->UpdateLeftBlock(HMix);
   this->Solve();
   this->SaveRightBlock(States);
   this->ShowInfo('P');

   while (C != LastSite)
   {
      this->TruncateAndShiftRight(States);
      this->Solve();
      this->ShowInfo('R');
   }

   SweepNumber++;
}

void
iDMRG::SweepLeft(StatesInfo const& States, double HMix, bool NoUpdate)
{
   if (!NoUpdate)
      this->UpdateRightBlock(HMix);
   this->Solve();
   this->SaveLeftBlock(States);
   this->ShowInfo('Q');

   while (C != FirstSite)
   {
      this->TruncateAndShiftLeft(States);
      this->Solve();
      this->ShowInfo('L');
   }

   SweepNumber++;
}

void
iDMRG::Finish(StatesInfo const& States)
{
   CHECK(C == LastSite);

   // The final truncation.
   // This is actually quite important to get a translationally invariant wavefunction
   RealDiagonalOperator Lambda;
   MatrixOperator U;
   std::tie(Lambda, U) = SubspaceExpandBasis2(*C, *H, LeftHamiltonian.back(), MixingInfo, States, Info,
                                                RightHamiltonian.front());

   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }

   this->ShowInfo('F');

   (*C) = prod(*C, Solve_D_U_DInv(Lambda, U*herm(SaveU2), SaveLambda2));

   CHECK_EQUAL(Psi.Basis1(), delta_shift(Psi.Basis2(), QShift));
}

void
iDMRG::Solve()
{
   Solver_.Solve(*C, LeftHamiltonian.back(), *H, RightHamiltonian.front());
}

void
iDMRG::ShowInfo(char c)
{
   std::cout << c
             << " Sweep=" << SweepNumber
             << " Energy=";
   if (Solver_.is_complex())
      std::cout << Solver_.LastEnergy();
   else
      std::cout << Solver_.LastEnergyReal();
   std::cout << " States=" << Info.KeptStates()
             << " TruncError=" << Info.TruncationError()
             << " Entropy=" << Info.KeptEntropy()
             << " Fidelity=" << Solver_.LastFidelity()
      //                         << " FidelityAv=" << FidelityAv.value()
             << " Iter=" << Solver_.LastIter()
             << " Tol=" << Solver_.LastTol()
             << '\n';
}

int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true, false);
   try
   {
      // flush cout if we write to cerr
      std::cout.tie(&std::cerr);

      int NumIter = 20;
      int MinIter = 4;
      int MinStates = 1;
      std::string States = "100";
      int NumSteps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;
      std::string FName;
      std::string HamStr;
      bool Force = false;
      bool TwoSite = true;
      bool OneSite = false;
      int WavefuncUnitCellSize = 0;
      double MixFactor = 0.02;
      double RandomMixFactor = 0.0;
      bool NoFixedPoint = false;
      bool NoOrthogonalize = false;
      bool Create = false;
      bool ExactDiag = false;
      double FidelityScale = 1.0;
      int Verbose = 0;
      bool Quiet = false;
      bool DoRandom = false; // true if we want to start an iteration from a random centre matrix
      std::string TargetState;
      std::vector<std::string> BoundaryState;
      double EvolveDelta = 0.0;
      double InitialFidelity = 1E-7;
      //double ArnoldiTol = 1E-14;
      double GMRESTol = 1E-13;    // tolerance for GMRES for the initial H matrix elements.
                                  // ** 2016-01-25: 1e-14 seems too small, failed to converge with final tol 1.5e-14, increasing to 2e-14
                                  // ** 2016-02-23: increasing again to 3e-14 on an example that fails to converge
                                  // ** 2016-03-16: increased to 1e-13.  It perhaps needs to scale with the unit cell size?
      double MaxTol = 4E-4;  // never use an eigensolver tolerance larger than this
      double MinTol = 1E-16; // lower bound for the eigensolver tolerance - seems we dont really need it
      double HMix = 0;  // Hamiltonian length-scale mixing factor
      double ShiftInvertEnergy = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian, of the form lattice:operator")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
         ("force,f", prog_opt::bool_switch(&Force), "Allow overwriting output files")
         ("two-site,2", prog_opt::bool_switch(&TwoSite), "Modify two sites at once (default)")
         ("one-site,1", prog_opt::bool_switch(&OneSite), "Modify one site at a time")
#if defined(ENABLE_ONE_SITE_SCHEME)
         ("onesiteboundary", prog_opt::bool_switch(&UseOneSiteScheme), "Modify one site also at the boundary")
#endif
         ("states,m", prog_opt::value(&States),
          FormatDefault("number of states, or a StatesList", States).c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
          FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("mix-factor", prog_opt::value(&MixFactor),
          FormatDefault("Mixing coefficient for the density matrix", MixFactor).c_str())
         ("random-mix-factor", prog_opt::value(&RandomMixFactor),
          FormatDefault("Random mixing for the density matrix", RandomMixFactor).c_str())
         ("hmix", prog_opt::value(&HMix),
          FormatDefault("Hamiltonian mixing factor", HMix).c_str())
         ("evolve", prog_opt::value(&EvolveDelta),
          "Instead of Lanczos, do imaginary time evolution with this timestep")
         ("random,a", prog_opt::bool_switch(&Create),
          "Create a new wavefunction starting from a random state")
         ("unitcell,u", prog_opt::value(&WavefuncUnitCellSize),
          "Only if --create is specified, the size of the wavefunction unit cell")
         ("startrandom", prog_opt::bool_switch(&DoRandom),
          "Start the first iDMRG iteration from a random centre matrix")
         ("exactdiag,e", prog_opt::bool_switch(&ExactDiag),
          "Start from an effective exact diagonalization of the unit cell")
         ("target,q", prog_opt::value(&TargetState),
          "the target quantum number per unit cell")
         ("boundary", prog_opt::value(&BoundaryState),
          "use this boundary quantum number for initializing the unit cell "
          "(useful for integer spin chains, can be used multiple times)")
         ("create,b", prog_opt::bool_switch(&NoFixedPoint),
          "Construct a new wavefunction from a random state or single-cell diagonalization")
         ("steps,s", prog_opt::value<int>(&NumSteps),
          FormatDefault("Number of DMRG steps to perform", NumSteps).c_str())
         ("no-orthogonalize", prog_opt::bool_switch(&NoOrthogonalize),
          "Don't orthogonalize the wavefunction before saving")
         ("maxiter", prog_opt::value<int>(&NumIter),
          FormatDefault("Maximum number of Lanczos iterations per step (Krylov subspace size)", NumIter).c_str())
         ("miniter", prog_opt::value<int>(&MinIter),
          FormatDefault("Minimum number of Lanczos iterations per step", MinIter).c_str())
         ("maxtol", prog_opt::value(&MaxTol),
          FormatDefault("Maximum tolerance of the eigensolver", MaxTol).c_str())
         ("fidelityscale", prog_opt::value(&FidelityScale),
          FormatDefault("The tolerance of the eigensolver is min(maxtol, fidelityscale * sqrt(fidelity))",
                        FidelityScale).c_str())
         ("initialfidelity", prog_opt::value(&InitialFidelity),
          FormatDefault("Initial value for the fidelity to set the eigensolver tolerance, for the first iteration",
                        InitialFidelity).c_str())
         ("seed", prog_opt::value<unsigned long>(), "random seed")
         ("gmrestol", prog_opt::value(&GMRESTol),
          FormatDefault("tolerance for GMRES linear solver for the initial H matrix elements", GMRESTol).c_str())
	 ("solver", prog_opt::value<std::string>(),
	  "Eigensolver to use.  Supported values are lanzcos [default], arnoldi, shift-invert")
	 ("shift-invert-energy", prog_opt::value(&ShiftInvertEnergy),
	  "For the shift-invert solver, the target energy")

         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity (can be used more than once)")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::positional_options_description p;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      std::cout << "Starting iDMRG.  Hamiltonian = " << HamStr << '\n';
      std::cout << "Wavefunction = " << FName << std::endl;

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX)
         : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      if (vm.count("one-site"))
         TwoSite = !OneSite;

      bool StartFromFixedPoint = !NoFixedPoint; // we've reversed the option

      if (!StartFromFixedPoint && !ExactDiag)
         Create = true;  // start from random state

      // The main MPWavefunction object.  We use this for initialization (if we are starting from
      // an existing wavefunction), and it will be the final wavefunction that we save to disk.
      MPWavefunction Wavefunction;

      // The parameters for the iDMRG that we need to initialize
      LinearWavefunction Psi;
      QuantumNumber QShift;

      RealDiagonalOperator R;
      MatrixOperator UR;

      // Initialize the filesystem

      if (ExactDiag || Create)
      {
         pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
      }
      else
      {
         long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
         pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize);
         Wavefunction = *PsiPtr;

         InfiniteWavefunctionLeft StartingWavefunction = Wavefunction.get<InfiniteWavefunctionLeft>();

         std::tie(Psi, R) = get_left_canonical(StartingWavefunction);
         UR = MatrixOperator::make_identity(R.Basis2());
         QShift = StartingWavefunction.qshift();
      }

      // Hamiltonian
      InfiniteLattice Lattice;
      TriangularMPO HamMPO;

      // get the Hamiltonian from the attributes, if it wasn't supplied
      if (HamStr.empty())
      {
         if (Wavefunction.Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: no Hamiltonian specified, use -H option or set wavefunction attribute Hamiltonian.\n";
            return 1;
         }
         HamStr = Wavefunction.Attributes()["Hamiltonian"].as<std::string>();
      }
      else
         Wavefunction.Attributes()["Hamiltonian"] = HamStr;

      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);
      int const UnitCellSize = Lattice.GetUnitCell().size();
      if (WavefuncUnitCellSize == 0)
         WavefuncUnitCellSize = UnitCellSize;

      optimize(HamMPO);

      // load the wavefunction
      if (ExactDiag)
      {
         QShift = QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), TargetState);
         std::cout << "Target quantum number = " << QShift << '\n';

         std::vector<BasisList> BL = ExtractLocalBasis1(HamMPO.data());
         std::vector<BasisList> FullBL = BL;
         while (int(FullBL.size()) < WavefuncUnitCellSize)
            std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));
         if (WavefuncUnitCellSize != int(FullBL.size()))
         {
            std::cout << "mp-idmrg: fatal: the wavefunction unit cell must be a multiple of the Hamiltonian unit cell.\n";
            return 1;
         }
         std::cout << "Creating exact diagonalization basis.  Wvaefunction unit cell size = "
                   << WavefuncUnitCellSize << '\n';

         QuantumNumbers::QuantumNumberList BoundaryQ;
         if (BoundaryState.empty())
         {
            BoundaryQ.push_back(QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList()));
         }
         else
         {
            for (unsigned i = 0; i < BoundaryState.size(); ++i)
            {
               std::cout << "Adding boundary quantum number " << BoundaryState[i] << '\n';
               BoundaryQ.push_back(QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), BoundaryState[i]));
            }
         }

         QuantumNumbers::QuantumNumberList LeftBoundary;
         for (unsigned i = 0; i < BoundaryQ.size(); ++i)
         {
            LeftBoundary.push_back(transform_targets(QShift, BoundaryQ[i])[0]);
         }

         if (LeftBoundary.size() == 0)
         {
            std::cerr << "fatal: the target quntum number is incompatible with the boundary quantum number"
               " for this unit cell.\n";
            return 1;
         }

         VectorBasis B1(HamMPO.front().GetSymmetryList());
         for (unsigned i = 0; i < LeftBoundary.size(); ++i)
         {
            B1.push_back(LeftBoundary[i], 1);
         }
         VectorBasis B2(HamMPO.front().GetSymmetryList());
         for (unsigned i = 0; i < BoundaryQ.size(); ++i)
         {
            B2.push_back(BoundaryQ[i], 1);
         }

         Psi.push_back(ConstructFromLeftBasis(FullBL[0], B1));
         for (int i = 1; i < WavefuncUnitCellSize; ++i)
         {
            Psi.push_back(ConstructFromLeftBasis(FullBL[i], Psi.get_back().Basis2()));
         }

         UR = MakeRandomMatrixOperator(Psi.Basis2(), B2);

         // adjust for periodic basis
         StateComponent x = prod(Psi.get_back(), UR);
         std::tie(R, UR);
         MatrixOperator X = TruncateBasis2(x); // the Basis2 is already 1-dim.  This just orthogonalizes x
         MatrixOperator U;
         SingularValueDecomposition(X, U, R, UR);
         x = prod(x, U);
         Psi.set_back(x);


         //      L = delta_shift(R, QShift);
      }
      else if (Create)
      {
         QShift = QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), TargetState);
         std::cout << "Target quantum number = " << QShift << '\n';

         std::cout << "Creating wavefunction.  Wavefunction unit cell size = " << WavefuncUnitCellSize << '\n';
         if (WavefuncUnitCellSize % UnitCellSize != 0)
         {
            std::cout << "mp-idmrg: fatal: the wavefunction unit cell must be a multiple of the Hamiltonian unit cell.\n";
            return 1;
         }
         std::vector<BasisList> BL = ExtractLocalBasis1(HamMPO.data());
         std::vector<BasisList> FullBL = BL;
         while (int(FullBL.size()) < WavefuncUnitCellSize)
            std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));

         QuantumNumber LBoundary, RBoundary;
         if (BoundaryState.empty())
         {
            RBoundary = QuantumNumber(HamMPO[0].GetSymmetryList());
            LBoundary = QShift;
         }
         else
         {
            RBoundary = QuantumNumber(HamMPO[0].GetSymmetryList(), BoundaryState[0]);
            std::cout << "Right boundary quantum number is " << RBoundary << '\n';
            if (BoundaryState.size() > 1)
            {
               std::cout << "WARNING: ignoring addititional boundary quantum numbers in random wavefunction\n";
            }
            QuantumNumbers::QuantumNumberList QL = transform_targets(QShift, RBoundary);
            if (QL.size() > 1)
            {
               PANIC("Don't know how to handle non-scalar non-abelian target state")(RBoundary)(QShift);
            }
            LBoundary = QL[0];
            std::cout << "Left boundary quantum number is " << LBoundary << '\n';
         }
         Psi = CreateRandomWavefunction(FullBL, LBoundary, 3, RBoundary);
         MatrixOperator X = left_orthogonalize(MatrixOperator::make_identity(Psi.Basis1()), Psi);
         MatrixOperator U;
         SingularValueDecomposition(X, U, R, UR);
         Psi.set_back(prod(Psi.get_back(), U));
      }

      WavefuncUnitCellSize = Psi.size();
      std::cout << "Wavefunction unit cell size = " << WavefuncUnitCellSize << '\n';
      if (WavefuncUnitCellSize % HamMPO.size() != 0)
      {
         std::cout << "mp-idmrg: fatal: the wavefunction unit cell must be a multiple of the Hamiltonian unit cell.\n";
         return 1;
      }

      if (vm.count("evolve"))
         std::cout << "Evolving with timestep " << EvolveDelta << '\n';

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = 0;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      StatesList MyStates(States);
      if (vm.count("steps") && MyStates.size() == 1)
      {
         MyStates.Repeat(NumSteps);
      }
      std::cout << MyStates << '\n';

      std::complex<double> InitialEnergy = 0.0;

      // replicate the HamMPO until it has the same size as the unit cell
      HamMPO = repeat(HamMPO, WavefuncUnitCellSize / HamMPO.size());
      CHECK_EQUAL(int(HamMPO.size()), WavefuncUnitCellSize);

      // Check that the local basis for the wavefunction and hamiltonian are compatible
      // Check that the local basis for the wavefunction and hamiltonian are compatible
      if (ExtractLocalBasis(Psi) != ExtractLocalBasis1(HamMPO))
      {
         std::cerr << "fatal: Hamiltonian is defined on a different local basis to the wavefunction.\n";
         return 1;
      }

      if (ExtractLocalBasis1(HamMPO) != ExtractLocalBasis2(HamMPO))
      {
         std::cerr << "fatal: Hamiltonian has different domain and co-domain.\n";
         return 1;
      }


      // Get the initial Hamiltonian matrix elements
      StateComponent BlockHamL = Initial_E(HamMPO, Psi.Basis1());
      if (StartFromFixedPoint)
      {
         std::cout << "Solving fixed-point Hamiltonian..." << std::endl;
         MatrixOperator Rho = delta_shift(scalar_prod(R, herm(R)), QShift);
         InitialEnergy = SolveSimpleMPO_Left(BlockHamL, Psi, QShift, HamMPO,
                                             Rho, GMRESTol, Verbose);
         std::cout << "Starting energy (left eigenvalue) = " << InitialEnergy << std::endl;

         //BlockHamL = delta_shift(BlockHamL, QShift);
      }

      StateComponent BlockHamR = Initial_F(HamMPO, Psi.Basis2());
      if (StartFromFixedPoint)
      {
         LinearWavefunction PsiR;
         MatrixOperator U;
         RealDiagonalOperator D;
         std::tie(U, D, PsiR) = get_right_canonical(Wavefunction.get<InfiniteWavefunctionLeft>());

         //TRACE(norm_frob(delta_shift(R,QShift)*U - U*D));

         // D and R are the 'same' matrix, but the basis may be in a different order
         //      TRACE(1.0-inner_prod(MatrixOperator(delta_shift(R,QShift)),MatrixOperator(D)));
         // TRACE(norm_frob(MatrixOperator(delta_shift(R,QShift))-MatrixOperator(D)));
         //TRACE(R)(D);

         MatrixOperator L = D;

#if 0
         L = triple_prod(U,L,herm(U));
         PsiR.set_front(prod(U, PsiR.get_front()));
#else
         PsiR.set_back(prod(PsiR.get_back(), delta_shift(U, adjoint(QShift))));
#endif

         BlockHamR = Initial_F(HamMPO, PsiR.Basis2());

         // check that we are orthogonalized
#if !defined(NDEBUG)
         MatrixOperator X = MatrixOperator::make_identity(PsiR.Basis2());
         X = inject_right(X, PsiR);
         CHECK(norm_frob(X - MatrixOperator::make_identity(PsiR.Basis1())) < 1E-12)(X);
#endif

         MatrixOperator Rho = L;
         Rho = D;
         //         MatrixOperator Rho = R;
         Rho = scalar_prod(Rho, herm(Rho));
#if !defined(NDEBUG)
         MatrixOperator XX = Rho;
         XX = inject_left(XX, PsiR);
         CHECK(norm_frob(delta_shift(XX,QShift) - Rho) < 1E-12)(norm_frob(delta_shift(XX,QShift) - Rho) )(XX)(Rho);
#endif

         //      BlockHamL.back() = triple_prod(herm(U), BlockHamL.back(), U);

         // We obtained Rho from the left side, so we need to delta shift to the right basis
         Rho = delta_shift(Rho, adjoint(QShift));

         std::complex<double> Energy = SolveSimpleMPO_Right(BlockHamR, PsiR, QShift, HamMPO,
                                                            Rho, GMRESTol, Verbose);
         std::cout << "Starting energy (right eigenvalue) = " << Energy << std::endl;

         //TRACE(norm_frob(delta_shift(MatrixOperator(R),QShift) - triple_prod(U,L,herm(U))));
         //      TRACE(MatrixOperator(R) - triple_prod(U,L,herm(U)));

#if 1
         U = delta_shift(U, adjoint(QShift));
         BlockHamR = prod(prod(U, BlockHamR), herm(U));

#if 0
         BlockHamL = triple_prod(herm(U), BlockHamL, U);
         Psi.set_back(prod(Psi.get_back(), U));
         Psi.set_front(prod(adjoint(U), Psi.get_front()));
         R = D;

#if 1
         StateComponent BlockHamLCheck = BlockHamL;

         BlockHamL = Initial_E(HamMPO , Psi.Basis1());
         Rho = scalar_prod(R, herm(R));
         InitialEnergy = MPO_EigenvaluesLeft(BlockHamL, Psi, QShift, HamMPO, Rho);
         std::cout << "Starting energy (left eigenvalue) = " << InitialEnergy << std::endl;

         BlockHamL = delta_shift(BlockHamL, QShift);

         TRACE(inner_prod(BlockHamLCheck - BlockHamL, Rho));
#endif
#endif
#endif

         //      TRACE(norm_frob(BlockHamL.back() - BlockHamR.front()));
         //TRACE(inner_prod(BlockHamL.back() - BlockHamR.front(), Rho));

         // Check the energy
#if 0
         MatrixOperator Cn = delta_shift(D, adjoint(QShift)); // center matrix
         Cn = U*Cn*adjoint(U);
         TRACE(norm_frob(Cn));
         Cn *= 1.0 / norm_frob(Cn);
         MatrixOperator Cnp = operator_prod(BlockHamL, Cn, herm(BlockHamR));
         std::complex<double> e = inner_prod(Cn, Cnp);
         TRACE(e);
         TRACE(norm_frob(e*Cn - Cnp));
         TRACE(inner_prod(e*Cn - Cnp, Cn));
         TRACE(norm_frob(conj(e)*Cn - Cnp));
#endif



      }

      // initialization complete.

      UR = delta_shift(UR, QShift);

      // Construct the iDMRG object
      iDMRG idmrg(Psi, R, UR, QShift, HamMPO, BlockHamL,
                  BlockHamR, InitialEnergy, Verbose);

      idmrg.MixingInfo.MixFactor = MixFactor;
      idmrg.MixingInfo.RandomMixFactor = RandomMixFactor;
      idmrg.SetInitialFidelity(InitialFidelity);
      idmrg.Solver().MaxTol = MaxTol;
      idmrg.Solver().MinTol = MinTol;
      idmrg.Solver().MinIter = MinIter;
      idmrg.Solver().MaxIter = NumIter;
      idmrg.Solver().FidelityScale = FidelityScale;
      idmrg.Solver().Verbose = Verbose;
      idmrg.Solver().EvolveDelta = EvolveDelta;
      if (vm.count("solver"))
      {
	 idmrg.Solver().SetSolver(vm["solver"].as<std::string>());
      }
      idmrg.Solver().SetShiftInvertEnergy(ShiftInvertEnergy);

      int ReturnCode = 0;

      try
      {
         bool First = true;
         for (int i = 0; i < MyStates.size(); ++i)
         {
            SInfo.MaxStates = MyStates[i].NumStates;

            if (i % 2 == 0)
            {
               idmrg.SweepLeft(SInfo, HMix, First);
               First = false;
            }
            else
            {
               idmrg.SweepRight(SInfo, HMix);
            }
         }
         idmrg.Finish(SInfo);

      }
      catch (ProcControl::Checkpoint& c)
      {
         ReturnCode = c.ReturnCode();
         std::cerr << "Early termination: "
                   << c.Reason() << '\n';
         EarlyTermination = true;
      }
      catch (...)
      {
         throw;      // if we got some other exception, don't even try and recover
      }

      // finished the iterations.
      std::cout << "Orthogonalizing wavefunction...\n";
      Wavefunction.Wavefunction() = InfiniteWavefunctionLeft::Construct(idmrg.Wavefunction(), idmrg.QShift);

      // any other attributes?
      Wavefunction.Attributes()["LastEnergy"] = idmrg.Solver().LastEnergy();
      Wavefunction.SetDefaultAttributes();

      // History log
      Wavefunction.AppendHistory(EscapeCommandline(argc, argv));

      // save wavefunction
      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Wavefunction));
      pheap::ShutdownPersistent(P);

      ProcControl::Shutdown();
      return ReturnCode;
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      if (e.Why == "File exists")
         std::cerr << "Note: use --force (-f) option to overwrite.\n";
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
