// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/itdvp.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "itdvp.h"
//#include "lanczos-exponential.h"
#include "triangular_mpo_solver.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/exponential.h"

template <typename VectorType, typename MultiplyFunctor>
VectorType LanczosExponential(VectorType const& x,
                              MultiplyFunctor MatVecMultiply, int& Iterations,
                              std::complex<double> const& Theta, double& ErrTol)
{
   typedef std::complex<double> complex;

   std::vector<VectorType> v(Iterations+1);
   std::vector<double> Alpha(Iterations+1);
   std::vector<double> Beta(Iterations);
   LinearAlgebra::Matrix<double> T(Iterations+1, Iterations+1, 0.0);
   LinearAlgebra::Vector<complex> C(1, 1.0);
   double Err = 1.0;
   double const BetaTol = 1e-14;

   v[0] = x * (1.0/norm_frob(x));

   VectorType w = MatVecMultiply(v[0]);
   Alpha[0] = real(inner_prod(w, v[0]));
   w = w - Alpha[0] * v[0];

   T(0, 0) = Alpha[0];

   int i;
   for (i = 1; i <= Iterations && Err > ErrTol; ++i)
   {
      Beta[i-1] = norm_frob(w);

      if (Beta[i-1] < BetaTol)
      {
         DEBUG_TRACE("Beta hit tolerance")(Beta[i-1]);
         break;
      }

      v[i] = w * (1.0/Beta[i-1]);

      w = MatVecMultiply(v[i]);
      Alpha[i] = real(inner_prod(w, v[i]));
      w = w - Alpha[i] * v[i] - Beta[i-1] * v[i-1];

      T(i, i) = Alpha[i];
      T(i-1, i) = T(i, i-1) = Beta[i-1];

      LinearAlgebra::Matrix<complex> P = Theta * T(LinearAlgebra::range(0, i+1),
                                                   LinearAlgebra::range(0, i+1));
      LinearAlgebra::Matrix<complex> X = LinearAlgebra::Exponentiate(1.0, P);
      LinearAlgebra::Vector<complex> CNew = X(LinearAlgebra::all, 0);

      Err = norm_2(direct_sum(C, complex(0.0)) - CNew);

      C = CNew;
   }

   Iterations = i-1;
   ErrTol = Err;

   // Calculate the time-evolved vector.
   VectorType Out = v[0] * C[0];
   for (int j = 1; j <= Iterations; ++j)
   {
      Out += v[j] * C[j];
   }

   // Normalize Out.
   Out *= 1.0/norm_frob(Out);

   return Out;
}

struct HEff1
{
   HEff1(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_)
      : E(E_), H(H_), F(F_)
   {
   }

   StateComponent operator()(StateComponent const& x) const
   {
      StateComponent R = operator_prod_inner(H, E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
};

struct HEff2
{
   HEff2(StateComponent const& E_, StateComponent const& F_)
      : E(E_), F(F_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator R = operator_prod(E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   StateComponent const& F;
};

iTDVP::iTDVP(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_,
             std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
             double GMRESTol_, double FidelityTol_, int MaxSweeps_, StatesInfo SInfo_,
             int Verbose_)
   : Hamiltonian(Ham_), Timestep(Timestep_), MaxIter(MaxIter_), ErrTol(ErrTol_),
     GMRESTol(GMRESTol_), FidelityTol(FidelityTol_), MaxSweeps(MaxSweeps_),
     SInfo(SInfo_), Verbose(Verbose_)
{
   QShift = Psi_.qshift();
   LambdaR = Psi_.lambda_r();

   // Initialize Psi and Ham.
   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;

   H = Hamiltonian.begin();

   BlockHamL = Initial_E(Hamiltonian, Psi_.Basis1());
   SolveSimpleMPO_Left(BlockHamL, Psi_, Hamiltonian, GMRESTol, Verbose-1);
   HamL.push_back(BlockHamL);

   for (InfiniteWavefunctionLeft::const_mps_iterator I = Psi_.begin(); I != Psi_.end(); ++I)
   {
      if (Verbose > 1)
         std::cout << "Site " << (HamL.size()) << std::endl;
      HamL.push_back(contract_from_left(*H, herm(*I), HamL.back(), *I));
      Psi.push_back(*I);
      MaxStates = std::max(MaxStates, (*I).Basis2().total_dimension());
      ++H;
   }

   // Calculate initial right Hamiltonian.
   LinearWavefunction PsiR;
   MatrixOperator U;
   RealDiagonalOperator D;
   std::tie(U, D, PsiR) = get_right_canonical(Psi_);

   PsiR.set_back(prod(PsiR.get_back(), delta_shift(U, adjoint(QShift))));

   BlockHamR = Initial_F(Hamiltonian, PsiR.Basis2());
   MatrixOperator Rho = scalar_prod(D, herm(D));
   Rho = delta_shift(Rho, adjoint(QShift));
   SolveSimpleMPO_Right(BlockHamR, PsiR, QShift, Hamiltonian, Rho, GMRESTol, Verbose-1);

   HamR.push_front(BlockHamR);

   // Initialize to the right-most site.
   HamL.pop_back();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size() - 1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;
}

InfiniteWavefunctionLeft iTDVP::Wavefunction() const
{
   return InfiniteWavefunctionLeft::ConstructFromOrthogonal(Psi, LambdaR, QShift);
}

void iTDVP::IterateLeft()
{
   // Evolve current site.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, 0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }

   // Perform SVD to right-orthogonalize current site.
   MatrixOperator M = ExpandBasis1(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   *C = prod(Vh, *C);

   // Update the effective Hamiltonian.
   HamR.push_front(contract_from_right(herm(*H), *C, HamR.front(), herm(*C)));

   // Evolve the UD term backwards in time.
   Iter = MaxIter;
   Err = ErrTol;
   MatrixOperator UD = U*D;

   UD = LanczosExponential(UD, HEff2(HamL.back(), HamR.front()), Iter, -0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " UD Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }

   // Move to the next site.
   --Site;
   --H;
   --C;

   *C = prod(*C, UD);

   HamL.pop_back();
}

void iTDVP::EvolveLeftmostSite()
{
   // Evolve current site.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, 0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void iTDVP::IterateRight()
{
   // Perform SVD to left-orthogonalize current site.
   MatrixOperator M = ExpandBasis2(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   *C = prod(*C, U);

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));

   // Evolve the DVh term backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;
   MatrixOperator DVh = D*Vh;

   DVh = LanczosExponential(DVh, HEff2(HamL.back(), HamR.front()), Iter, -0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " DVh Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }

   // Move to the next site.
   ++Site;
   ++H;
   ++C;

   *C = prod(DVh, *C);

   HamR.pop_front();

   // Evolve current site.
   Iter = MaxIter;
   Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, 0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void iTDVP::OrthogonalizeLeftmostSite()
{
   // Right-orthogonalize current site.
   MatrixOperator M = ExpandBasis1(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(U*Vh, *C);
   LambdaR = (U*D)*herm(U);

   // Update right block Hamiltonian.
   BlockHamR = contract_from_right(herm(*H), *C, HamR.front(), herm(*C));
   HamR.back() = BlockHamR;
}

void iTDVP::OrthogonalizeRightmostSite()
{
   // Left-orthogonalize current site.
   MatrixOperator M = ExpandBasis2(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(*C, U*Vh);
   LambdaR = herm(Vh)*(D*Vh);

   // Update left block Hamiltonian.
   BlockHamL = contract_from_left(*H, herm(*C), HamL.back(), *C);
   HamL.front() = BlockHamL;
}

void iTDVP::EvolveLambdaR()
{
   int Iter = MaxIter;
   double Err = ErrTol;
   LambdaR = LanczosExponential(LambdaR, HEff2(BlockHamL, BlockHamR), Iter, -0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " LambdaR Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void iTDVP::Evolve()
{
   ++TStep;

   // Sweep left to evolve the unit cell for half a time step.
   SweepL = 0;
   FidelityL = 1.0;

   LinearWavefunction PsiOld = Psi;
   std::deque<StateComponent> HamLOld = HamL;

   do {
      ++SweepL;

      LinearWavefunction PsiPrev = Psi;
      MatrixOperator LambdaRPrev = LambdaR;

      Psi = PsiOld;
      C = Psi.end();
      --C;
      H = Hamiltonian.end();
      --H;
      Site = RightStop;

      HamL = HamLOld;
      HamR = std::deque<StateComponent>(1, BlockHamR);

      if (SweepL > 1)
         this->EvolveLambdaR();

      *C = prod(*C, LambdaR);

      while (Site > LeftStop)
         this->IterateLeft();

      this->EvolveLeftmostSite();

      this->OrthogonalizeLeftmostSite();

      if (SweepL > 1)
      {
         MatrixOperator Rho = scalar_prod(herm(LambdaRPrev), LambdaR);
         FidelityL = 1.0 - sum(SingularValues(inject_left(Rho, Psi, PsiPrev)));

         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepL=" << SweepL
                      << " FidelityL=" << FidelityL
                      << std::endl;
         }
      }
      else
      {
         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepL=" << SweepL
                      << std::endl;
         }
      }
   }
   while (FidelityL > FidelityTol && SweepL < MaxSweeps);

   // Sweep right to evolve the unit cell for half a time step.
   SweepR = 0;
   FidelityR = 1.0;

   PsiOld = Psi;
   std::deque<StateComponent> HamROld = HamR;

   do {
      ++SweepR;

      LinearWavefunction PsiPrev = Psi;
      MatrixOperator LambdaRPrev = LambdaR;

      Psi = PsiOld;
      C = Psi.begin();
      H = Hamiltonian.begin();
      Site = LeftStop;

      HamL = std::deque<StateComponent>(1, BlockHamL);
      HamR = HamROld;

      if (SweepR > 1)
         this->EvolveLambdaR();

      *C = prod(LambdaR, *C);

      this->EvolveLeftmostSite();

      while (Site < RightStop)
         this->IterateRight();

      this->OrthogonalizeRightmostSite();

      if (SweepR > 1)
      {
         MatrixOperator Rho = scalar_prod(LambdaR, herm(LambdaRPrev));
         FidelityR = 1.0 - sum(SingularValues(inject_right(Rho, Psi, PsiPrev)));

         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepR=" << SweepR
                      << " FidelityR=" << FidelityR
                      << std::endl;
         }
      }
      else
      {
         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepR=" << SweepR
                      << std::endl;
         }
      }
   }
   while (FidelityR > FidelityTol && SweepR < MaxSweeps);

   this->CalculateEps();
}

#if 0
void iTDVP::CalculateEps()
{
   X = std::deque<StateComponent>();
   Y = std::deque<StateComponent>();

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   // Create a local copy of Psi so we can perform a single right-to-left sweep
   // without having to go back to the right.
   LinearWavefunction PsiLocal = Psi;
   LinearWavefunction::iterator CLocal = Psi.end();
   --CLocal;
   BasicTriangularMPO::const_iterator HLocal = Hamiltonian.end();
   --HLocal;
   std::deque<StateComponent> HamLLocal = HamL;
   std::deque<StateComponent> HamRLocal = HamR;
   int SiteLocal = RightStop;

   *CLocal = prod(*CLocal, LambdaR);

   // Calcaulate the left half of epsilon_2 for the right end of the unit cell.
   X.push_front(contract_from_left(*HLocal, herm(NullSpace2(*CLocal)), HamLLocal.back(), *CLocal));

   while (SiteLocal > LeftStop)
   {
      // Perform SVD to right-orthogonalize current site.
      MatrixOperator M = ExpandBasis1(*CLocal);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      *CLocal = prod(Vh, *CLocal);

      // Calculate the right half of epsilon_2.
      Y.push_front(contract_from_right(herm(*HLocal), NullSpace1(*CLocal), HamRLocal.front(), herm(*CLocal)));

      // Update the effective Hamiltonian.
      HamRLocal.push_front(contract_from_right(herm(*HLocal), *CLocal, HamRLocal.front(), herm(*CLocal)));

      // Move to the next site.
      --SiteLocal;
      --HLocal;
      --CLocal;

      *CLocal = prod(*CLocal, U*D);

      HamLLocal.pop_back();

      // Calculate the left half of epsilon_2.
      X.push_front(contract_from_left(*HLocal, herm(NullSpace2(*CLocal)), HamLLocal.back(), *CLocal));

      // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
      double Eps1Sq = norm_frob_sq(scalar_prod(X.front(), herm(HamRLocal.front())));
      double Eps2Sq = norm_frob_sq(scalar_prod(X.front(), herm(Y.front())));
      Eps1SqSum += Eps1Sq;
      Eps2SqSum += Eps2Sq;

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << SiteLocal
                   << " Eps1Sq=" << Eps1Sq
                   << " Eps2Sq=" << Eps2Sq
                   << std::endl;
      }
   }

   // Perform SVD to right-orthogonalize the leftmost site.
   MatrixOperator M = ExpandBasis1(*CLocal);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   *CLocal = prod(U*Vh, *CLocal);

   // Calculate the right half of epsilon_2.
   Y.push_back(contract_from_right(herm(*HLocal), NullSpace1(*CLocal), HamRLocal.front(), herm(*CLocal)));

   // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
   double Eps1Sq = norm_frob_sq(scalar_prod(X.back(), herm(HamRLocal.back())));
   double Eps2Sq = norm_frob_sq(scalar_prod(X.back(), herm(Y.back())));
   Eps1SqSum += Eps1Sq;
   Eps2SqSum += Eps2Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << RightStop
                << " Eps1Sq=" << Eps1Sq
                << " Eps2Sq=" << Eps2Sq
                << std::endl;
   }
}
#endif

void iTDVP::CalculateEps()
{
   X = std::deque<StateComponent>();
   Y = std::deque<StateComponent>();

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   *C = prod(*C, LambdaR);

   // Right-orthogonalize the unit cell.
   while (Site > LeftStop)
   {
      // Perform SVD to right-orthogonalize current site.
      MatrixOperator M = ExpandBasis1(*C);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      *C = prod(Vh, *C);

      // Update the effective Hamiltonian.
      HamR.push_front(contract_from_right(herm(*H), *C, HamR.front(), herm(*C)));

      // Move to the next site.
      --Site;
      --H;
      --C;

      *C = prod(*C, U*D);

      HamL.pop_back();
   }

   // Calcaulate the right half of epsilon_2 for the left end of the unit cell.
   Y.push_back(contract_from_right(herm(*H), NullSpace1(*C), HamR.front(), herm(*C)));

   while (Site < RightStop)
   {
      // Perform SVD to left-orthogonalize current site.
      MatrixOperator M = ExpandBasis2(*C);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      *C = prod(*C, U);

      // Calculate the left half of epsilon_2.
      X.push_back(contract_from_left(*H, herm(NullSpace2(*C)), HamL.back(), *C));

      // Update the effective Hamiltonian.
      HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));

      // Move to the next site.
      ++Site;
      ++H;
      ++C;

      *C = prod(D*Vh, *C);

      HamR.pop_front();

      // Calculate the right half of epsilon_2.
      Y.push_back(contract_from_right(herm(*H), NullSpace1(*C), HamR.front(), herm(*C)));

      // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
      double Eps1Sq = norm_frob_sq(scalar_prod(HamL.back(), herm(Y.back())));
      double Eps2Sq = norm_frob_sq(scalar_prod(X.back(), herm(Y.back())));
      Eps1SqSum += Eps1Sq;
      Eps2SqSum += Eps2Sq;

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << Site
                   << " Eps1Sq=" << Eps1Sq
                   << " Eps2Sq=" << Eps2Sq
                   << std::endl;
      }
   }

   // Perform SVD to left-orthogonalize the rightmost site.
   MatrixOperator M = ExpandBasis2(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(*C, U*Vh);
   LambdaR = herm(Vh)*(D*Vh);

   // Update left block Hamiltonian.
   BlockHamL = contract_from_left(*H, herm(*C), HamL.back(), *C);
   HamL.front() = BlockHamL;

   // Calculate the left half of epsilon_2.
   X.push_back(contract_from_left(*H, herm(NullSpace2(*C)), HamL.back(), *C));

   Y.push_back(Y.front());
   Y.pop_front();

   // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
   double Eps1Sq = norm_frob_sq(scalar_prod(HamL.front(), herm(Y.back())));
   double Eps2Sq = norm_frob_sq(scalar_prod(X.back(), herm(Y.back())));
   Eps1SqSum += Eps1Sq;
   Eps2SqSum += Eps2Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << LeftStop
                << " Eps1Sq=" << Eps1Sq
                << " Eps2Sq=" << Eps2Sq
                << std::endl;
   }
}

void iTDVP::ExpandRightBond()
{
   // Take the truncated SVD of P_2 H|Psi>.
   CMatSVD SL(scalar_prod(X.front(), herm(Y.front())));
   TruncationInfo Info;
   StatesInfo SInfoLocal = SInfo;
   // Subtract the current bond dimension from the number of additional states to be added.
   SInfoLocal.MinStates = std::max(0, SInfoLocal.MinStates - (*C).Basis2().total_dimension());
   SInfoLocal.MaxStates = std::max(0, SInfoLocal.MaxStates - (*C).Basis2().total_dimension());
   CMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfoLocal, Info);

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SL.ConstructMatrices(SL.begin(), Cutoff, U, D, Vh);

   // Construct new basis.
   SumBasis<VectorBasis> NewBasis((*C).Basis2(), U.Basis2());
   // Construct a unitary to regularize the new basis.
   MatrixOperator UReg = Regularize(NewBasis);

   MaxStates = std::max(MaxStates, NewBasis.total_dimension());

   // Add the new states to CExpand.
   StateComponent CExpand = tensor_row_sum(*COld, prod(NullSpace2(*COld), U), NewBasis);
   CExpand = prod(CExpand, herm(UReg));

   // Add a zero tensor of the same dimensions to C.
   StateComponent Z = StateComponent((*C).LocalBasis(), (*C).Basis1(), U.Basis2());
   *C = tensor_row_sum(*C, Z, NewBasis);
   *C = prod(*C, herm(UReg));

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " NewStates=" << Info.KeptStates()
                << " TotalStates=" << NewBasis.total_dimension()
                << std::endl;
   }

   // Move to next site, handling the rightmost site separately.
   if (Site < RightStop)
   {
      HamL.push_back(contract_from_left(*H, herm(CExpand), HamLOld.front(), CExpand));

      HamLOld.pop_front();
      ++C;
      ++COld;
   }
   else
   {
      BlockHamL = contract_from_left(*H, herm(CExpand), HamLOld.front(), CExpand);
      HamL.front() = BlockHamL;

      C = Psi.begin();

      // Add zeros to LambdaR so the bonds match. 
      MatrixOperator Z = MatrixOperator(Vh.Basis1(), LambdaR.Basis2());
      LambdaR = tensor_col_sum(LambdaR, Z, NewBasis);
      LambdaR = UReg * LambdaR;
   }

   ++Site;
   ++H;
   X.pop_front();
   Y.pop_front();

   // Add zeros to the current site so the left bond matches. 
   Z = StateComponent((*C).LocalBasis(), Vh.Basis1(), (*C).Basis2());
   *C = tensor_col_sum(*C, Z, NewBasis);
   *C = prod(UReg, *C);
}

void iTDVP::ExpandBonds()
{
   C = Psi.begin();
   H = Hamiltonian.begin();
   Site = LeftStop;
   HamLOld = HamL;
   HamL = std::deque<StateComponent>(1, BlockHamL);

   LinearWavefunction PsiOld = Psi;
   COld = PsiOld.begin();

   while (Site <= RightStop)
      this->ExpandRightBond();

   C = Psi.end();
   --C;
   H = Hamiltonian.end();
   --H;
   Site = RightStop;
}
