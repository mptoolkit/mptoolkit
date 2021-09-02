// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/itdvp-parallel.cpp
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

#include "itdvp-parallel.h"
//#include "lanczos-exponential.h"
#include "triangular_mpo_solver.h"
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

   // Normalise Out.
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
             double GMRESTol_, double OverlapTol_, int MaxSweeps_, StatesInfo SInfo_,
             int Verbose_)
   : Hamiltonian(Ham_), Timestep(Timestep_), MaxIter(MaxIter_), ErrTol(ErrTol_),
     GMRESTol(GMRESTol_), OverlapTol(OverlapTol_), MaxSweeps(MaxSweeps_),
     SInfo(SInfo_), Verbose(Verbose_)
{
   QShift = Psi_.qshift();

   // Initialise left- and right-canonical PsiAL and PsiAR
#if 0
   RealDiagonalOperator DR, DL;
   MatrixOperator U;
   std::tie(PsiAL, DR) = get_left_canonical(Psi_);
   std::tie(U, DL, PsiAR) = get_right_canonical(Psi_);
   PsiAR.set_front(prod(U, PsiAR.get_front()));
   PsiAR.set_back(prod(PsiAR.get_back(), delta_shift(U, adjoint(QShift))));
#else
   RealDiagonalOperator DR;
   std::tie(PsiAL, DR) = get_left_canonical(Psi_);
#endif

   // Initialise AC and C tensors.
   InfiniteWavefunctionLeft::const_mps_iterator A = Psi_.begin();
   InfiniteWavefunctionLeft::const_lambda_iterator Lambda = Psi_.lambda_begin();
   ++Lambda;
   while (A != Psi_.end())
   {
      PsiAC.push_back(prod(*A, *Lambda));
      PsiC.push_back(*Lambda);

      MaxStates = std::max(MaxStates, (*A).Basis1().total_dimension());

      ++A, ++Lambda;
   }
#if 0
   LinearWavefunction::const_iterator AR = PsiAR.begin();
   InfiniteWavefunctionLeft::const_lambda_iterator Lambda = Psi_.lambda_begin();
   while (AR != PsiAR.end())
   {
      PsiAC.push_back(prod(*Lambda, *AR));
      ++Lambda;
      PsiC.push_back(*Lambda);

      MaxStates = std::max(MaxStates, (*AR).Basis1().total_dimension());

      ++AR;
   }
   LinearWavefunction::const_iterator AL = PsiAL.begin();
   LinearWavefunction::const_iterator AR = PsiAR.begin();
   InfiniteWavefunctionLeft::const_lambda_iterator LambdaL = Psi_.lambda_begin();
   InfiniteWavefunctionLeft::const_lambda_iterator LambdaR = Psi_.lambda_begin();
   ++LambdaR;
   while (AR != PsiAR.end())
   {
      PsiAC.push_back(0.5*prod(*LambdaL, *AR) + 0.5*prod(*AL, *LambdaR));
      PsiC.push_back(*LambdaR);

      MaxStates = std::max(MaxStates, (*AL).Basis1().total_dimension());

      ++AL, ++AR, ++LambdaL, ++LambdaR;
   }
#endif

   PsiAR = PsiAL;
   this->FindAR();

   // Calculate initial left and right Hamiltonians.
   HamL = std::deque<StateComponent>(1, Initial_E(Hamiltonian, PsiAL.Basis1()));
   HamR = std::deque<StateComponent>(1, Initial_F(Hamiltonian, PsiAR.Basis2()));
   this->UpdateHam();
}

InfiniteWavefunctionLeft iTDVP::Wavefunction() const
{
   return InfiniteWavefunctionLeft::Construct(PsiAL, QShift);
}

void iTDVP::UpdateHam()
{
   //MatrixOperator Rho = delta_shift(scalar_prod(PsiC.back(), herm(PsiC.back())), QShift);
   MatrixOperator Rho = PsiC.back() * herm(PsiC.back());
   //StateComponent BlockHamL = HamL.front();
   StateComponent BlockHamL = Initial_E(Hamiltonian, PsiAL.Basis1());
   SolveSimpleMPO_Left(BlockHamL, PsiAL, QShift, Hamiltonian, Rho, GMRESTol, Verbose-1);
   HamL = std::deque<StateComponent>(1, BlockHamL);

   LinearWavefunction::const_iterator AL = PsiAL.begin();
   BasicTriangularMPO::const_iterator H = Hamiltonian.begin();
   while (AL != PsiAL.end())
   {
      HamL.push_back(contract_from_left(*H, herm(*AL), HamL.back(), *AL));
      ++AL, ++H;
   }

   //Rho = delta_shift(scalar_prod(PsiC.back(), herm(PsiC.back())), adjoint(QShift));
   Rho = herm(PsiC.back()) * PsiC.back();
   //StateComponent BlockHamR = HamR.back();
   StateComponent BlockHamR = Initial_F(Hamiltonian, PsiAR.Basis2());
   SolveSimpleMPO_Right(BlockHamR, PsiAR, QShift, Hamiltonian, Rho, GMRESTol, Verbose-1);
   HamR = std::deque<StateComponent>(1, BlockHamR);

   LinearWavefunction::const_iterator AR = PsiAR.end();
   H = Hamiltonian.end();
   while (AR != PsiAR.begin())
   {
      --AR, --H;
      HamR.push_front(contract_from_right(herm(*H), *AR, HamR.front(), herm(*AR)));
   }
}

void iTDVP::FindAL()
{
   EpsLSqSum = 0.0;

   std::deque<StateComponent>::iterator AC = PsiAC.begin();
   std::deque<MatrixOperator>::iterator C = PsiC.begin();
   LinearWavefunction::iterator AL = PsiAL.begin();
   int Site = 0;

   while (AC != PsiAC.end())
   {
#if 0
      *AL = prod(*AC, herm(*C));

      MatrixOperator M = ExpandBasis2(*AL);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);
      
      *AL = prod(*AL, U*Vh);
#else
      *AL = *AC;

      MatrixOperator M = ExpandBasis2(*AL);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      *AL = prod(*AL, U*Vh);

      SingularValueDecomposition(*C, U, D, Vh);

      *AL = prod(*AL, herm(U*Vh));
#endif

      double EpsLSq = norm_frob_sq(*AC - prod(*AL, *C));
      EpsLSqSum += EpsLSq;

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << Site
                   << " AL EpsLSq=" << EpsLSq
                   << std::endl;
      }

      ++AC, ++C, ++AL, ++Site;
   }
}

void iTDVP::FindAR()
{
   EpsRSqSum = 0.0;

   std::deque<StateComponent>::iterator AC = PsiAC.begin();
   std::deque<MatrixOperator>::iterator C = PsiC.end();
   --C;
   LinearWavefunction::iterator AR = PsiAR.begin();
   int Site = 0;

   while (AC != PsiAC.end())
   {
#if 0
      *AR = prod(herm(*C), *AC);

      MatrixOperator M = ExpandBasis1(*AR);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);
      
      *AR = prod(U*Vh, *AR);
#else
      *AR = *AC;

      MatrixOperator M = ExpandBasis1(*AR);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      *AR = prod(U*Vh, *AR);

      SingularValueDecomposition(*C, U, D, Vh);

      *AR = prod(herm(U*Vh), *AR);
#endif
      
      double EpsRSq = norm_frob_sq(*AC - prod(*C, *AR));
      EpsRSqSum += EpsRSq;

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << Site
                   << " AR EpsRSq=" << EpsRSq
                   << std::endl;
      }

      ++AC, ++C, ++AR, ++Site;
      if (C == PsiC.end())
         C = PsiC.begin();
   }
}

void iTDVP::Evolve()
{
   ++TStep;

   // Evolve the AC matrices forward in time.
   std::deque<StateComponent>::iterator AC = PsiAC.begin();
   std::deque<StateComponent>::iterator HL = HamL.begin();
   std::deque<StateComponent>::iterator HR = HamR.begin();
   ++HR;
   BasicTriangularMPO::const_iterator H = Hamiltonian.begin();
   int Site = 0;
   while (AC != PsiAC.end())
   {
      int Iter = MaxIter;
      double Err = ErrTol;

      *AC = LanczosExponential(*AC, HEff1(*HL, *H, *HR), Iter, Timestep, Err);

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << Site
                   << " AC Iter=" << Iter
                   << " Err=" << Err
                   << std::endl;
      }

      ++AC, ++HL, ++HR, ++H, ++Site;
   }

   // Evolve the C matrices forward in time.
   std::deque<MatrixOperator>::iterator C = PsiC.begin();
   HL = HamL.begin();
   HR = HamR.begin();
   ++HL, ++HR;
   Site = 0;
   while (C != PsiC.end())
   {
      int Iter = MaxIter;
      double Err = ErrTol;

      *C = LanczosExponential(*C, HEff2(*HL, *HR), Iter, Timestep, Err);

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << Site
                   << " C Iter=" << Iter
                   << " Err=" << Err
                   << std::endl;
      }

      ++C, ++HL, ++HR, ++Site;
   }

   this->FindAL();
   this->FindAR();

   // TEST: reform AC after each timestep.
#if 0
   if (TStep % 2 == 0)
   {
      std::deque<StateComponent>::iterator AC = PsiAC.begin();
      std::deque<MatrixOperator>::iterator C = PsiC.begin();
      AL = PsiAL.begin();
      while (AC != PsiAC.end())
      {
         *AC = prod(*AL, *C);
         ++AC, ++C, ++AL;
      }
   }
   else
   {
      std::deque<StateComponent>::iterator AC = PsiAC.begin();
      std::deque<MatrixOperator>::iterator C = PsiC.end();
      --C;
      AR = PsiAR.begin();
      while (AC != PsiAC.end())
      {
         *AC = prod(*C, *AR);
         ++AC, ++C, ++AR;
         if (C == PsiC.end())
            C = PsiC.begin();
      }
   }
   AC = PsiAC.begin();
   std::deque<MatrixOperator>::iterator CL = PsiC.end();
   --CL;
   std::deque<MatrixOperator>::iterator CR = PsiC.begin();
   AL = PsiAL.begin();
   AR = PsiAR.begin();
   while (AC != PsiAC.end())
   {
      *AC = 0.5*prod(*AL, *CR) + 0.5*prod(*CL, *AR);
      ++AC, ++CL, ++CR, ++AL, ++AR;
      if (CL == PsiC.end())
         CL = PsiC.begin();
   }
#endif

   this->UpdateHam();

   this->CalculateEps();
}

void iTDVP::CalculateEps()
{
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::deque<MatrixOperator>::iterator C = PsiC.end();
   --C;
   LinearWavefunction::iterator AL = PsiAL.end();
   --AL;
   LinearWavefunction::iterator AR = PsiAR.begin();
   std::deque<StateComponent>::iterator HL = HamL.end();
   --HL, --HL;
   std::deque<StateComponent>::iterator HR = HamR.begin();
   ++HR;
   BasicTriangularMPO::const_iterator H = Hamiltonian.end();
   --H;
   int Site = 0;

   while (AR != PsiAR.end())
   {
      StateComponent X = contract_from_left(*H, herm(NullSpace2(*AL)), *HL, *AL);
      ++H;
      if (H == Hamiltonian.end())
         H = Hamiltonian.begin();
      StateComponent Y = contract_from_right(herm(*H), NullSpace1(*AR), *HR, herm(prod(*C, *AR)));

      double Eps1Sq = norm_frob_sq(scalar_prod(*HL, herm(Y)));
      double Eps2Sq = norm_frob_sq(scalar_prod(X, herm(Y)));
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

      ++C, ++AL, ++AR, ++HL, ++HR, ++Site;
      if (AL == PsiAL.end())
      {
         C = PsiC.begin();
         AL = PsiAL.begin();
         HL = HamL.begin();
      }
   }
}

void iTDVP::ExpandBonds()
{
   MaxStates = 1;

   std::deque<StateComponent>::iterator AC = PsiAC.begin();
   std::deque<MatrixOperator>::iterator C = PsiC.end();
   --C;
   LinearWavefunction::iterator AL = PsiAL.end();
   --AL;
   LinearWavefunction::iterator AR = PsiAR.begin();
   std::deque<StateComponent>::iterator HL = HamL.end();
   --HL, --HL;
   std::deque<StateComponent>::iterator HR = HamR.begin();
   ++HR;
   BasicTriangularMPO::const_iterator H = Hamiltonian.end();
   --H;
   int Site = 0;

   std::deque<VectorBasis> AddBasis;

   while (AC != PsiAC.end())
   {
      StateComponent NL = NullSpace2(*AL);
      StateComponent NR = NullSpace1(*AR);

      StateComponent X = contract_from_left(*H, herm(NL), *HL, *AL);
      ++H;
      if (H == Hamiltonian.end())
         H = Hamiltonian.begin();
      StateComponent Y = contract_from_right(herm(*H), NR, *HR, herm(prod(*C, *AR)));

      // Calculate the truncated SVD of XY to find the most relevant states to be added.
      CMatSVD SL(scalar_prod(X, herm(Y)));
      TruncationInfo Info;
      StatesInfo SInfoLocal = SInfo;
      // Subtract the current bond dimension from the number of additional states to be added.
      SInfoLocal.MinStates = std::max(0, SInfoLocal.MinStates - (*AL).Basis2().total_dimension());
      SInfoLocal.MaxStates = std::max(0, SInfoLocal.MaxStates - (*AL).Basis2().total_dimension());
      CMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfoLocal, Info);

      MatrixOperator U, Vh;
      RealDiagonalOperator D;
      SL.ConstructMatrices(SL.begin(), Cutoff, U, D, Vh);

      // Add the new states to AL and AR.
      AddBasis.push_back(U.Basis2());
      SumBasis<VectorBasis> NewBasisL((*AL).Basis2(), AddBasis.back());
      SumBasis<VectorBasis> NewBasisR((*AR).Basis1(), AddBasis.back());

      MaxStates = std::max(MaxStates, NewBasisL.total_dimension());

      *AL = tensor_row_sum(*AL, prod(NL, U), NewBasisL);
      *AR = tensor_col_sum(*AR, prod(Vh, NR), NewBasisR);

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Bond=" << Site
                   << " NewStates=" << Info.KeptStates()
                   << " TotalStates=" << NewBasisL.total_dimension()
                   << std::endl;
      }

      ++AC, ++C, ++AL, ++AR, ++HL, ++HR, ++Site;
      if (AL == PsiAL.end())
      {
         C = PsiC.begin();
         AL = PsiAL.begin();
         HL = HamL.begin();
      }
   }

   // Fill the rest of AC, C, AL and AR with zeros.
   AC = PsiAC.begin();
   C = PsiC.begin();
   AL = PsiAL.begin();
   AR = PsiAR.begin();
   std::deque<VectorBasis>::iterator BasisL = AddBasis.begin();
   std::deque<VectorBasis>::iterator BasisR = AddBasis.begin();
   ++BasisR;
   Site = 0;

   while (AC != PsiAC.end())
   {
      if (BasisR == AddBasis.end())
         BasisR = AddBasis.begin();

      SumBasis<VectorBasis> NewBasisL((*AC).Basis1(), *BasisL);
      SumBasis<VectorBasis> NewBasisR((*AC).Basis2(), *BasisR);

      StateComponent Zeros1 = StateComponent((*AC).LocalBasis(), *BasisL, *BasisR);
      *AC = tensor_sum(*AC, Zeros1, NewBasisL, NewBasisR);

      SumBasis<VectorBasis> NewBasisLC((*C).Basis1(), *BasisR);

      MatrixOperator Zeros2 = MatrixOperator(*BasisR, *BasisR);
      *C = tensor_sum(*C, Zeros2, NewBasisLC, NewBasisR);

      StateComponent Zeros3 = StateComponent((*AL).LocalBasis(), *BasisL, (*AL).Basis2());
      *AL = tensor_col_sum(*AL, Zeros3, NewBasisL);

      StateComponent Zeros4 = StateComponent((*AR).LocalBasis(), (*AR).Basis1(), *BasisR);
      *AR = tensor_row_sum(*AR, Zeros4, NewBasisR);

      ++AC, ++C, ++AL, ++AR, ++BasisL, ++BasisR, ++Site;
   }

   this->UpdateHam();
}
