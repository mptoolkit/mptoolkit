// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/tdvp.cpp
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

#include "tdvp.h"
//#include "lanczos-exponential.h"
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
         break;

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

TDVP::TDVP(BasicTriangularMPO const& Ham_, std::complex<double> Timestep_,
           int MaxIter_, double ErrTol_, StatesInfo SInfo_, int Verbose_)
   : Hamiltonian(Ham_), Timestep(Timestep_), MaxIter(MaxIter_), ErrTol(ErrTol_),
     SInfo(SInfo_), Verbose(Verbose_)
{
}

TDVP::TDVP(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_,
           std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
           StatesInfo SInfo_, int Verbose_)
   : TDVP(Ham_, Timestep_, MaxIter_, ErrTol_, SInfo_, Verbose_)
{
   // Initialize Psi and Ham.
   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;
   H = Hamiltonian.begin();
   HamL.push_back(Initial_E(Hamiltonian, Psi_.Basis1()));
   for (FiniteWavefunctionLeft::const_mps_iterator I = Psi_.begin(); I != Psi_.end(); ++I)
   {
      if (Verbose > 1)
         std::cout << "Site " << (HamL.size()) << std::endl;
      HamL.push_back(contract_from_left(*H, herm(*I), HamL.back(), *I));
      Psi.push_back(*I);
      MaxStates = std::max(MaxStates, (*I).Basis2().total_dimension());
      ++H;
   }
   HamR.push_front(Initial_F(Hamiltonian, Psi_.Basis2()));
   Psi.set_back(prod(Psi.get_back(), Psi_.lambda_r()));

   // Initialize to the right-most site.
   HamL.pop_back();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size() - 1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;
}

FiniteWavefunctionLeft
TDVP::Wavefunction() const
{
   return FiniteWavefunctionLeft::Construct(Psi);
}

std::complex<double>
TDVP::Energy() const
{
   return inner_prod(contract_from_left(*H, herm(*C), HamL.back(), *C), HamR.front());
}

void
TDVP::IterateLeft()
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

void
TDVP::EvolveLeftmostSite()
{
   // Evolve current site.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, Timestep, Err);

   // Calculate error measure epsilon_1 and add to sum.
   StateComponent Y = contract_from_right(herm(*H), NullSpace1(*C), HamR.front(), herm(*C));
   double Eps1Sq = norm_frob_sq(scalar_prod(HamL.back(), herm(Y)));
   Eps1SqSum += Eps1Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << " Eps1Sq=" << Eps1Sq
                << std::endl;
   }
}

void
TDVP::IterateRight()
{
   // Perform SVD to left-orthogonalize current site.
   MatrixOperator M = ExpandBasis2(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   *C = prod(*C, U);

   // Calculate the left half of epsilon_2 (see below).
   StateComponent X = contract_from_left(*H, herm(NullSpace2(*C)), HamL.back(), *C);

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

   // Perform SVD to right-orthogonalize current site for NullSpace1.
   M = ExpandBasis1(*C);

   SingularValueDecomposition(M, U, D, Vh);

   StateComponent CRightOrtho = prod(Vh, *C);
   *C = prod(U*D*Vh, *C);

   // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
   StateComponent Y = contract_from_right(herm(*H), NullSpace1(CRightOrtho), HamR.front(), herm(*C));
   double Eps1Sq = norm_frob_sq(scalar_prod(HamL.back(), herm(Y)));
   double Eps2Sq = norm_frob_sq(scalar_prod(X, herm(Y)));
   Eps1SqSum += Eps1Sq;
   Eps2SqSum += Eps2Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << " Eps1Sq=" << Eps1Sq
                << " Eps2Sq=" << Eps2Sq
                << std::endl;
   }
}

void
TDVP::Evolve()
{
   ++TStep;
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   while (Site > LeftStop)
      this->IterateLeft();

   this->EvolveLeftmostSite();

   while (Site < RightStop)
      this->IterateRight();
}

void
TDVP::ExpandLeftBond()
{
   // Construct the projection of H|Psi> onto the space of two-site variations.
   LinearWavefunction::iterator CNext = C;
   --CNext;
   BasicTriangularMPO::const_iterator HNext = H;
   --HNext;

   HamL.pop_back();

   StateComponent NL = NullSpace2(*CNext);

   // Perform SVD to right-orthogonalize current site for NullSpace1.
   MatrixOperator M = ExpandBasis1(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   StateComponent CRightOrtho = prod(Vh, *C);
   *C = prod(U*D*Vh, *C);

   StateComponent NR = NullSpace1(CRightOrtho);

   StateComponent X = contract_from_left(*HNext, herm(NL), HamL.back(), *CNext);
   StateComponent Y = contract_from_right(herm(*H), NR, HamR.front(), herm(*C));

   // Take the truncated SVD of P_2 H|Psi>.
   CMatSVD SL(scalar_prod(X, herm(Y)));
   TruncationInfo Info;
   StatesInfo SInfoLocal = SInfo;
   // Subtract the current bond dimension from the number of additional states to be added.
   SInfoLocal.MinStates = std::max(0, SInfoLocal.MinStates - (*C).Basis1().total_dimension());
   SInfoLocal.MaxStates = std::max(0, SInfoLocal.MaxStates - (*C).Basis1().total_dimension());
   CMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfoLocal, Info);

   SL.ConstructMatrices(SL.begin(), Cutoff, U, D, Vh);

   // Construct new basis.
   SumBasis<VectorBasis> NewBasis((*CNext).Basis2(), U.Basis2());
   // Construct a unitary to regularize the new basis.
   MatrixOperator UReg = Regularize(NewBasis);

   MaxStates = std::max(MaxStates, NewBasis.total_dimension());

   // Add the new states to CNext, and add zeros to C.
   *CNext = tensor_row_sum(*CNext, prod(NL, U), NewBasis);
   *CNext = prod(*CNext, herm(UReg));

   StateComponent Z = StateComponent((*C).LocalBasis(), Vh.Basis1(), (*C).Basis2());
   *C = tensor_col_sum(*C, Z, NewBasis);
   *C = prod(UReg, *C);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " NewStates=" << Info.KeptStates()
                << " TotalStates=" << NewBasis.total_dimension()
                << std::endl;
   }

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*HNext, herm(*CNext), HamL.back(), *CNext));
}

void
TDVP::EvolveExpand()
{
   ++TStep;
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   while (Site > LeftStop)
   {
      if ((*C).Basis1().total_dimension() < SInfo.MaxStates)
         this->ExpandLeftBond();
      this->IterateLeft();
   }

   this->EvolveLeftmostSite();

   while (Site < RightStop)
      this->IterateRight();
}

void
TDVP::IterateLeft2()
{
   // Form two-site centre block.
   --Site;

   LinearWavefunction::iterator CPrev = C;
   --C;
   StateComponent C2 = local_tensor_prod(*C, *CPrev);

   BasicTriangularMPO::const_iterator HPrev = H;
   --H;
   OperatorComponent H2 = local_tensor_prod(*H, *HPrev);

   HamL.pop_back();

   // Evolve two-site centre block.
   int Iter = MaxIter;
   double Err = ErrTol;

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, 0.5*Timestep, Err);

   // Perform SVD on new C2.
   AMatSVD SL(C2, Tensor::ProductBasis<BasisList, BasisList>((*C).LocalBasis(), (*CPrev).LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfo, Info);
   RealDiagonalOperator D;
   SL.ConstructMatrices(SL.begin(), Cutoff, *C, D, *CPrev);
   *C = prod(*C, D);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Sites=" << Site << "," << Site+1
                << " C2 Iter=" << Iter
                << " Err=" << Err
                << " States=" << Info.KeptStates()
		<< " TruncErr=" << Info.TruncationError()
                << std::endl;
   }

   TruncErrSum += Info.TruncationError();

   // Update the effective Hamiltonian.
   HamR.push_front(contract_from_right(herm(*HPrev), *CPrev, HamR.front(), herm(*CPrev)));

   // Evolve the current site backwards in time.
   Iter = MaxIter;
   Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, -0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void
TDVP::EvolveLeftmostSite2()
{
   // Form two-site centre block.
   LinearWavefunction::iterator CPrev = C;
   --CPrev;
   StateComponent C2 = local_tensor_prod(*CPrev, *C);

   BasicTriangularMPO::const_iterator HPrev = H;
   --HPrev;
   OperatorComponent H2 = local_tensor_prod(*HPrev, *H);

   HamL.pop_back();

   // Evolve two-site centre block.
   int Iter = MaxIter;
   double Err = ErrTol;

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, Timestep, Err);

   // Perform SVD on new C2.
   AMatSVD SL(C2, Tensor::ProductBasis<BasisList, BasisList>((*CPrev).LocalBasis(), (*C).LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfo, Info);
   RealDiagonalOperator D;
   SL.ConstructMatrices(SL.begin(), Cutoff, *CPrev, D, *C);
   *C = prod(D, *C);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Sites=" << Site-1 << "," << Site
                << " C2 Iter=" << Iter
                << " Err=" << Err
                << " States=" << Info.KeptStates()
		<< " TruncErr=" << Info.TruncationError()
                << std::endl;
   }

   MaxStates = Info.KeptStates();
   TruncErrSum += Info.TruncationError();

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*HPrev, herm(*CPrev), HamL.back(), *CPrev));
}

void
TDVP::IterateRight2()
{
   // Evolve the current site backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, -0.5*Timestep, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }

   // Form two-site centre block.
   ++Site;

   LinearWavefunction::iterator CPrev = C;
   ++C;
   StateComponent C2 = local_tensor_prod(*CPrev, *C);

   BasicTriangularMPO::const_iterator HPrev = H;
   ++H;
   OperatorComponent H2 = local_tensor_prod(*HPrev, *H);

   HamR.pop_front();

   // Evolve two-site centre block.
   Iter = MaxIter;
   Err = ErrTol;

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, 0.5*Timestep, Err);

   // Perform SVD on new C2.
   AMatSVD SL(C2, Tensor::ProductBasis<BasisList, BasisList>((*CPrev).LocalBasis(), (*C).LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfo, Info);
   RealDiagonalOperator D;
   SL.ConstructMatrices(SL.begin(), Cutoff, *CPrev, D, *C);
   *C = prod(D, *C);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Sites=" << Site-1 << "," << Site
                << " C2 Iter=" << Iter
                << " Err=" << Err
                << " States=" << Info.KeptStates()
		<< " TruncErr=" << Info.TruncationError()
                << std::endl;
   }

   MaxStates = std::max(MaxStates, Info.KeptStates());
   TruncErrSum += Info.TruncationError();

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*HPrev, herm(*CPrev), HamL.back(), *CPrev));
}

void
TDVP::Evolve2()
{
   ++TStep;
   TruncErrSum = 0.0;

   while (Site > LeftStop + 1)
      this->IterateLeft2();

   this->EvolveLeftmostSite2();

   while (Site < RightStop)
      this->IterateRight2();

   this->CalculateEps();
}

void
TDVP::CalculateEps()
{
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   // Create a local copy so we can perform a single right-to-left sweep
   // without having to go back to the right.
   LinearWavefunction PsiLocal = Psi;
   LinearWavefunction::iterator CLocal = PsiLocal.end();
   --CLocal;
   BasicTriangularMPO::const_iterator HLocal = Hamiltonian.end();
   --HLocal;
   std::deque<StateComponent> HamLLocal = HamL;
   std::deque<StateComponent> HamRLocal = HamR;
   int SiteLocal = RightStop;

   // Perform SVD to left-orthogonalize current site for NullSpace2.
   MatrixOperator M = ExpandBasis2(*CLocal);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   StateComponent CLeftOrtho = prod(*CLocal, U);
   *CLocal = prod(*CLocal, U*D*Vh);

   // Calculate error measure epsilon_1 and add to sum.
   StateComponent X = contract_from_left(*HLocal, herm(NullSpace2(CLeftOrtho)), HamLLocal.back(), *CLocal);
   double Eps1Sq = norm_frob_sq(scalar_prod(X, herm(HamRLocal.front())));
   Eps1SqSum += Eps1Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << SiteLocal
                << " Eps1Sq=" << Eps1Sq
                << std::endl;
   }

   while (SiteLocal > LeftStop)
   {
      // Perform SVD to right-orthogonalize current site.
      M = ExpandBasis1(*CLocal);

      SingularValueDecomposition(M, U, D, Vh);

      *CLocal = prod(Vh, *CLocal);

      // Calculate the right half of epsilon_2 (see below).
      StateComponent Y = contract_from_right(herm(*HLocal), NullSpace1(*CLocal), HamRLocal.front(), herm(*CLocal));

      // Update the effective Hamiltonian.
      HamRLocal.push_front(contract_from_right(herm(*HLocal), *CLocal, HamRLocal.front(), herm(*CLocal)));

      // Move to the next site.
      --SiteLocal;
      --HLocal;
      --CLocal;

      *CLocal = prod(*CLocal, U*D);

      HamLLocal.pop_back();

      // Perform SVD to left-orthogonalize current site for NullSpace2.
      M = ExpandBasis2(*CLocal);

      SingularValueDecomposition(M, U, D, Vh);

      CLeftOrtho = prod(*CLocal, U);
      *CLocal = prod(*CLocal, U*D*Vh);

      // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
      StateComponent X = contract_from_left(*HLocal, herm(NullSpace2(CLeftOrtho)), HamLLocal.back(), *CLocal);
      Eps1Sq = norm_frob_sq(scalar_prod(X, herm(HamRLocal.front())));
      double Eps2Sq = norm_frob_sq(scalar_prod(X, herm(Y)));
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
}
