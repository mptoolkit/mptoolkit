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
#include "lanczos-exponential-new.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "linearalgebra/eigen.h"
#include "lattice/infinite-parser.h"
#include "parser/parser.h"

std::complex<double> const I(0.0, 1.0);

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

std::map<std::string, Composition>
Compositions = {
   {"secondorder", Composition(2, "Standard second-order symmetric composition", {0.5, 0.5})},
   {"triplejump4", Composition(4, "Fourth-order triple jump composition",
         {0.6756035959798289, 0.6756035959798289,
         -0.8512071919596577, -0.8512071919596577,
         0.6756035959798289, 0.6756035959798289})},
   {"suzukifractal4", Composition(4, "Fourth-order Suzuki fractal composition",
         {0.20724538589718786, 0.20724538589718786,
         0.20724538589718786, 0.20724538589718786,
         -0.3289815435887514, -0.3289815435887514,
         0.20724538589718786, 0.20724538589718786,
         0.20724538589718786, 0.20724538589718786})},
   {"mclachlan4-10", Composition(4, "Fourth-order 10-term McLachlan composition",
         {0.08926945422647525, 0.31073054577352477,
         -0.40806658841042026, 0.3080665884104203,
         0.2, 0.2,
         0.3080665884104203, -0.40806658841042026,
         0.31073054577352477, 0.08926945422647525})},
   {"symmetric4-10", Composition(4, "Symmetric fourth-order 10-term Barthel-Zhang composition",
         {0.12843317950293848, 0.12843317950293848,
         0.3388120161527937, 0.3388120161527937,
         -0.4344903913114644, -0.4344903913114644,
         0.3388120161527937, 0.3388120161527937,
         0.12843317950293848, 0.12843317950293848})},
   {"optimized4-10", Composition(4, "Optimized fourth-order 10-term Barthel-Zhang composition",
         {0.09584850274120368, 0.3306761585746725,
         -0.40878731749631037, 0.2883920480412131,
         0.19387060813922113, 0.19387060813922113,
         0.2883920480412131, -0.40878731749631037,
         0.3306761585746725, 0.09584850274120368})},
   {"lazy-10", Composition(4, "Lazy 10-term composition",
         {0.1, 0.1,
         0.1, 0.1,
         0.1, 0.1,
         0.1, 0.1,
         0.1, 0.1})}
};

Hamiltonian::Hamiltonian(std::string HamStr, int Size_,
                         std::string Magnus_, std::string TimeVar_)
   : Size(Size_), Magnus(Magnus_), TimeVar(TimeVar_)
{
   std::tie(HamOperator, Lattice) = ParseOperatorStringAndLattice(HamStr);

   // Attempt to convert the operator into an MPO.  If it works, then the
   // Hamiltonian is time-independent. If it fails, assume it failed because
   // the time variable isn't defined yet.  If there is some other error in the
   // operator, we'll catch it later.
   try
   {
      HamMPO = ParseTriangularOperator(Lattice, HamOperator);
      if (Size != 0)
         HamMPO = repeat(HamMPO, Size / HamMPO.size());
      TimeDependent = false;
   }
   catch (Parser::ParserError& e)
   {
      //if (Verbose > 1)
      //   std::cerr << "Parser error converting the Hamiltonian to an MPO - assuming the Hamiltonian is time-dependent." << std::endl;
      TimeDependent = true;
   }
   catch (...)
   {
      throw;
   }
}

BasicTriangularMPO
Hamiltonian::operator()(std::complex<double> t, std::complex<double> dt) const
{
   if (TimeDependent == false)
      return HamMPO;
   else
   {
      BasicTriangularMPO HamMPO_;

      if (Magnus == "2")
      {
         std::complex<double> Time = t + 0.5*dt;
         HamMPO_ = ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, Time}});
      }
      else if (Magnus == "4")
      {
         std::complex<double> t1 = t + (0.5 - std::sqrt(3.0)/6.0) * dt;
         std::complex<double> t2 = t + (0.5 + std::sqrt(3.0)/6.0) * dt;
         BasicTriangularMPO A1 = ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, t1}});
         BasicTriangularMPO A2 = ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, t2}});
         HamMPO_ = 0.5 * (A1 + A2) + (1.0 / 12.0) * (A1*A2 - A2*A1);
      }

      if (Size != 0)
         HamMPO_ = repeat(HamMPO_, Size / HamMPO_.size());

      return HamMPO_;
   }
}

void
Hamiltonian::set_size(int Size_)
{
   Size = Size_;

   if (Size != 0 && TimeDependent == false)
      HamMPO = repeat(HamMPO, Size / HamMPO.size());
}

TDVP::TDVP(Hamiltonian const& Ham_, std::complex<double> InitialTime_,
           std::complex<double> Timestep_, Composition Comp_, int MaxIter_,
           double ErrTol_, StatesInfo SInfo_, bool Epsilon_, int Verbose_)
   : Ham(Ham_), InitialTime(InitialTime_), Timestep(Timestep_), Comp(Comp_),
     MaxIter(MaxIter_), ErrTol(ErrTol_), SInfo(SInfo_), Epsilon(Epsilon_),
     Verbose(Verbose_)
{
}

TDVP::TDVP(FiniteWavefunctionLeft const& Psi_, Hamiltonian const& Ham_,
           std::complex<double> InitialTime_, std::complex<double> Timestep_,
           Composition Comp_, int MaxIter_, double ErrTol_, StatesInfo SInfo_,
           bool Epsilon_, int Verbose_)
   : TDVP(Ham_, InitialTime_, Timestep_, Comp_, MaxIter_, ErrTol_, SInfo_,
          Epsilon_, Verbose_)
{
   // Initialize Psi and Ham.
   Time = InitialTime;
   HamMPO = Ham(InitialTime);

   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;
   H = HamMPO.begin();
   HamL.push_back(Initial_E(HamMPO, Psi_.Basis1()));
   for (FiniteWavefunctionLeft::const_mps_iterator I = Psi_.begin(); I != Psi_.end(); ++I)
   {
      if (Verbose > 1)
         std::cout << "Site " << (HamL.size()) << std::endl;
      HamL.push_back(contract_from_left(*H, herm(*I), HamL.back(), *I));
      Psi.push_back(*I);
      MaxStates = std::max(MaxStates, (*I).Basis2().total_dimension());
      ++H;
   }
   HamR.push_front(Initial_F(HamMPO, Psi_.Basis2()));
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
TDVP::EvolveCurrentSite(std::complex<double> Tau)
{
   // Evolve current site.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, -I*Tau, Err);

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
TDVP::IterateLeft(std::complex<double> Tau)
{
   // Perform SVD to right-orthogonalize current site.
   MatrixOperator M = ExpandBasis1(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   *C = prod(Vh, *C);

   // Update the effective Hamiltonian.
   HamR.push_front(contract_from_right(herm(*H), *C, HamR.front(), herm(*C)));

   // Evolve the UD term backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;
   MatrixOperator UD = U*D;

   UD = LanczosExponential(UD, HEff2(HamL.back(), HamR.front()), Iter, I*Tau, Err);

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
TDVP::IterateRight(std::complex<double> Tau)
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

   DVh = LanczosExponential(DVh, HEff2(HamL.back(), HamR.front()), Iter, I*Tau, Err);

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
}

void
TDVP::SweepLeft(std::complex<double> Tau)
{
   HamMPO = Ham(Time, Tau);
   H = HamMPO.end();
   --H;

   while (Site > LeftStop)
   {
      this->EvolveCurrentSite(Tau);
      this->IterateLeft(Tau);
   }

   this->EvolveCurrentSite(Tau);
}

void
TDVP::SweepRight(std::complex<double> Tau)
{
   HamMPO = Ham(Time, Tau);
   H = HamMPO.begin();

   this->EvolveCurrentSite(Tau);

   while (Site < RightStop)
   {
      this->IterateRight(Tau);
      this->EvolveCurrentSite(Tau);
   }
}

void
TDVP::CalculateEps1()
{
   // Calculate error measure epsilon_1 and add to sum.
   StateComponent Y = contract_from_right(herm(*H), NullSpace1(*C), HamR.front(), herm(*C));
   double Eps1Sq = norm_frob_sq(scalar_prod(HamL.back(), herm(Y)));
   Eps1SqSum += Eps1Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " Eps1Sq=" << Eps1Sq
                << std::endl;
   }
}

void
TDVP::CalculateEps12()
{
   // Calculate the left half of epsilon_2.
   StateComponent HamLBack = HamL.back();
   HamL.pop_back();

   LinearWavefunction::const_iterator CPrev = C;
   --CPrev;
   BasicTriangularMPO::const_iterator HPrev = H;
   --HPrev;

   StateComponent X = contract_from_left(*HPrev, herm(NullSpace2(*CPrev)), HamL.back(), *CPrev);

   HamL.push_back(HamLBack);

   // Perform SVD to right-orthogonalize current site for NullSpace1.
   MatrixOperator M = ExpandBasis1(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

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
                << " Eps1Sq=" << Eps1Sq
                << " Eps2Sq=" << Eps2Sq
                << std::endl;
   }
}

void
TDVP::SweepRightFinal(std::complex<double> Tau)
{
   HamMPO = Ham(Time, Tau);
   H = HamMPO.begin();

   this->EvolveCurrentSite(Tau);
   this->CalculateEps1();

   while (Site < RightStop)
   {
      this->IterateRight(Tau);
      this->EvolveCurrentSite(Tau);
      this->CalculateEps12();
   }
}

void
TDVP::Evolve()
{
   Time = InitialTime + ((double) TStep)*Timestep;
   ++TStep;

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::vector<double>::const_iterator Gamma = Comp.Gamma.cbegin();
   std::vector<double>::const_iterator GammaEnd = Comp.Gamma.cend();
   --GammaEnd;

   this->SweepLeft((*Gamma)*Timestep);
   Time += (*Gamma)*Timestep;
   ++Gamma;

   while(Gamma != GammaEnd)
   {
      this->SweepRight((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;

      this->SweepLeft((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;
   }

   if (Epsilon)
      this->SweepRightFinal((*Gamma)*Timestep);
   else
      this->SweepRight((*Gamma)*Timestep);

   Time += (*Gamma)*Timestep;
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
TDVP::SweepLeftExpand(std::complex<double> Tau)
{
   HamMPO = Ham(Time, Tau);
   H = HamMPO.end();
   --H;

   while (Site > LeftStop)
   {
      if ((*C).Basis1().total_dimension() < SInfo.MaxStates)
         this->ExpandLeftBond();
      this->EvolveCurrentSite(Tau);
      this->IterateLeft(Tau);
   }

   this->EvolveCurrentSite(Tau);
}

void
TDVP::EvolveExpand()
{
   Time = InitialTime + ((double) TStep)*Timestep;
   ++TStep;

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::vector<double>::const_iterator Gamma = Comp.Gamma.cbegin();
   std::vector<double>::const_iterator GammaEnd = Comp.Gamma.cend();
   --GammaEnd;

   this->SweepLeftExpand((*Gamma)*Timestep);
   Time += (*Gamma)*Timestep;
   ++Gamma;

   while(Gamma != GammaEnd)
   {
      this->SweepRight((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;

      this->SweepLeftExpand((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;
   }
   
   if (Epsilon)
      this->SweepRightFinal((*Gamma)*Timestep);
   else
      this->SweepRight((*Gamma)*Timestep);

   Time += (*Gamma)*Timestep;
}

void
TDVP::IterateLeft2(std::complex<double> Tau)
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

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, -I*Tau, Err);

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

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, I*Tau, Err);

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
TDVP::EvolveLeftmostSite2(std::complex<double> Tau)
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

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, -I*Tau, Err);

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
TDVP::IterateRight2(std::complex<double> Tau)
{
   // Evolve the current site backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, I*Tau, Err);

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

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, -I*Tau, Err);

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
TDVP::SweepLeft2(std::complex<double> Tau)
{
   HamMPO = Ham(Time, Tau);
   H = HamMPO.end();
   --H;

   while (Site > LeftStop + 1)
      this->IterateLeft2(Tau);

   this->EvolveLeftmostSite2(Tau);
}

void
TDVP::SweepRight2(std::complex<double> Tau)
{
   HamMPO = Ham(Time, Tau);
   H = HamMPO.begin();
   ++H;

   this->EvolveLeftmostSite2(Tau);

   while (Site < RightStop)
      this->IterateRight2(Tau);
}

void
TDVP::Evolve2()
{
   Time = InitialTime + ((double) TStep)*Timestep;
   ++TStep;

   TruncErrSum = 0.0;

   std::vector<double>::const_iterator Gamma = Comp.Gamma.cbegin();

   while(Gamma != Comp.Gamma.cend())
   {
      this->SweepLeft2((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;

      this->SweepRight2((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;
   }
   
   if (Epsilon)
      this->CalculateEps();
}

void
TDVP::CalculateEps()
{
   HamMPO = Ham(Time);
   H = HamMPO.end();
   --H;

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   // Create a local copy so we can perform a single right-to-left sweep
   // without having to go back to the right.
   LinearWavefunction PsiLocal = Psi;
   LinearWavefunction::iterator CLocal = PsiLocal.end();
   --CLocal;
   BasicTriangularMPO::const_iterator HLocal = HamMPO.end();
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
