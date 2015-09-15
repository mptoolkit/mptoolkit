// -*- C++ -*-

// iDMRG with single-site subspace expansion

#include "mpo/triangular_mpo.h"
#include "mps/infinitewavefunction.h"
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
#include "mp-algorithms/arnoldi.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mp-algorithms/stateslist.h"
#include "mps/operator_actions.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

#if !defined(NDEBUG)
#include "mp-algorithms/triangular_mpo_solver.h"
#endif

#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "lattice/infinite-parser.h"

namespace prog_opt = boost::program_options;

using statistics::moving_average;
using statistics::moving_exponential;

bool EarlyTermination = false;  // we set this to true if we get a checkpoint

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
      MatrixOperator Result = delta_shift(In, QShift);
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
       {
	  Result = operator_prod(herm(*I), Result, *I);
       }
      Result -= inner_prod(Proj, Result) * Ident;
      return In - Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator Proj;
   MatrixOperator Ident;
};

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
      Result -= inner_prod(Proj, Result) * Ident;
      return In - Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator const& Proj;
   MatrixOperator const& Ident;
};

std::complex<double>
MPO_EigenvaluesLeft(StateComponent& Guess, LinearWavefunction const& Psi,
		    QuantumNumber const& QShift, TriangularMPO const& Op,
		    MatrixOperator const& Rho)
{
   ProductLeft Prod(Psi, Op, QShift);
   Guess = Initial_E(Op, DeltaShift(Psi.Basis1(), adjoint(QShift)));
   MatrixOperator Ident = Guess.front();
   for (int i = 0; i < int(Guess.size())-1; ++i)
   {
      Guess.back() *= 0.0;
      Guess = Prod(Guess);
      Guess.front() = Ident;
   }
   // calculate the energy
   double Energy = inner_prod(Guess.back(), Rho).real();

   MatrixOperator H0 = Guess.back() - Energy*Guess.front();
   // Now we want the fixed point of H = U(H) + H_0
   // where U(H) is the shift one site.
   // Let F(H) = H - U(H).  Then we want the solution of F(H) = H_0

   // solve for the first component
   SubProductLeftProject ProdL(Psi, QShift, Rho, Ident);

   int m = 30;
   int max_iter = 10000;
   double tol = 1e-14;
   int Res = GmRes(Guess.back(), ProdL, H0, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>(),1);
   CHECK_EQUAL(Res, 0);

   // remove the spurious constant term from the energy
   DEBUG_TRACE("Spurious part")(inner_prod(Guess.back(), Rho));
   Guess.back() -= inner_prod(Rho, Guess.back()) * Guess.front();

#if 0
#if !defined(NDEBUG)
   KMatrixPolyType CheckEMat = SolveMPO_Left(Psi, QShift, Op, delta_shift(Ident, QShift), Rho, 1);
   ComplexPolyType EValues = ExtractOverlap(CheckEMat[std::complex<double>(1.0,0.0)], delta_shift(Rho, QShift));
   TRACE(EValues);
   TRACE(CheckEMat[std::complex<double>(1.0,0.0)]);
#endif
#endif

   return Energy;
}

std::complex<double>
MPO_EigenvaluesRight(StateComponent& Guess, LinearWavefunction const& Psi,
		     QuantumNumber const& QShift, TriangularMPO const& Op,
		     MatrixOperator const& Rho)
{
   ProductRight Prod(Psi, Op, QShift);
   Guess = Initial_F(Op, Psi.Basis2());
   MatrixOperator Ident = Guess.back();
   for (int i = 0; i < int(Guess.size())-1; ++i)
   {
      Guess.front() *= 0.0;
      Guess = Prod(Guess);
      Guess.back() = Ident;
   }
   // calculate the energy
   double Energy = inner_prod(Guess.front(), Rho).real();

   MatrixOperator H0 = Guess.front() - Energy*Guess.back();

   // Now we want the fixed point of H = U(H) + H_0
   // where U(H) is the shift one site.
   // Let F(H) = H - U(H).  Then we want the solution of F(H) = H_0

   SubProductRightProject ProdR(Psi, QShift, Rho, Ident);

   int m = 30;
   int max_iter = 10000;
   double tol = 1e-14;
   int Res = GmRes(Guess.front(), ProdR, H0, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>(),1);
   CHECK_EQUAL(Res, 0);

   // remove the spurious constant term from the energy
   Guess.front() =  Guess.front() - inner_prod(Rho, Guess.front()) * Guess.back();

   return Energy;
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

class LocalEigensolver
{
   public:
      LocalEigensolver();

      void SetInitialFidelity(int UnitCellSize, double f);

      // Apply the solver
      double Solve(StateComponent& C,
		   StateComponent const& LeftBlockHam,
		   OperatorComponent const& H,
		   StateComponent const& RightBlockHam);

      // Eigensolver parameters
      // Eigensolver tolerance is min(sqrt(AverageFidelity()) * FidelityScale, MaxTol)
      double FidelityScale;
      double MaxTol;
      double MinTol;

      int MinIter; // Minimum number of iterations to perform (unless the eigensolver breaks down)
      int MaxIter; // Stop at this number, even if the eigensolver hasn't converged

      int Verbose;

      // information on the state of the solver
      double LastEnergy() const { return LastEnergy_; }
      double LastFidelity() const { return LastFidelity_; }
      double LastTol() const { return LastTol_; }
      double LastIter() const { return LastIter_; }
      double AverageFidelity() const { return FidelityAv_.value(); }

   private:
      // information on the last solver application
      double LastFidelity_;
      double LastEnergy_;
      double LastTol_;
      int LastIter_;

      // state information
      moving_exponential<double> FidelityAv_;
};

LocalEigensolver::LocalEigensolver()
   : FidelityScale(0), MaxTol(0), MinTol(0), MinIter(0), MaxIter(0), Verbose(0)
{
}

void
LocalEigensolver::SetInitialFidelity(int UnitCellSize, double f)
{
   FidelityAv_ = moving_exponential<double>(exp(log(0.25)/UnitCellSize));
   FidelityAv_.push(f);
   TRACE(FidelityAv_.value());
}

double
LocalEigensolver::Solve(StateComponent& C,
			StateComponent const& LeftBlockHam,
			OperatorComponent const& H,
			StateComponent const& RightBlockHam)
{
   DEBUG_CHECK_EQUAL(C.Basis1(), LeftBlockHam.Basis2());
   DEBUG_CHECK_EQUAL(C.Basis2(), RightBlockHam.Basis1());
   DEBUG_CHECK_EQUAL(LeftBlockHam.LocalBasis(), H.Basis1());
   DEBUG_CHECK_EQUAL(RightBlockHam.LocalBasis(), H.Basis2());
   DEBUG_CHECK_EQUAL(C.LocalBasis(), H.LocalBasis2());

   StateComponent ROld = C;
   LastTol_ = std::min(std::sqrt(this->AverageFidelity()) * FidelityScale, MaxTol);
   LastTol_ = std::max(LastTol_, MinTol);
   //LastTol_ = std::min(this->AverageFidelity() * FidelityScale, MaxTol);
   LastIter_ = MaxIter;
   if (Verbose > 2)
   {
      std::cerr << "Starting eigensolver.  Initial guess vector has dimensions "
		<< C.Basis1().total_dimension() << " x " << C.LocalBasis().size()
		<< " x " << C.Basis2().total_dimension() << '\n';
   }
   LastEnergy_ = Lanczos(C, MPSMultiply(LeftBlockHam, H, RightBlockHam),
			 LastIter_, LastTol_, MinIter, Verbose);

   LastFidelity_ = std::max(1.0 - norm_frob(inner_prod(ROld, C)), 0.0);
   FidelityAv_.push(LastFidelity_);
   return LastEnergy_;
   
}

struct MixInfo
{
   double MixFactor;
   double RandomMixFactor;
};

// Apply subspace expansion / truncation on the left (C.Basis1()).
// Returns a matrix Lambda (not diagonal!)
// Postcondition: Lambda' C' = C (up to truncation!)
MatrixOperator
SubspaceExpandBasis1(StateComponent& C, OperatorComponent const& H, StateComponent const& RightHam,
		     MixInfo const& Mix, StatesInfo const& States, TruncationInfo& Info,
		     StateComponent const& LeftHam)
{
   // truncate - FIXME: this is the s3e step
   MatrixOperator Lambda = ExpandBasis1(C);

   MatrixOperator Rho = scalar_prod(herm(Lambda), Lambda);
   if (Mix.MixFactor > 0)
   {
      StateComponent RH = contract_from_right(herm(H), C, RightHam, herm(C));
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
      //      MatrixOperator RhoMix = operator_prod(herm(RH), Rho, RH);
      Rho += (Mix.MixFactor / trace(RhoMix)) * RhoMix;
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
   Lambda = Lambda * herm(U);
   C = prod(U, C);

   return Lambda;
}

// Apply subspace expansion / truncation on the right (C.Basis2()).
// Returns Lambda matrix (not diagonal!)
// Postcondition: C' Lambda' = C (up to truncation!)
MatrixOperator
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
      Rho += (Mix.MixFactor / trace(RhoMix)) * RhoMix;
   }
   if (Mix.RandomMixFactor > 0)
   {
      MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
      RhoMix = herm(RhoMix) * RhoMix;
      Rho += (Mix.RandomMixFactor / trace(RhoMix)) * RhoMix;
   }
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
					 States,
					 Info);
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);
   
   Lambda = U * Lambda;
   C = prod(C, herm(U));

   return Lambda;

}



class iDMRG
{
   public:
      iDMRG(LinearWavefunction const& Psi_, QuantumNumber const& QShift_, TriangularMPO const& Hamiltonian_,
	    StateComponent const& LeftHam, StateComponent const& RightHam,
	    MatrixOperator const& LambdaR,
	    std::complex<double> InitialEnergy = 0.0, int Verbose = 0);

      void SetMixInfo(MixInfo const& m);

      void SetInitialFidelity(double f)
      {
	 Solver_.SetInitialFidelity(Psi.size(), f);
      }

      LocalEigensolver& Solver() { return Solver_; }

      void SweepRight(StatesInfo const& SInfo);
      void SweepLeft(StatesInfo const& SInfo,  bool NoUpdate = false);

      // call after a SweepRight() to make Psi an infinite wavefunction
      void Finish(StatesInfo const& SInfo);

      // returns the wavefunction, which is in center-regular form
      // (all MPS orthogonalized except for one site, and Psi.Basis1() == Psi.Basis2())
      LinearWavefunction const& Wavefunction() const { return Psi; }

      // Above functions are implemented in terms of:

      // construct initial hamiltonian, called by the constructor
      void Initialize(MatrixOperator const& LambdaR, std::complex<double> InitialEnergy);

      void Solve();
      
      // at the left boundary of the unit cell, construct and save the RightHamiltonian for later use
      void SaveRightBlock(StatesInfo const& States);

      // at the right boundary of the unit cell, construct and save the LeftHamiltonian for later use
      void SaveLeftBlock(StatesInfo const& States);
      
      // at the left boundary of the unit cell, update the LeftHamiltonian and C to SaveLeftHamiltonian
      // At the end of this step, the wavefunction is in a 'regular' form with
      // C.Basis1() == Psi.Basis2()
      void UpdateLeftBlock();

      void UpdateRightBlock();

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
      MatrixOperator SaveLambda2;

      StateComponent SaveRightHamiltonian;
      MatrixOperator SaveLambda1;

      MixInfo    MixingInfo;
      TruncationInfo Info;
};

iDMRG::iDMRG(LinearWavefunction const& Psi_, QuantumNumber const& QShift_, TriangularMPO const& Hamiltonian_,
	     StateComponent const& LeftHam, StateComponent const& RightHam,
	     MatrixOperator const& LambdaR,
	     std::complex<double> InitialEnergy, int Verbose_)
   : Hamiltonian(Hamiltonian_), Psi(Psi_), QShift(QShift_), 
     LeftHamiltonian(1, LeftHam), RightHamiltonian(1, RightHam), Verbose(Verbose_)
{
   this->Initialize(LambdaR, InitialEnergy);
}

void
iDMRG::Initialize(MatrixOperator const& LambdaR, std::complex<double> InitialEnergy)
{
   // the starting point is at the right hand side, so fill out the LeftHamiltonian
   // fill out the LeftHamiltonian
   C = Psi.begin();
   H = Hamiltonian.begin();

   FirstSite = Psi.begin();
   LastSite = Psi.end();
   --LastSite;

   // Generate the left Hamiltonian matrices at each site, up to (but not including) the last site.
   while (C != LastSite)
   {
      TRACE(LeftHamiltonian.back().Basis2());
      LeftHamiltonian.push_back(contract_from_left(*H, herm(*C), LeftHamiltonian.back(), *C));
      ++C;
      ++H;
   }

   // The InitialEnergy is the energy per unit cell.  At this point,
   // the Hamiltonian we have represents the sum of interactions across 
   // one unit cell, but the boundary interactions are counted at both boundaries.
   // So this gives an energy that is (UnitCellSize+1) times the energy per site.
   // To compensate, we subtract a term off the RightHamiltonian
   // so that the initial energy is correct.

   if (InitialEnergy != 0.0)
   {
      if (Verbose > 2)
	 std::cerr << "Correcting initial energy to " << InitialEnergy << '\n';

      StateComponent R = operator_prod_inner(*H, LeftHamiltonian.back(), *C, herm(RightHamiltonian.front()));
      std::complex<double> E = inner_prod(*C, R) / inner_prod(*C, *C);
      RightHamiltonian.front().front() += (InitialEnergy - E) * RightHamiltonian.front().back();
   }

   // We are going to start by sweeping left, so we need to initialize the
   // left blocks such that UpdateLeftBlock() will function
   SaveLeftHamiltonian = LeftHamiltonian.front();
   SaveRightHamiltonian = RightHamiltonian.front();
   SaveLambda2 = LambdaR; //MatrixOperator::make_identity(Psi.Basis1());
   SaveLambda1 = delta_shift(LambdaR, QShift); // MatrixOperator::make_identity(Psi.Basis2());

   this->CheckConsistency();
}

void
iDMRG::UpdateLeftBlock()
{
   if (Verbose > 2)
   {
      std::cerr << "Updating left Hamiltonian.  Old basis has " 
		<< LeftHamiltonian.back().Basis1().total_dimension() 
		<< " states, new basis has states = " 
		<< SaveLeftHamiltonian.Basis1().total_dimension() << "\n";
   }

   LeftHamiltonian = std::deque<StateComponent>(1, delta_shift(SaveLeftHamiltonian, QShift));

   // Subtract off the energy
   LeftHamiltonian.back().back() -= Solver_.LastEnergy() * LeftHamiltonian.back().front();

   // do an SVD of SaveLambda1 so that we can invert it
   MatrixOperator U,D,Vh;
   SingularValueDecomposition(SaveLambda1, U, D, Vh);

   (*C) = prod(delta_shift(SaveLambda2, QShift) * herm(Vh) * InvertDiagonal(D) * herm(U), *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   this->CheckConsistency();
}

void
iDMRG::UpdateRightBlock()
{
   if (Verbose > 2)
   {
      std::cerr << "Updating right Hamiltonian.  Old basis has " 
		<< RightHamiltonian.front().Basis1().total_dimension() 
		<< " states, new basis has states = " 
		<< SaveRightHamiltonian.Basis1().total_dimension() << "\n";
   }

   RightHamiltonian = std::deque<StateComponent>(1, delta_shift(SaveRightHamiltonian, adjoint(QShift)));

   // Subtract off the energy
   RightHamiltonian.front().front() -= Solver_.LastEnergy() * RightHamiltonian.front().back();

   // do an SVD of SaveLambda2 so that we can invert it
   MatrixOperator U,D,Vh;
   SingularValueDecomposition(SaveLambda2, U, D, Vh);

   (*C) = prod(*C, herm(Vh) * InvertDiagonal(D) * herm(U) * delta_shift(SaveLambda1, adjoint(QShift)));

   // normalize
   *C *= 1.0 / norm_frob(*C);

   this->CheckConsistency();
}

void
iDMRG::SaveLeftBlock(StatesInfo const& States)
{
   // When we save the block, we need to end up with
   // C.Basis2() == SaveLeftHamiltonian.Basis()
   CHECK(C == LastSite);
   StateComponent L = *C;
   SaveLambda2 = SubspaceExpandBasis2(L, *H, LeftHamiltonian.back(),
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
   SaveLambda1 = SubspaceExpandBasis1(R, *H, RightHamiltonian.front(),
				      MixingInfo, States, Info, LeftHamiltonian.back());
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
   MatrixOperator Lambda = SubspaceExpandBasis1(*C, *H, RightHamiltonian.front(), MixingInfo, States, Info,
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

   *C = prod(*C, Lambda);

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
   MatrixOperator Lambda = SubspaceExpandBasis2(*C, *H, LeftHamiltonian.back(), MixingInfo, States, Info,
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

   *C = prod(Lambda, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   this->CheckConsistency();
}

void
iDMRG::SweepRight(StatesInfo const& States)
{
   this->UpdateLeftBlock();
   this->ShowInfo('P');
   this->Solve();
   this->SaveRightBlock(States);

   while (C != LastSite)
   {
      this->TruncateAndShiftRight(States);
      this->ShowInfo('R');
      this->Solve();
   }
}

void
iDMRG::SweepLeft(StatesInfo const& States, bool NoUpdate)
{
   if (!NoUpdate)
      this->UpdateRightBlock();
   this->ShowInfo('Q');
   this->Solve();
   this->SaveLeftBlock(States);

   while (C != FirstSite)
   {
      this->TruncateAndShiftLeft(States);
      this->ShowInfo('L');
      this->Solve();
   }
}

void
iDMRG::Finish(StatesInfo const& States)
{
   CHECK(C == LastSite);

   // The final truncation.
   // This is actually quite important to get a translationally invariant wavefunction
   MatrixOperator Lambda = SubspaceExpandBasis2(*C, *H, LeftHamiltonian.back(), MixingInfo, States, Info,
						RightHamiltonian.front());
   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }

   this->ShowInfo('F');

   MatrixOperator U,D,Vh;
   SingularValueDecomposition(SaveLambda2, U, D, Vh);
   (*C) = prod(*C, Lambda * herm(Vh) * InvertDiagonal(D) * herm(U));

   CHECK_EQUAL(Psi.Basis1(), DeltaShift(Psi.Basis2(), QShift));
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
	     << " Energy=" << Solver_.LastEnergy()
	     << " States=" << Info.KeptStates()
	     << " TruncError=" << Info.TruncationError()
	     << " Entropy=" << Info.KeptEntropy()
	     << " Fidelity=" << Solver_.LastFidelity()
      //			 << " FidelityAv=" << FidelityAv.value()
	     << " Iter=" << Solver_.LastIter()
	     << " Tol=" << Solver_.LastTol()
	     << '\n';
}

int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true);
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
      double EigenCutoff = -1;
      std::string FName;
      std::string HamStr;
      std::string CouplingFile;
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
      bool DoRandom = false; // true if we want to start an iteration from a random centre matrix
      std::string TargetState;
      std::vector<std::string> BoundaryState;
      double EvolveDelta = 0.0;
      double InitialFidelity = 1E-7;
      double MaxTol = 4E-4;  // never use an eigensolver tolerance larger than this
      double MinTol = 1E-16; // lower bound for the eigensolver tolerance - seems we dont really need it

      pvalue_ptr<InfiniteWavefunction> PsiPtr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian, of the form lattice:operator")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
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
	 ("mix-factor,f", prog_opt::value(&MixFactor),
	  FormatDefault("Mixing coefficient for the density matrix", MixFactor).c_str())
	 ("random-mix-factor", prog_opt::value(&RandomMixFactor),
	  FormatDefault("Random mixing for the density matrix", RandomMixFactor).c_str())
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
	 ("bootstrap,b", prog_opt::bool_switch(&NoFixedPoint),
	  "boostrap iterations by starting from a single unit cell, "
	  "instead of obtaining the fixed point Hamiltonian "
	  "('bootstrap' is necessary if the wavefunction is not orthonormal)")
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
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
	  ;

      prog_opt::options_description opt;
      opt.add(desc);
      
      prog_opt::variables_map vm;
      prog_opt::positional_options_description p;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || HamStr.empty())
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-idmrg [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX)
         : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      if (vm.count("one-site"))
	 TwoSite = !OneSite;

      bool StartFromFixedPoint = !NoFixedPoint; // we've reversed the option

      // Initialize the filesystem
      InfiniteWavefunction Psi;

      if (ExactDiag || Create)
      {
	 pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize());
      }
      else
      {
	 long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
	 PsiPtr = pheap::OpenPersistent(FName, CacheSize);
	 Psi = *PsiPtr;
      }

      // Hamiltonian
      TriangularMPO HamMPO;
      InfiniteLattice Lattice;
      boost::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);
      int const UnitCellSize = Lattice.GetUnitCell().size();
      if (WavefuncUnitCellSize == 0)
	 WavefuncUnitCellSize = UnitCellSize;

      optimize(HamMPO);

      // load the wavefunction
      if (ExactDiag)
      {
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

	 QuantumNumbers::QuantumNumber q(HamMPO[0].GetSymmetryList(), TargetState);
	 std::cout << "Target quantum number = " << q << '\n';

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

	 //         CHECK_EQUAL(num_transform_targets(q, BoundaryQ), 1)
	 //            ("The boundary quantum number is incompatible with the target quantum number");

         QuantumNumbers::QuantumNumberList LeftBoundary;
	 for (unsigned i = 0; i < BoundaryQ.size(); ++i)
	 {
	    LeftBoundary.push_back(transform_targets(q, BoundaryQ[i])[0]);
	 }

	 LinearWavefunction W;
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
	 W.push_back(ConstructFromLeftBasis(FullBL[0], B1));
	 for (int i = 1; i < WavefuncUnitCellSize; ++i)
	 {
	    W.push_back(ConstructFromLeftBasis(FullBL[i], W.get_back().Basis2()));
	 }

	 Psi.Psi = W;
	 Psi.QShift = q;
	 Psi.C_old = MatrixOperator::make_identity(B1);
	 Psi.C_right = MakeRandomMatrixOperator(Psi.Psi.Basis2(), B2);
	 // adjust for periodic basis
	 StateComponent x = prod(Psi.Psi.get_back(), Psi.C_right);
	 Psi.C_right = TruncateBasis2(x); // the Basis2 is already 1-dim.  This just orthogonalizes x
	 Psi.Psi.set_back(x);
      }
      else if (Create)
      {
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

	 QuantumNumbers::QuantumNumber q(HamMPO[0].GetSymmetryList(), TargetState);
	 std::cout << "Target quantum number = " << q << '\n';

	 QuantumNumber LBoundary, RBoundary;
	 if (BoundaryState.empty())
	 {
	    RBoundary = QuantumNumber(HamMPO[0].GetSymmetryList());
	    LBoundary = q;
	 }
	 else
	 {
	    RBoundary = QuantumNumber(HamMPO[0].GetSymmetryList(), BoundaryState[0]);
	    std::cout << "Rignt boundary quantum number is " << RBoundary << '\n';
	    if (BoundaryState.size() > 1)
	    {
	       std::cout << "WARNING: ignoring addititional boundary quantum numbers in random wavefunction\n";
	    }
	    QuantumNumbers::QuantumNumberList QL = transform_targets(q, RBoundary);
	    if (QL.size() > 1)
	    {
	       PANIC("Don't know how to handle non-scalar non-abelian target state")(RBoundary)(q);
	    }
	    LBoundary = QL[0];
	    std::cout << "Left boundary quantum number is " << LBoundary << '\n';
	 }
	 LinearWavefunction W = CreateRandomWavefunction(FullBL, LBoundary, 3, RBoundary);
	 Psi.QShift = q;
	 Psi.C_old = MatrixOperator::make_identity(W.Basis2());
	 MatrixOperator C = MatrixOperator::make_identity(W.Basis1());
         C = left_orthogonalize(C, W);
         Psi.Psi = W;
	 Psi.C_right = Psi.C_old;
	 Psi.C_old = delta_shift(Psi.C_old, q);
      }

      WavefuncUnitCellSize = Psi.Psi.size();
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
      local_basis_compatible_or_abort(Psi.Psi, HamMPO);
      
      // Get the initial Hamiltonian matrix elements
      LinearWavefunction Lin = Psi.Psi; // get_orthogonal_wavefunction(Psi);
      //      StateComponent BlockHamL = Initial_E(HamMPO.front() , Lin.Basis2());
      StateComponent BlockHamL = Initial_E(HamMPO , Psi.C_right.Basis2());
      if (StartFromFixedPoint)
      {
         MatrixOperator Rho = scalar_prod(Psi.C_right, herm(Psi.C_right));
	 //MatrixOperator Rho = scalar_prod(Psi.C_old, herm(Psi.C_old));

         //TRACE(norm_frob_sq(SubProductLeft(Lin, Psi.QShift)(MatrixOperator::make_identity(Rho.Basis1()))));
         //TRACE(norm_frob_sq(SubProductRight(Lin, Psi.QShift)(Rho)));

	 InitialEnergy = MPO_EigenvaluesLeft(BlockHamL, Lin, Psi.QShift, HamMPO, Rho);
	 std::cout << "Starting energy (left eigenvalue) = " << InitialEnergy << '\n';
      }

      LinearWavefunction LinR = get_right_orthogonal_wavefunction(Psi);
      StateComponent BlockHamR = Initial_F(HamMPO, LinR.Basis2());
      if (StartFromFixedPoint)
      {
         MatrixOperator Rho = scalar_prod(herm(Psi.C_right), Psi.C_right);
	 //MatrixOperator Rho = scalar_prod(herm(Psi.C_old), Psi.C_old);

         //TRACE(norm_frob_sq(SubProductLeft(LinR, Psi.QShift)(Rho)));

	 std::complex<double> Energy = MPO_EigenvaluesRight(BlockHamR, LinR, Psi.QShift, HamMPO, Rho);
	 std::cout << "Starting energy (right eigenvalue) = " << Energy << '\n';
      }

      // initialization complete

      Lin.set_back(prod(Lin.get_back(), Psi.C_right));

      // Construct the iDMRG object
      iDMRG idmrg(Lin, Psi.QShift, HamMPO, delta_shift(BlockHamL,Psi.QShift), 
		  BlockHamR, Psi.C_right, InitialEnergy, Verbose);
      
      idmrg.MixingInfo.MixFactor = MixFactor;
      idmrg.MixingInfo.RandomMixFactor = RandomMixFactor;
      idmrg.SetInitialFidelity(InitialFidelity);
      idmrg.Solver().MaxTol = MaxTol;
      idmrg.Solver().MinTol = MinTol;
      idmrg.Solver().MinIter = MinIter;
      idmrg.Solver().MaxIter = NumIter;
      idmrg.Solver().FidelityScale = FidelityScale;
      idmrg.Solver().Verbose = Verbose;

      int ReturnCode = 0;

      try
      {
	 bool First = true;
	 for (int i = 0; i < MyStates.size(); ++i)
	 {
	    SInfo.MaxStates = MyStates[i].NumStates;

	    if (i % 2 == 0)
	    {
	       idmrg.SweepLeft(SInfo, First);
	       First = false;
	    }
	    else
	    {
	       idmrg.SweepRight(SInfo);
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
      InfiniteWavefunction iPsi;
      iPsi.QShift = idmrg.QShift;
      iPsi.Psi = idmrg.Wavefunction();
      CHECK_EQUAL(iPsi.Psi.Basis1(), DeltaShift(iPsi.Psi.Basis2(), iPsi.QShift));
      std::cerr << "Orthogonalizing wavefunction...\n";
      orthogonalize_linear(iPsi);

      PsiPtr = new InfiniteWavefunction(iPsi);
      pheap::ShutdownPersistent(PsiPtr);

      ProcControl::Shutdown();
      return ReturnCode;
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
