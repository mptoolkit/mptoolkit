// -*- C++ -*- $Id$

// variant of iDMRG where we keep intact the unit cell.
// This prohibits relfection symmetry.

#include "mpo/triangular_operator.h"
#include "mps/infinitewavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"
#include "common/statistics.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"

#include "models/spin-su2.h"
#include "models/spin-u1.h"
#include "models/spin-z2.h"
#include "models/spin.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"
#include "models/hubbard-u1u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/bosehubbard-spinless.h"
#include "models/bosehubbard-spinless-u1.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

namespace prog_opt = boost::program_options;

using statistics::moving_average;
using statistics::moving_exponential;

// some global parameters

// MaxTol is the maximum acceptable value for the Tol parameter in the eigensolver.
// The actual tolerance is min(MaxTol, FidelityScale * sqrt(Fidelity))
double MaxTol = 4E-4;
double FidelityScale = 0.1;

bool EarlyTermination = false;  // we set this to true if we get a checkpoint

struct ProductLeft
{
   typedef StateComponent result_type;
   typedef StateComponent argument_type;

   ProductLeft(LinearWavefunction const& Psi_, TriangularOperator const& Op_,
	       QuantumNumber const& QShift_)
      : Psi(Psi_), Op(Op_), QShift(QShift_)
   {
   }

   StateComponent operator()(StateComponent const& In) const
   {
      StateComponent Guess = delta_shift(In, QShift);
      TriangularOperator::const_iterator HI = Op.begin();
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I, ++HI)
      {
	 Guess = operator_prod(herm(*HI), herm(*I), Guess, *I);
      }
      DEBUG_CHECK(HI == Op.end());
      return Guess; //delta_shift(Guess, adjoint(QShift));
   }

   LinearWavefunction const& Psi;
   TriangularOperator const& Op;
   QuantumNumber QShift;
};

struct ProductRight
{
   typedef StateComponent result_type;
   typedef StateComponent argument_type;

   ProductRight(LinearWavefunction const& Psi_, TriangularOperator const& Op_,
		QuantumNumber const& QShift_)
      : Psi(Psi_), Op(Op_), QShift(QShift_)
   {
   }

   StateComponent operator()(StateComponent const& In) const
   {
      StateComponent Guess = In;
      LinearWavefunction::const_iterator I = Psi.end();
      TriangularOperator::const_iterator HI = Op.end();
      while (I != Psi.begin())
      {
	 --I;
	 --HI;
	 Guess = operator_prod(*HI, *I, Guess, herm(*I));
      }
      DEBUG_CHECK(HI == Op.begin());
      return delta_shift(Guess, adjoint(QShift));
   }

   LinearWavefunction const& Psi;
   TriangularOperator const& Op;
   QuantumNumber QShift;
};

struct FrontProductLeft
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   FrontProductLeft(LinearWavefunction const& Psi_, TriangularOperator const& Op_,
		    StateComponent const& E_, double Energy_)
      : Psi(Psi_), Op(Op_), E(E_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      StateComponent Guess = E;
      Guess.front() = In;
      TriangularOperator::const_iterator OpIter = Op.begin();
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I, ++OpIter)
      {
	 Guess = operator_prod(herm(*OpIter), herm(*I), Guess, *I);
      }
      return Guess.front() - Energy * Guess.back();
   }

   LinearWavefunction const& Psi;
   TriangularOperator const& Op;
   StateComponent const& E;
   double Energy;
};

struct SubProductLeft
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductLeft(LinearWavefunction const& Psi_, QuantumNumber const& QShift_)
      : Psi(Psi_), QShift(QShift_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = delta_shift(In, QShift);
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
       {
	  Result = operator_prod(herm(*I), Result, *I);
       }
      return In - Result;
   }

   LinearWavefunction const& Psi;
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
      Result -= inner_prod(Result, Proj) * Ident;
      return In - Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator Proj;
   MatrixOperator Ident;
};

struct SubProductRight
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductRight(LinearWavefunction const& Psi_, QuantumNumber const& QShift_)
      : Psi(Psi_), QShift(QShift_)
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
      return In - delta_shift(Result, adjoint(QShift));
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
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
      Result -= inner_prod(Result, Proj) * Ident;
      return In - Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator const& Proj;
   MatrixOperator const& Ident;
};

std::complex<double>
MPO_EigenvaluesLeft(StateComponent& Guess, LinearWavefunction const& Psi,
		    QuantumNumber const& QShift, TriangularOperator const& Op,
		    MatrixOperator const& Rho)
{
   ProductLeft Prod(Psi, Op, QShift);
   Guess = Initial_E(Op, DeltaShift(Psi.Basis1(), adjoint(QShift)));
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

   // solve for the first component
   SubProductLeftProject ProdL(Psi, QShift, Rho, Ident);

   int m = 30;
   int max_iter = 1000;
   double tol = 1e-14;
   GmRes(Guess.front(), ProdL, H0, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());

   // remove the spurious constant term from the energy
   DEBUG_TRACE("Spurious part")(inner_prod(Guess.front(), Rho));
   Guess.front() -= inner_prod(Guess.front(), Rho) * Guess.back();

   return Energy;
}

std::complex<double>
MPO_EigenvaluesRight(StateComponent& Guess, LinearWavefunction const& Psi,
		     QuantumNumber const& QShift, TriangularOperator const& Op,
		     MatrixOperator const& Rho)
{
   ProductRight Prod(Psi, Op, QShift);
   Guess = Initial_F(Op, Psi.Basis2());
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

   SubProductRightProject ProdR(Psi, QShift, Rho, Ident);

   int m = 30;
   int max_iter = 1000;
   double tol = 1e-14;
   GmRes(Guess.back(), ProdR, H0, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());

   // remove the spurious constant term from the energy
   Guess.back() =  Guess.back() - inner_prod(Guess.back(), Rho) * Guess.front();

   return Energy;
}

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(StateComponent const& Left_,
		      StateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }

   StateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(StateComponent const& Left_,
				       StateComponent const& Right_)
   : Left(Left_), Right(Right_)
{
}

bool ExpandL = true, ExpandR = true;

//
// Sweep to the left.  Given the wavefunction, the center matrix C_r that sits
// at the right hand edge, the queue of left block Hamiltonian operators
// and the right block Hamiltonian that sits on the right of C_r, sweep to the left.
// This removes Psi.size() components from LeftBlockHam and adds them to the RightBlockHam.
// On entry, the LeftBlockHam must be at least Psi.size() elements long,
// and RightBlockHam must be at least 1 element. The wavefunction is assumed to be already
// obtained at the current position; the first action is to rotate the wavefunction to the left.
// Returns the new center matrix for the left hand position, which has the basis expanded.
MatrixOperator
DoDMRGSweepLeft(LinearWavefunction& Psi,
		MatrixOperator const& C_r,
		TriangularOperator const& Ham,
		std::deque<StateComponent>& LeftBlockHam,
		std::deque<StateComponent>& RightBlockHam,
		StatesInfo const& SInfo, int MinIter, int NumIter,
		moving_exponential<double>& FidelityAv, bool TwoSite, double MixFactor,
                double EvolveDelta)
{
   LinearWavefunction Result;
   LinearWavefunction::const_iterator I = Psi.end();
   TriangularOperator::const_iterator H = Ham.end();
      --I; --H;

   StateComponent R = prod(*I, C_r);
   MatrixOperator C = ExpandBasis1Used(R, *H);
   RightBlockHam.push_front(operator_prod(*H, R, RightBlockHam.front(), herm(R)));
   LeftBlockHam.pop_back();
   while (I != Psi.begin())
   {
      // Expand the left matrix
      --I;
      --H;
      StateComponent L = *I;
      if (TwoSite)
      {
	 L = prod(L, C);
	 C = ExpandBasis2Used(L, *H);

	 LeftBlockHam.pop_back();
	 LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), LeftBlockHam.back(), L));
      }

      // apply the solver
      int Iterations = NumIter;
      double Tol = std::min(std::sqrt(FidelityAv.value()) * FidelityScale, MaxTol);
      MatrixOperator COld = C;
      double COldNorm = norm_frob(COld);
      double Energy;
      if (EvolveDelta == 0.0)
      {
         Energy = Lanczos(C, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
                          Iterations,
                          Tol, MinIter);
      }
      else
      {
         MatrixOperator CNew = operator_prod(LeftBlockHam.back(), C, herm(RightBlockHam.front()));
         Energy = inner_prod(C, CNew).real() / COldNorm;
         C = C - EvolveDelta * CNew; // imaginary time evolution step
         C *= 1.0 / norm_frob(C);
         Iterations = 1;
      }

      double Fidelity = std::max(1.0 - norm_frob(inner_prod(COld, C)) / COldNorm, 0.0);
      FidelityAv.push(Fidelity);

      // truncate
      MatrixOperator Rho = scalar_prod(herm(C), C);
      if (MixFactor > 0)
      {
	 MatrixOperator RhoMix = operator_prod(herm(RightBlockHam.front()), Rho, RightBlockHam.front());
	 Rho += (MixFactor / trace(RhoMix)) * RhoMix;
      }
      DensityMatrix<MatrixOperator> DM(Rho);
      TruncationInfo Info;
      DensityMatrix<MatrixOperator>::const_iterator DMPivot =
         TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
                                            SInfo,
                                            Info);
      MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);

      std::cout << "L Energy=" << Energy
		<< " States=" << Info.KeptStates()
		<< " TruncError=" << Info.TruncationError()
		<< " Entropy=" << Info.KeptEntropy()
                << " Fidelity=" << Fidelity
         //			 << " FidelityAv=" << FidelityAv.value()
		<< " Iter=" << Iterations
		<< " Tol=" << Tol
                << '\n';

      C = C * herm(U);
      R = prod(U, R);
      RightBlockHam.front() = triple_prod(U, RightBlockHam.front(), herm(U));

      // shift left
      Result.push_front(R);
      R = prod(L, C);
      C = ExpandBasis1Used(R, *H);
      RightBlockHam.push_front(operator_prod(*H, R, RightBlockHam.front(), herm(R)));
      LeftBlockHam.pop_back();
   }

   // cleanup
   Result.push_front(R);
   Psi = Result;
   return C;
}

MatrixOperator
DoDMRGSweepRight(MatrixOperator const& C_l,
		 LinearWavefunction& Psi,
		 TriangularOperator const& Ham,
		 std::deque<StateComponent>& LeftBlockHam,
		 std::deque<StateComponent>& RightBlockHam,
		 StatesInfo const& SInfo, int MinIter, int NumIter,
		 moving_exponential<double>& FidelityAv, bool TwoSite, double MixFactor,
                 double EvolveDelta)
{
   LinearWavefunction Result;

   LinearWavefunction::const_iterator I = Psi.begin();
   TriangularOperator::const_iterator H = Ham.begin();

   StateComponent L = prod(C_l, *I);
   MatrixOperator C = ExpandBasis2Used(L, *H);
   LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), LeftBlockHam.back(), L));
   RightBlockHam.pop_front();

   ++I; ++H;

   while (I != Psi.end())
   {
      StateComponent R = *I;
      if (TwoSite)
      {
	 R = prod(C, R);
	 C = ExpandBasis1Used(R, *H);
	 RightBlockHam.pop_front();
	 RightBlockHam.push_front(operator_prod(*H, R, RightBlockHam.front(), herm(R)));
      }

      // apply the solver
      int Iterations = NumIter;
      double Tol = std::min(std::sqrt(FidelityAv.value()) * FidelityScale, MaxTol);
      MatrixOperator COld = C;
      double COldNorm = norm_frob(COld);
      double Energy;
      if (EvolveDelta == 0.0)
      {
         Energy = Lanczos(C, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
                          Iterations,
                          Tol, MinIter);
      }
      else
      {
         MatrixOperator CNew = operator_prod(LeftBlockHam.back(), C, herm(RightBlockHam.front()));
         Energy = inner_prod(C, CNew).real() / COldNorm;
         C = C - EvolveDelta * CNew; // imaginary time evolution step
         C *= 1.0 / norm_frob(C);
         Iterations = 1;
      }

      double Fidelity = std::max(1.0 - norm_frob(inner_prod(COld, C)) / COldNorm, 0.0);
      FidelityAv.push(Fidelity);

      // truncate
      MatrixOperator Rho = scalar_prod(C, herm(C));
      if (MixFactor > 0)
      {
	 MatrixOperator RhoMix = operator_prod(LeftBlockHam.back(), Rho, herm(LeftBlockHam.back()));
	 Rho += (MixFactor / trace(RhoMix)) * RhoMix;
      }
      DensityMatrix<MatrixOperator> DM(Rho);
      TruncationInfo Info;
      DensityMatrix<MatrixOperator>::const_iterator DMPivot =
         TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
                                            SInfo,
                                            Info);
      MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);

      std::cout << "R Energy=" << Energy
		<< " States=" << Info.KeptStates()
		<< " TruncError=" << Info.TruncationError()
		<< " Entropy=" << Info.KeptEntropy()
                << " Fidelity=" << Fidelity
         //			 << " FidelityAv=" << FidelityAv.value()
		<< " Iter=" << Iterations
		<< " Tol=" << Tol
                << '\n';

      C = U * C;
      L = prod(L, herm(U));
      LeftBlockHam.back() = triple_prod(U, LeftBlockHam.back(), herm(U));

      // shift right
      Result.push_back(L);
      L = prod(C, R);
      C = ExpandBasis2Used(L, *H);
      LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), LeftBlockHam.back(), L));
      RightBlockHam.pop_front();

      ++I;
      ++H;
   }

   // cleanup
   Result.push_back(L);
   Psi = Result;
   return C;
}

int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true);
   try
   {
      int NumIter = 20;
      int MinIter = 4;
      int MinStates = 1;
      int MaxStates = 100000;
      //double MixFactor = 0.01;
      //bool TwoSite = false;
      int NumSteps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      //bool NoVariance = false;
      //bool UseDGKS = false;
      std::string FName;
      std::string HamStr;
      double Lambda = 1.0;
      double J = 1.0;
      double J2 = 0.0;
      double B = 0;
      double D = 0;
      double U = 0;
      double Jz = 0.0;
      double t = 1.0;
      double tprime = 1.0;
      double delta = 0.0;
      double Theta = 0.0;
      double Beta = 0.0;
      double mu = 0.0;
      half_int Spin = 0.5;
      double Jleg = 1.0;
      double Jcross = 1.0;
      int NLegs = 1;
      int NMax = 3;  // maximum number of particles for Bose-Hubbard mode
      bool TwoSite = true;
      bool OneSite = false;
      double MixFactor = 0.0;
      bool NoFixedPoint = false;
      int Verbose = 0;
      bool NoOrthogonalize = false;
      bool Create = false;
      bool ExactDiag = false;
      int UnitCellSize;
#if defined(ENABLE_ONE_SITE_SCHEME)
      bool UseOneSiteScheme = false;
#endif
      std::string TargetState;
      std::string BoundaryState;
      double EvolveDelta = 0.0;
      double InitialFidelity = 1E-7;
      int KagomeUnitCell = 24;
      int Length = 100;
      //      bool TwoSiteTurn = true;  // we can control whether we want two sites at the turning points separately

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian.  Valid choices: itf, itf-z2, xxx-su2, xxx-u1, xxx, tj-zigzag-u1su2, "
          "tj-zigzag-u1, sf-zigzag-u1, klm-u1su2, klm-u1, bh, bh2, bh-u1, bh2-u1, kagome-su2")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
	 ("two-site,2", prog_opt::bool_switch(&TwoSite), "Modify two sites at once (default)")
	 ("one-site,1", prog_opt::bool_switch(&OneSite), "Modify one site at a time")
#if defined(ENABLE_ONE_SITE_SCHEME)
         ("onesiteboundary", prog_opt::bool_switch(&UseOneSiteScheme), "Modify one site also at the boundary")
#endif
	 ("max-states,m", prog_opt::value<int>(&MaxStates),
	  FormatDefault("Maximum number of states to keep", MaxStates).c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
	  FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
	 ("mix-factor,f", prog_opt::value(&MixFactor),
	  FormatDefault("Mixing coefficient for the density matrix", MixFactor).c_str())
         ("evolve", prog_opt::value(&EvolveDelta),
          "Instead of Lanczos, do imaginary time evolution with this timestep")
	 ("random,a", prog_opt::bool_switch(&Create),
	  "Create a new wavefunction starting from a random state")
	 ("exactdiag,e", prog_opt::bool_switch(&ExactDiag),
	  "Start from an effective exact diagonalization of the unit cell")
	 ("unitcell,u", prog_opt::value(&UnitCellSize),
	  "Only if --create is specified, the size of the unit cell")
	 ("target,q", prog_opt::value(&TargetState),
	  "the target quantum number per unit cell")
	 ("boundary", prog_opt::value(&BoundaryState),
	  "use this boundary quantum number for initializing the unit cell "
          "(useful for integer spin chains)")
	 ("bootstrap,b", prog_opt::bool_switch(&NoFixedPoint),
	  "boostrap iterations by starting from a single unit cell, "
	  "instead of obtaining the fixed point Hamiltonian "
	  "('bootstrap' is necessary if the wavefunction is not orthonormal)")
	 ("steps,s", prog_opt::value<int>(&NumSteps),
	  FormatDefault("Number of DMRG steps to perform", NumSteps).c_str())
	 ("no-orthogonalize", prog_opt::bool_switch(&NoOrthogonalize),
	  "Don't orthogonalize the wavefunction before saving")
         ("length", prog_opt::value(&Length), "Length of the boundary segment")
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
	 ("spin", prog_opt::value(&Spin),
	  FormatDefault("spin (for xxx,xxz,xyz hamiltonians)", Spin).c_str())
	 ("J", prog_opt::value(&J),
	  FormatDefault("nearest-neighbor exchange J (for xxx,itf, etc)", J).c_str())
	 ("J2", prog_opt::value(&J2),
	  FormatDefault("next-nearest-neighbor exchange J2 (for xxx)", J2).c_str())
	 ("D", prog_opt::value(&D),
	  FormatDefault("single-ion anisotropy (for xxx-u1 and xxx)", D).c_str())
	 ("U", prog_opt::value(&U),
	  FormatDefault("coulomb repulsion", U).c_str())
	 ("B", prog_opt::value(&B),
	  FormatDefault("magnetic field (for xxx)", B).c_str())
	 ("Jz", prog_opt::value(&Jz),
	  FormatDefault("Jz coupling (for Kondo)", Jz).c_str())
	 ("Jleg", prog_opt::value(&Jleg),
	  FormatDefault("Jleg coupling (for Kagome strop)", Jleg).c_str())
	 ("Jcross", prog_opt::value(&Jcross),
	  FormatDefault("Jcross coupling (for Kagome strip)", Jcross).c_str())
	 ("mu", prog_opt::value(&mu),
	  FormatDefault("Chemical potential (bose-hubbard)", mu).c_str())
	 ("kagome-cell", prog_opt::value(&KagomeUnitCell),
	  FormatDefault("Unit cell for kagome with plaquette (for Kagome strip with field, kagome-field-su2)", 
                        KagomeUnitCell).c_str())
	 ("nlegs", prog_opt::value(&NLegs),
	  FormatDefault("Number of legs (for triangular ladder)", NLegs).c_str())
	 ("tprime", prog_opt::value(&tprime),
	  FormatDefault("next-nearest-neighbor hopping t' (for tj-zigzag, sf-zigzag)", tprime).c_str())
         ("nmax", prog_opt::value(&NMax),
          FormatDefault("Maximum number of particles (for bose-hubbard model)", NMax).c_str())
	 ("delta", prog_opt::value(&delta),
	  FormatDefault("Zigzag ladder potential imbalance (for tj-zigzag, sf-zigzag)", delta).c_str())
	 ("theta", prog_opt::value(&Theta),
	  FormatDefault("theta (for biquadratic xxx)", Theta).c_str())
	 ("Beta", prog_opt::value(&Beta),
	  FormatDefault("Beta (for biquadratic xxx)", Beta).c_str())
	 ("lambda", prog_opt::value(&Lambda),
	  FormatDefault("transverse field strength (for itf hamiltonian)", Lambda).c_str())
         ("seed", prog_opt::value<unsigned long>(), "random seed")
	  ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || HamStr.empty())
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-idmrg [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX) 
         : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      if (vm.count("one-site"))
	 TwoSite = !OneSite;

      bool StartFromFixedPoint = !NoFixedPoint; // we've reversed the option

      // Hamiltonian
      TriangularOperator HamMPO;
      if (HamStr == "itf")
      {
	 std::cout << "Hamiltonian is transverse-field Ising, J=" << J << ", Lambda=" << Lambda << "\n";
	 SiteBlock Site = CreateSpinSite(0.5);
	 TriangularOperator Ham;
	 Ham = J * 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
	    + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
         HamMPO = Ham;
      }
      else if (HamStr == "itf-z2")
      {
	 std::cout << "Hamiltonian is transverse-field Ising with Z2, J=" << J << ", Lambda=" << Lambda << "\n";
	 SiteBlock Site = CreateZ2SpinSite(0.5);
	 TriangularOperator Ham;
	 Ham = J * 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
	    + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "xxx-su2")
      {
         if (vm.count("theta"))
         {
            J = cos(Theta * math_const::pi);
            Beta = sin(Theta * math_const::pi);
         }
	 std::cout << "Hamiltonian is XXX model with spin S=" << Spin << ", theta="<<Theta
		   << ", J=" << J << ",beta=" << Beta << ", J2=" << J2 << '\n';
	 SiteBlock Site = CreateSU2SpinSite(Spin);
	 TriangularOperator Ham;
	 Ham = J*TriangularTwoSite(-sqrt(3.0)*Site["S"], Site["S"], Site["I"].TransformsAs());
	 // The Beta*1.2 here is an SU(2) factor, because we use Q.Q instead of (S.S)^2
	 if (Beta != 0.0)
	    Ham = Ham + (Beta*1.2) * TriangularTwoSite(sqrt(5.0)*Site["Q"], Site["Q"], Site["I"].TransformsAs());
	 if (J2 != 0.0)
	    Ham = Ham + J2 * TriangularThreeSite(-sqrt(3.0)*Site["S"],
						 Site["I"], Site["S"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "xxx-u1")
      {
	 std::cout << "Hamiltonian is XXX model with spin S=" << Spin
		   << ", J=" << J << ", Jz=" << Jz << ", J2=" << J2 << ", D=" << D << ", delta=" << delta << '\n';
	 SiteBlock Site = CreateU1SpinSite(Spin);
	 TriangularOperator Ham;
	 Ham = J * (Jz * TriangularTwoSite(Site["Sz"], Site["Sz"], Site["I"].TransformsAs())
                    + 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"], Site["I"].TransformsAs())
                             + TriangularTwoSite(Site["Sm"], Site["Sp"], Site["I"].TransformsAs())));
	 // The Beta*1.2 here is an SU(2) factor, because we use Q.Q instead of (S.S)^2
	 if (J2 != 0.0)
	    Ham = Ham + J2 * (delta * TriangularThreeSite(Site["Sz"], Site["I"], Site["Sz"])
                              + 0.5 * (TriangularThreeSite(Site["Sp"], Site["I"], Site["Sm"])
                                       + TriangularThreeSite(Site["Sm"], Site["I"], Site["Sp"])));
         if (D != 0)
            Ham = Ham + D * TriangularOneSite(Site["Sz2"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "xxx")
      {
	 std::cout << "Hamiltonian is XXX model with spin S=" << Spin
		   << ", J=" << J << ", Jz=" << Jz << ", J2=" << J2 << ", D=" << D << ", B=" << B << '\n';
	 SiteBlock Site = CreateSpinSite(Spin);
	 TriangularOperator Ham;
	 Ham = (J*Jz)*TriangularTwoSite(Site["Sz"], Site["Sz"], Site["I"].TransformsAs())
            + 0.5 * J * (TriangularTwoSite(Site["Sp"], Site["Sm"], Site["I"].TransformsAs())
                                + TriangularTwoSite(Site["Sm"], Site["Sp"], Site["I"].TransformsAs()));
	 if (J2 != 0.0)
	    Ham = Ham + J2 * (delta * TriangularThreeSite(Site["Sz"], Site["I"], Site["Sz"])
                              + 0.5 * (TriangularThreeSite(Site["Sp"], Site["I"], Site["Sm"])
                                       + TriangularThreeSite(Site["Sm"], Site["I"], Site["Sp"])));
         if (D != 0)
            Ham = Ham + D * TriangularOneSite(Site["Sz2"]);
         if (B != 0)
            Ham = Ham - B * TriangularOneSite(Site["Sz"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "tj-zigzag-u1su2")
      {
	 std::cout << "Hamiltonian is t-J model with t=" << t << ", t'=" << tprime
		   << ", delta=" << delta << '\n';
	 SiteBlock Site = CreateU1SU2tJSite();
	 double tSqrt2 = (-sqrt(2.0)) * t;  // the -sqrt(2) is an SU(2) factor
	 double tprimeSqrt2 = (-sqrt(2.0)) * tprime;  // the -sqrt(2) is an SU(2) factor
	 TriangularOperator H1 = -tSqrt2 * (TriangularTwoSite(Site["CHP"], Site["C"])
					+ TriangularTwoSite(Site["CP"], Site["CH"]))
	    + (-tprimeSqrt2) * (TriangularThreeSite(Site["CHP"], Site["P"], Site["C"])
				+  TriangularThreeSite(Site["CP"], Site["P"], Site["CH"]));
         TriangularOperator H2 = H1;
         H1 += (delta / 2) * TriangularOneSite(Site["N"]);
         H2 += (-delta / 2) * TriangularOneSite(Site["N"]);
         std::vector<OperatorComponent> HamList;
	 HamList.push_back(H1[0]);
	 HamList.push_back(H2[0]);
         HamMPO = TriangularOperator(HamList);
      }
      else if (HamStr == "tj-zigzag-u1")
      {
	 std::cout << "Hamiltonian is t-J model with t=" << t << ", t'=" << tprime
		   << ", delta=" << delta << '\n';
	 SiteBlock Site = CreateU1tJSite();
	 TriangularOperator H1 = -t * (TriangularTwoSite(Site["CHupP"], Site["Cup"])
				   - TriangularTwoSite(Site["CupP"], Site["CHup"])
				   + TriangularTwoSite(Site["CHdownP"], Site["Cdown"])
				   - TriangularTwoSite(Site["CdownP"], Site["CHdown"]))
	    + (-tprime) * (TriangularThreeSite(Site["CHupP"], Site["P"], Site["Cup"])
			   -  TriangularThreeSite(Site["CupP"], Site["P"], Site["CHup"])
			   + TriangularThreeSite(Site["CHdownP"], Site["P"], Site["Cdown"])
			   -  TriangularThreeSite(Site["CdownP"], Site["P"], Site["CHdown"]));
         TriangularOperator H2 = H1;
         H1 += (delta / 2) * TriangularOneSite(Site["N"]);
         H2 += (-delta / 2) * TriangularOneSite(Site["N"]);
         std::vector<OperatorComponent> HamList;
	 HamList.push_back(H1[0]);
	 HamList.push_back(H2[0]);
         HamMPO = TriangularOperator(HamList);
      }
      else if (HamStr == "sf-zigzag-u1")
      {
	 std::cout << "Hamiltonian is spinless fermions with t=" << t << ", t'=" << tprime
		   << ", delta=" << delta << '\n';
	 SiteBlock Site = CreateU1SpinlessFermion();
	 TriangularOperator H1 = -t * (TriangularTwoSite(Site["CHP"], Site["C"])
				   - TriangularTwoSite(Site["CP"], Site["CH"]))
	    + (-tprime) * (TriangularThreeSite(Site["CHP"], Site["P"], Site["C"])
			   -  TriangularThreeSite(Site["CP"], Site["P"], Site["CH"]));
         TriangularOperator H2 = H1;
         H1 += (delta / 2) * TriangularOneSite(Site["N"]);
         H2 += (-delta / 2) * TriangularOneSite(Site["N"]);
         std::vector<OperatorComponent> HamList;
	 HamList.push_back(H1[0]);
	 HamList.push_back(H2[0]);
         HamMPO = TriangularOperator(HamList);
      }
      else if (HamStr == "klm-u1su2")
      {
	 double J = J2;
	 std::cout << "Hamiltonian is Kondo Lattice with J=" << J << '\n';
	 SiteBlock Site = CreateU1SU2KondoSite();
	 TriangularOperator Ham;
	 double tSqrt2 = (-sqrt(2.0)) * t;  // the -sqrt(2) is an SU(2) factor
         Ham = -tSqrt2 * (TriangularTwoSite(Site["CHP"], Site["C"])
                          + TriangularTwoSite(Site["CP"], Site["CH"]));
         Ham = Ham + J * TriangularOneSite(Site["ScSf"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "klm-u1")
      {
	 double J = J2;
	 std::cout << "Hamiltonian is U(1) Kondo Lattice with J=" << J << ", Jz=" << Jz << '\n';
	 SiteBlock Site = CreateU1KondoSite();
	 TriangularOperator Ham;
         Ham =  -t * (TriangularTwoSite(Site["CHupP"], Site["Cup"])
                      - TriangularTwoSite(Site["CupP"], Site["CHup"])
                      + TriangularTwoSite(Site["CHdownP"], Site["Cdown"])
                      - TriangularTwoSite(Site["CdownP"], Site["CHdown"]));
         if (J != 0)
            Ham = Ham + J * TriangularOneSite(Site["ScSf"]);
         if (Jz != 0)
            Ham = Ham + Jz * TriangularOneSite(Site["SczSfz"]);

	 HamMPO = Ham;
      }
      else if (HamStr == "hubbard-so4")
      {
	 std::cout << "Hamiltonian is Hubbard model with t=" << t << ", U = " << U << '\n';
	 SiteBlock Site = CreateSO4HubbardSiteA();
	 TriangularOperator H1 = -2.0 * t * TriangularTwoSite(Site["C_A"], Site["C_B"]);
	 TriangularOperator H2 = -2.0 * t * TriangularTwoSite(Site["C_B"], Site["C_A"]);
	 H1 = H1 + 0.25 * U * TriangularOneSite(Site["P"]);
	 H2 = H2 + 0.25 * U * TriangularOneSite(Site["P"]);
         std::vector<OperatorComponent> HamList;
	 HamList.push_back(H1[0]);
	 HamList.push_back(H2[0]);
         HamMPO = TriangularOperator(HamList);
      }
      else if (HamStr == "hubbard-u1")
      {
	 std::cout << "Hamiltonian is Hubbard model U(1) with t=" << t << ", U = " << U << '\n';
	 SiteBlock Site = CreateU1HubbardSite();
	 TriangularOperator Ham;
         Ham =  -t * (TriangularTwoSite(Site["CHupP"], Site["Cup"])
                      - TriangularTwoSite(Site["CupP"], Site["CHup"])
                      + TriangularTwoSite(Site["CHdownP"], Site["Cdown"])
                      - TriangularTwoSite(Site["CdownP"], Site["CHdown"]));
         if (U != 0)
            Ham = Ham + (U*0.25) * TriangularOneSite(Site["P"]);

         HamMPO = Ham;
      }
      else if (HamStr == "bh")
      {
         // These parameters (e.g. mu) are defined in this form to match the Lieb-Liniger model, U = m g dz / hbar^2, energy scaled by hbar^2 / 2 m dz^2
	 std::cout << "Hamiltonian is spinless Bose-Hubbard, no symmetry, T=1, U=" << U << ", mu=" << mu << "- 2.0 \n";
	 SiteBlock Site = CreateBoseHubbardSpinlessSite(3);
	 TriangularOperator Ham;
	 Ham = -1.0 * TriangularTwoSite(Site["BH"], Site["B"]) - 1.0 * TriangularTwoSite(Site["B"], Site["BH"])
	    + U * TriangularOneSite(Site["N2"]) - (mu - 2.0) * TriangularOneSite(Site["N"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "bhj")
      {
	 std::cout << "Hamiltonian is spinless Bose-Hubbard, no symmetry, J/U=" << J << ", mu=" << mu << "\n";
	 SiteBlock Site = CreateBoseHubbardSpinlessSite(NMax);
	 TriangularOperator Ham;
	 Ham = -J * TriangularTwoSite(Site["BH"], Site["B"]) - J * TriangularTwoSite(Site["B"], Site["BH"])
	    + 0.5 * TriangularOneSite(Site["N2"]) - mu * TriangularOneSite(Site["N"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "bhj-u1")
      {
	 std::cout << "Hamiltonian is spinless Bose-Hubbard, U(1), J/U=" << J << ", mu=" << mu << "\n";
	 SiteBlock Site = CreateBoseHubbardSpinlessU1Site(NMax);
	 TriangularOperator Ham;
	 Ham = -J * TriangularTwoSite(Site["BH"], Site["B"]) - J * TriangularTwoSite(Site["B"], Site["BH"])
	    + 0.5 * TriangularOneSite(Site["N2"]) - mu * TriangularOneSite(Site["N"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "bh2")
      {
         // These parameters are defined in this form to match the Lieb-Liniger model, U = m g dz / hbar^2, energy scaled by hbar^2 / 2 m dz^2
	 std::cout << "Hamiltonian is spinless Bose-Hubbard with next-nearest neighbour tunnelling, no symmetry, T=4/3, Tprime=" << tprime << "/-12, U=" << U << ", mu=" << mu << " - 2.5 \n";
	 SiteBlock Site = CreateBoseHubbardSpinlessSite(3);
	 TriangularOperator Ham;
	 Ham = -4.0/3.0 * TriangularTwoSite(Site["BH"], Site["B"]) - 4.0/3.0 * TriangularTwoSite(Site["B"], Site["BH"])
	       + 1.0/12.0*tprime * TriangularThreeSite(Site["B"], Site["I"], Site["BH"]) + 1.0/12.0*tprime * TriangularThreeSite(Site["BH"], Site["I"], Site["B"])
	    + U * TriangularOneSite(Site["N2"]) - (mu - 2.5) * TriangularOneSite(Site["N"]);
	 HamMPO = Ham;
      }
      else if (HamStr == "bh-u1")
      {
         // These parameters are defined in this form to match the Lieb-Liniger model, U = m g dz / hbar^2, energy scaled by hbar^2 / 2 m dz^2
      	 std::cout << "Hamiltonian is spinless Bose-Hubbard, U1 symmetry, T=1, U=" << U << ", Nmax=" << NMax << "\n";
      	 SiteBlock Site = CreateBoseHubbardSpinlessU1Site(NMax);
      	 TriangularOperator Ham;
      	 Ham = -1.0 * TriangularTwoSite(Site["BH"], Site["B"]) - 1.0 * TriangularTwoSite(Site["B"], Site["BH"])
      	    + U * TriangularOneSite(Site["N2"]) + 2.0 * TriangularOneSite(Site["N"]);
      	 HamMPO = Ham;
      }
      else if (HamStr == "bh2-u1")
      {
         // These parameters are defined in this form to match the Lieb-Liniger model, U = m g dz / hbar^2, energy scaled by hbar^2 / 2 m dz^2
      	 std::cout << "Hamiltonian is spinless Bose-Hubbard with next-nearest neighbour tunnelling, U1 symmetry, T=4/3, Tprime=" << tprime << "/-12, U=" << U << ", Nmax=3\n";
      	 SiteBlock Site = CreateBoseHubbardSpinlessU1Site(3);
      	 TriangularOperator Ham;
      	 Ham = -4.0/3.0 * TriangularTwoSite(Site["BH"], Site["B"]) - 4.0/3.0 * TriangularTwoSite(Site["B"], Site["BH"])
      	       + 1.0/12.0*tprime * TriangularThreeSite(Site["B"], Site["I"], Site["BH"]) + 1.0/12.0*tprime * TriangularThreeSite(Site["BH"], Site["I"], Site["B"])
      	    + U * TriangularOneSite(Site["N2"]) + 2.5 * TriangularOneSite(Site["N"]);
      	 HamMPO = Ham;
      }
      else if (HamStr == "kagome-su2")
      {
	 std::cout << "Hamiltonian is Kagome strip, Jleg=" << Jleg << ", Jcross=" << Jcross << '\n';
	 SiteBlock Site = CreateSU2SpinSite(0.5);
         // 3-site unit cell
         std::vector<BasisList> Sites(3, Site["I"].Basis().Basis());
         TriangularOperator Ham = Jcross * (TwoPointOperator(Sites, 0, Site["S"], 1, Site["S"])
                                            + TwoPointOperator(Sites, 0, Site["S"], 2, Site["S"])
                                            + TwoPointOperator(Sites, 1, Site["S"], 3, Site["S"])
                                            + TwoPointOperator(Sites, 2, Site["S"], 3, Site["S"]));
         Ham += Jleg * (TwoPointOperator(Sites, 1, Site["S"], 4, Site["S"])
                        + TwoPointOperator(Sites, 2, Site["S"], 5, Site["S"]));
         // SU(2) factor
         Ham *= -sqrt(3.0);
         HamMPO = Ham;
      }
      else if (HamStr == "kagome-field-su2")
      {
	 std::cout << "Hamiltonian is Kagome strip with field, Jleg=" << Jleg << ", Jcross=" << Jcross << '\n';
	 SiteBlock Site = CreateSU2SpinSite(0.5);
         // 24-site unit cell
         std::vector<BasisList> Sites(KagomeUnitCell, Site["I"].Basis().Basis());
         TriangularOperator Ham;
         for (int k = 0; k < KagomeUnitCell / 3; ++k)
         {
            Ham += Jcross * (TwoPointOperator(Sites, k*3 + 0, Site["S"], k*3 + 1, Site["S"])
                                            + TwoPointOperator(Sites, k*3 + 0, Site["S"], k*3 + 2, Site["S"])
                                            + TwoPointOperator(Sites, k*3 + 1, Site["S"], k*3 + 3, Site["S"])
                                            + TwoPointOperator(Sites, k*3 + 2, Site["S"], k*3 + 3, Site["S"]));
            Ham += Jleg * (TwoPointOperator(Sites, k*3 + 1, Site["S"], k*3 + 4, Site["S"])
                           + TwoPointOperator(Sites, k*3 + 2, Site["S"], k*3 + 5, Site["S"]));
         }
         // pinning fields, the diamond configuration
         Ham += TwoPointOperator(Sites, 0, Site["S"], 1, Site["S"])
            + TwoPointOperator(Sites, 1, Site["S"], 3, Site["S"])
            + TwoPointOperator(Sites, 0, Site["S"], 2, Site["S"])
            + TwoPointOperator(Sites, 2, Site["S"], 3, Site["S"]);
         // SU(2) factor
         Ham *= -sqrt(3.0);
         HamMPO = Ham;
      }
      else
      {
	 std::cerr << "mp-idmrg: error: Hamiltonian parameter must be one of itf, xxx-su2, xxx-u1, "
	    "tj-zigzag-u1su2, tj-zigzag-u1, sf-zigzag-u1, klm-u1su2, klm-u1, bh, bh2, bh-u1, bh2-u1.\n";
	 exit(1);
      }

      // load the wavefunction
      InfiniteWavefunction Psi;
      if (ExactDiag)
      {
	 pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize());
	 std::vector<BasisList> BL = ExtractLocalBasis1(HamMPO.data());
	 std::vector<BasisList> FullBL = BL;
	 while (int(FullBL.size()) < UnitCellSize)
	    std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));
	 UnitCellSize = FullBL.size();
	 std::cout << "Creating exact diagonalization basis.  Unit cell size = " << UnitCellSize << '\n';

	 QuantumNumbers::QuantumNumber q(HamMPO[0].GetSymmetryList(), TargetState);
	 std::cout << "Target quantum number = " << q << '\n';

         QuantumNumbers::QuantumNumber BoundaryQ(HamMPO[0].GetSymmetryList());
         if (BoundaryState != "")
            BoundaryQ = QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), BoundaryState);

         CHECK_EQUAL(num_transform_targets(q, BoundaryQ), 1)
            ("The boundary quantum number is incompatible with the target quantum number");

         QuantumNumbers::QuantumNumber ActualBoundary = transform_targets(q, BoundaryQ)[0];

	 LinearWavefunction W;
	 VectorBasis B1(HamMPO.front().GetSymmetryList());
	 B1.push_back(ActualBoundary, 1);
	 VectorBasis B2(HamMPO.front().GetSymmetryList());
	 B2.push_back(BoundaryQ, 1);
	 W.push_back(ConstructFromLeftBasis(FullBL[0], B1));
	 for (int i = 1; i < UnitCellSize; ++i)
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
	 std::cout << "Creating wavefunction.  unit cell size = " << UnitCellSize << '\n';
	 pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize());
	 std::vector<BasisList> BL = ExtractLocalBasis1(HamMPO.data());
	 std::vector<BasisList> FullBL = BL;
	 while (int(FullBL.size()) < UnitCellSize)
	    std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));

	 QuantumNumbers::QuantumNumber q(HamMPO[0].GetSymmetryList(), TargetState);
	 std::cout << "Target quantum number = " << q << '\n';
	 LinearWavefunction W = CreateRandomWavefunction(FullBL, q, 3);
	 Psi.QShift = q;
	 Psi.C_old = MatrixOperator::make_identity(W.Basis2());
	 MatrixOperator C = MatrixOperator::make_identity(W.Basis1());
         C = left_orthogonalize(C, W);
         Psi.Psi = W;
	 Psi.C_right = Psi.C_old;
	 Psi.C_old = delta_shift(Psi.C_old, q);
      }
      else
      {
	 long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
	 pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize);
	 Psi = *PsiPtr;
      }

      UnitCellSize = Psi.Psi.size();
      std::cout << "Unit cell size = " << UnitCellSize << '\n';

      if (vm.count("evolve"))
         std::cout << "Evolving with timestep " << EvolveDelta << '\n';

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      // replicate the HamMPO until it has the same size as the unit cell
      HamMPO = repeat(HamMPO, UnitCellSize / HamMPO.size());
      CHECK_EQUAL(int(HamMPO.size()), UnitCellSize);

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

	 std::complex<double> Energy = MPO_EigenvaluesLeft(BlockHamL, Lin, Psi.QShift, HamMPO, Rho);
	 std::cout << "Starting energy (left eigenvalue) = " << Energy << '\n';
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

      // The initial wavefunction is left-orthogonalized, so the initial center matrix
      // is at the right hand side.  Set up the block Hamiltonians
      std::deque<StateComponent> LeftBlock, RightBlock;
      LeftBlock.push_back(BlockHamL);
      RightBlock.push_back(BlockHamR);

      //      LinearWavefunction MyPsi = Psi.Psi;
      LinearWavefunction MyPsi = Lin;
      MatrixOperator C = Psi.C_right;
      QuantumNumber QShift = Psi.QShift;  // quantum number shift per unit cell

      LinearWavefunction::const_iterator I = Psi.Psi.end();
      TriangularOperator::const_iterator HI = HamMPO.begin();

      DEBUG_TRACE(inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))));
      DEBUG_TRACE(inner_prod(LeftBlock.back().front(), scalar_prod(C, herm(C))));
      DEBUG_TRACE(inner_prod(RightBlock.front().back(), scalar_prod(herm(C), C)));

      // set the initial energy to zero, so that we get the correct energy per unit cell
      // even from the first iteration
#if 1
      double LastEnergy = inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))).real();
      DEBUG_TRACE(LastEnergy);
      // subtract the constant energy from the Hamiltonian component of the LeftBlock
      LeftBlock.back().front() -= LastEnergy * LeftBlock.back().back();
#else
      double LastEnergy = 0.0;
#endif

      I = MyPsi.begin();
      HI = HamMPO.begin();
      LeftBlock.back() = delta_shift(LeftBlock.back(), QShift);
      while (I != MyPsi.end())
      {
	 LeftBlock.push_back(operator_prod(herm(*HI), herm(*I), LeftBlock.back(), *I));
	 ++HI;
	 ++I;
      }

      DEBUG_TRACE(inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))));
      DEBUG_TRACE(inner_prod(LeftBlock.back().front(), scalar_prod(C, herm(C))));
      DEBUG_TRACE(inner_prod(RightBlock.front().back(), scalar_prod(herm(C), C)));

#if 0
      I = MyPsi.begin();
      HI = HamMPO.begin();
      LeftBlock.back() = delta_shift(LeftBlock.back(), QShift);
      while (I != MyPsi.end())
      {
	 LeftBlock.push_back(operator_prod(herm(*HI), herm(*I), LeftBlock.back(), *I));
	 ++HI;
	 ++I;
      }
      TRACE(inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))));
      TRACE(inner_prod(C, operator_prod(herm(LeftBlock.back()), C, RightBlock.front())));
      TRACE(inner_prod(LeftBlock.back().front(), scalar_prod(C, herm(C))));
      TRACE(inner_prod(RightBlock.front().back(), scalar_prod(herm(C), C)));
#endif

      // initialization of the blocks

      int ReturnCode = 0;

#if defined(ENABLE_ONE_SITE_SCHEME)
      if (UseOneSiteScheme)
      {
         OneSiteScheme(Psi, MyPsi, LastEnergy, C, HamMPO, QShift, LeftBlock, RightBlock, SInfo, NumIter, MixFactor, NumSteps, Verbose);
      }
      else
      {
#endif

      StateComponent SaveLeftBlock = LeftBlock.back();
      SaveLeftBlock = delta_shift(SaveLeftBlock, QShift);

      MatrixOperator PsiL = Psi.C_old;                // ** overwriten before being used

      MatrixOperator DiagonalL = Psi.C_old;                                              // **OK**
      MatrixOperator ExpanderL = MatrixOperator::make_identity(DiagonalL.Basis2());      // **OK**

      StateComponent SaveRightBlock = RightBlock.front();  // ** overwriten before being used

      //      MatrixOperator PsiR = delta_shift(Psi.C_old, adjoint(QShift));                    // **OK**
      MatrixOperator PsiR = Psi.C_right;                    // **OK**

      MatrixOperator DiagonalR;
      MatrixOperator ExpanderR;

      MatrixOperator Vh;


      SingularValueDecomposition(C, ExpanderR, DiagonalR, Vh);
      C = ExpanderR * DiagonalR;
      RightBlock.front() = triple_prod(Vh, RightBlock.front(), herm(Vh));

      // now do the DMRG
      //int ReturnCode = 0; // return code becomes non-zero if we have a checkpoint
      int NumIterationsCompleted = 0;
      std::cout << "Starting iDMRG...\n";

      // initial energy.  If we started from the fixed point, this should be
      // the same as the energy eigenvalues
      LastEnergy = inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))).real();

      DEBUG_TRACE(LastEnergy);

      //      moving_average<double> FidelityAv(UnitCellSize*2);
      moving_exponential<double> FidelityAv(exp(log(0.25)/UnitCellSize));
      FidelityAv.push(InitialFidelity); // initialization


      // Extend the left block
      TRACE(LeftBlock.size())(MyPsi.size());

      while (LeftBlock.size() != Length+1)
      {
         LinearWavefunction::const_iterator I = MyPsi.begin();
         TriangularOperator::const_iterator HI = HamMPO.begin();
         while (I != MyPsi.end())
         {
            LeftBlock.push_back(operator_prod(herm(*HI), herm(*I), LeftBlock.back(), *I));
            ++HI;
            ++I;
         }
      }

      // now make a block of size Length, to represent the bounary
      HamMPO = repeat(HamMPO, Length / HamMPO.size());
      CHECK_EQUAL(int(HamMPO.size()), Length);

      LinearWavefunction MyPsiRpt = MyPsi;
      TRACE(LeftBlock.size());
      while (MyPsi.size() < Length)
      {
         for (LinearWavefunction::const_iterator I = MyPsiRpt.begin(); I != MyPsiRpt.end(); ++I)
         {
            MyPsi.push_back(*I);
         }
      }

      // boundary
      VectorBasis Vac = VectorBasis(make_vacuum_basis(LeftBlock.front().GetSymmetryList()));
      StateComponent BoundaryE(LeftBlock.front().LocalBasis(), Vac, Vac);
      BoundaryE.back()(0,0) = LinearAlgebra::Matrix<double>(1,1, 1.0);
      StateComponent BoundaryPsi = MyPsi.get_front();
      BoundaryPsi = StateComponent(BoundaryPsi.LocalBasis(), Vac, BoundaryPsi.Basis2());
      for (int i = 0 ; i < BoundaryPsi.size(); ++i)
      {
         BoundaryPsi[i] = MakeRandomMatrixOperator(BoundaryPsi.Basis1(), BoundaryPsi.Basis2(), BoundaryPsi.LocalBasis()[i]);
      }
#if 0
      LeftBlock.front() = BoundaryE;
      MyPsi.set_front(BoundaryPsi);
#endif


#if 1
      // open BC on the right hand side
      BoundaryE = StateComponent(RightBlock.front().LocalBasis(), Vac, Vac);
      BoundaryE.front()(0,0) = LinearAlgebra::Matrix<double>(1,1, 1.0);
      RightBlock.front() = BoundaryE;
      BoundaryPsi = MyPsi.get_back();
      BoundaryPsi = StateComponent(BoundaryPsi.LocalBasis(), BoundaryPsi.Basis1(), Vac);
      for (int i = 0 ; i < BoundaryPsi.size(); ++i)
      {
         BoundaryPsi[i] = MakeRandomMatrixOperator(BoundaryPsi.Basis1(), BoundaryPsi.Basis2(), BoundaryPsi.LocalBasis()[i]);
      }
      MyPsi.set_back(BoundaryPsi);
      C = MakeRandomMatrixOperator(Vac, Vac);
      C(0,0)(0,0) = 1.0;
#endif

      try
      {

	 for (int i = 0; i < NumSteps; ++i)
	 {

            TRACE(inner_prod(RightBlock.front().back(), scalar_prod(herm(C), C)));

	    C = DoDMRGSweepLeft(MyPsi, C, HamMPO, LeftBlock, RightBlock, SInfo, MinIter, NumIter,
				FidelityAv, TwoSite, MixFactor, EvolveDelta);

            TRACE(inner_prod(LeftBlock.back().front(), scalar_prod(C, herm(C))));

	    double Energy;
	    double Fidelity;
	    int Iterations;
	    double Tol;
	    {
	       Iterations = NumIter;
	       Tol = std::min(std::sqrt(FidelityAv.value()) * FidelityScale, MaxTol);
	       C *= 1.0 / norm_frob(C);
	       MatrixOperator COld = C;

	       //	       TRACE(C.Basis1().total_degree())(C.Basis2().total_degree());

               if (EvolveDelta == 0.0)
               {
                  Energy = Lanczos(C, SuperblockMultiply(LeftBlock.back(), RightBlock.front()),
                                   Iterations,
                                   Tol);
               }
               else
               {
                  MatrixOperator CNew = operator_prod(LeftBlock.back(), C, herm(RightBlock.front()));
                  Energy = inner_prod(C, CNew).real();
                  C = C - EvolveDelta * CNew; // imaginary time evolution step
                  C *= 1.0 / norm_frob(C);
                  Iterations = 1;
               }

	       Fidelity = std::max(1.0 - norm_frob(inner_prod(COld, C)), 0.0);
	       FidelityAv.push(Fidelity);
	    }

	    LastEnergy = Energy;

	    PsiL = C;

	    {
	       // truncate the left block
	       MatrixOperator RhoL = scalar_prod(C, herm(C));
	       if (MixFactor > 0)
	       {
		  MatrixOperator RhoMix = operator_prod(LeftBlock.back(), RhoL, herm(LeftBlock.back()));
		  RhoL += (MixFactor / trace(RhoMix)) * RhoMix;
	       }
	       DensityMatrix<MatrixOperator> DML(RhoL);
	       TruncationInfo Info;
	       MatrixOperator TruncL = DML.ConstructTruncator(DML.begin(),
							      TruncateFixTruncationErrorAbsolute(DML.begin(),
												 DML.end(),
												 SInfo,
												 Info));
	       std::cout << "A Energy=" << Energy
			 << " States=" << Info.KeptStates()
			 << " TruncError=" << Info.TruncationError()
			 << " Entropy=" << Info.KeptEntropy()
			 << " Fidelity=" << Fidelity
                  //			 << " FidelityAv=" << FidelityAv.value()
			 << " Iter=" << Iterations
			 << " Tol=" << Tol
			 << '\n';

	       C = TruncL * C;
	       LeftBlock.back() = triple_prod(TruncL, LeftBlock.back(), herm(TruncL));

	    }


	    C = DoDMRGSweepRight(C, MyPsi, HamMPO, LeftBlock, RightBlock, SInfo, MinIter, NumIter,
				 FidelityAv, TwoSite, MixFactor, EvolveDelta);

            TRACE(inner_prod(RightBlock.front().back(), scalar_prod(herm(C), C)));


	    // solve
	    {
	       Iterations = NumIter;
	       Tol = std::min(std::sqrt(FidelityAv.value()) * FidelityScale, MaxTol);
	       C *= 1.0 / norm_frob(C);
	       MatrixOperator COld = C;
               if (EvolveDelta == 0.0)
               {
                  Energy = Lanczos(C, SuperblockMultiply(LeftBlock.back(), RightBlock.front()),
                                   Iterations,
                                   Tol);
               }
               else
               {
                  MatrixOperator CNew = operator_prod(LeftBlock.back(), C, herm(RightBlock.front()));
                  Energy = inner_prod(C, CNew).real();
                  C = C - EvolveDelta * CNew; // imaginary time evolution step
                  C *= 1.0 / norm_frob(C);
                  Iterations = 1;
               }
              
	       Fidelity = std::max(1.0 - norm_frob(inner_prod(COld, C)), 0.0);
	       FidelityAv.push(Fidelity);
	    }

	    LastEnergy = Energy;
	    PsiR = C;

	    // truncate the right block
	    {
	       MatrixOperator RhoR = scalar_prod(herm(C), C);
	       if (MixFactor > 0)
	       {
		  MatrixOperator RhoMix = operator_prod(herm(RightBlock.front()), RhoR, RightBlock.front());
		  RhoR += (MixFactor / trace(RhoMix)) * RhoMix;
	       }
	       DensityMatrix<MatrixOperator> DMR(RhoR);
	       TruncationInfo Info;
	       MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(),
							      TruncateFixTruncationErrorAbsolute(DMR.begin(),
												 DMR.end(),
												 SInfo,
												 Info));
	       std::cout << "B Energy=" << Energy
			 << " States=" << Info.KeptStates()
			 << " TruncError=" << Info.TruncationError()
			 << " Entropy=" << Info.KeptEntropy()
			 << " Fidelity=" << Fidelity
                  //			 << " FidelityAv=" << FidelityAv.value()
			 << " Iter=" << Iterations
			 << " Tol=" << Tol
			 << '\n';

	       C = C * herm(TruncR);
	       //MyPsi.set_back(prod(TruncR, MyPsi.get_back()));
	       RightBlock.front() = triple_prod(TruncR, RightBlock.front(), herm(TruncR));
	    }


	    ++NumIterationsCompleted;
	    ProcControl::TestAsyncCheckpoint();
	 }

      }
      catch (ProcControl::Checkpoint& c)
      {
	 ReturnCode = c.ReturnCode();
	 std::cerr << "Early termination after " << NumIterationsCompleted << " iterations: "
		   << c.Reason() << '\n';
	 EarlyTermination = true;
      }
      catch (...)
      {
	 throw;      // if we got some other exception, don't even try and recover
      }

      // finished the iterations.  apply the truncation to the left block so that DiagonalR
      // can be the center matrix
      MatrixOperator MapToOldBasis = delta_shift(ExpanderL, adjoint(QShift));
      //      MyPsi.set_back(prod(MyPsi.get_back(), ExpanderR*herm(MapToOldBasis)));
      //MyPsi.set_front(prod(ExpanderL, MyPsi.get_front()));

      if (Verbose >= 1)
	 std::cerr << "Saving wavefunction.\n";

      // convert back to an InfiniteWavefunction
      Psi.C_old = DiagonalL;
      Psi.C_right = PsiR * herm(MapToOldBasis); //triple_prod(MapToOldBasis, DiagonalR, herm(MapToOldBasis));
      Psi.Psi = LinearWavefunction();
      for (LinearWavefunction::const_iterator I = MyPsi.begin(); I != MyPsi.end(); ++I)
      {
	 Psi.Psi.push_back(*I);
      }

      Psi.QShift = QShift;


#if defined(ENABLE_ONE_SITE_SCHEME)
      } // end test for one site scheme
#endif

      DEBUG_CHECK_EQUAL(Psi.C_old.Basis1(), Psi.C_old.Basis2());
      DEBUG_CHECK_EQUAL(Psi.C_old.Basis2(), Psi.Psi.Basis1());
      DEBUG_CHECK_EQUAL(Psi.Psi.Basis2(), Psi.C_right.Basis1());
      DEBUG_CHECK_EQUAL(Psi.Psi.Basis1(), DeltaShift(Psi.C_right.Basis2(), Psi.QShift));

      // orthogonalize it
      if (EarlyTermination && !NoOrthogonalize)
      {
	 std::cerr << "mp-idmrg: warning: early termination, not orthogonalizing the wavefunction!\n";
      }
      else if (!NoOrthogonalize)
      {
	 std::cerr << "Orthogonalizing wavefunction...\n";
	 orthogonalize(Psi);
      }

      DEBUG_CHECK_EQUAL(Psi.C_old.Basis1(), Psi.C_old.Basis2());
      DEBUG_CHECK_EQUAL(Psi.C_old.Basis2(), Psi.Psi.Basis1());
      DEBUG_CHECK_EQUAL(Psi.Psi.Basis2(), Psi.C_right.Basis1());
      DEBUG_CHECK_EQUAL(Psi.Psi.Basis1(), DeltaShift(Psi.C_right.Basis2(), Psi.QShift));

      pvalue_ptr<InfiniteWavefunction> PsiPtr = new InfiniteWavefunction(Psi);
      pheap::ShutdownPersistent(PsiPtr);

      ProcControl::Shutdown();
      return ReturnCode;
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
