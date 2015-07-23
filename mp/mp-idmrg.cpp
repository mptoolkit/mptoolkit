// -*- C++ -*- $Id$ 

// variant of iDMRG where we keep intact the unit cell.
// This prohibits relfection symmetry.

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

// some global parameters

// MaxTol is the maximum acceptable value for the Tol parameter in the eigensolver.
// The actual tolerance is min(MaxTol, FidelityScale * sqrt(Fidelity))
double MaxTol = 4E-4;
double FidelityScale = 0.1;
int Verbose = 0;

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
	 Guess = operator_prod(herm(*HI), herm(*I), Guess, *I);
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
	 Guess = operator_prod(*HI, *I, Guess, herm(*I));
      }
      DEBUG_CHECK(HI == Op.begin());
      return delta_shift(Guess, adjoint(QShift));
   }

   LinearWavefunction const& Psi;
   TriangularMPO const& Op;
   QuantumNumber QShift;
};

struct FrontProductLeft
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   FrontProductLeft(LinearWavefunction const& Psi_, TriangularMPO const& Op_,
		    StateComponent const& E_, double Energy_)
      : Psi(Psi_), Op(Op_), E(E_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      StateComponent Guess = E;
      Guess.back() = In;
      TriangularMPO::const_iterator OpIter = Op.begin();
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I, ++OpIter)
      {
	 Guess = operator_prod(herm(*OpIter), herm(*I), Guess, *I);
      }
      return Guess.back() - Energy * Guess.front();
   }

   LinearWavefunction const& Psi;
   TriangularMPO const& Op;
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
      Result -= inner_prod(Proj, Result) * Ident;
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

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(StateComponent const& Left_,
		      StateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      //     TRACE(Left.front())(Right.back());
      return operator_prod_regular(Left, Psi, herm(Right));
      //return operator_prod(Left, Psi, herm(Right));
   }

   StateComponent const& Left;
   StateComponent const& Right;
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
		TriangularMPO const& Ham,
		std::deque<StateComponent>& LeftBlockHam,
		std::deque<StateComponent>& RightBlockHam,
		StatesInfo const& SInfo, int MinIter, int NumIter,
		moving_exponential<double>& FidelityAv, bool TwoSite, double MixFactor,
                double RandomMixFactor,
                double EvolveDelta)
{
   LinearWavefunction Result;
   LinearWavefunction::const_iterator I = Psi.end();
   TriangularMPO::const_iterator H = Ham.end();
      --I; --H;

   StateComponent R = prod(*I, C_r);
   MatrixOperator C = ExpandBasis1Used(R, *H);
   //MatrixOperator C = ExpandBasis1(R);
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
	 //C = ExpandBasis2(L);

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
                          Tol, MinIter, Verbose);
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
      if (RandomMixFactor > 0)
      {
	 MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
         RhoMix = herm(RhoMix) * RhoMix;
	 Rho += (RandomMixFactor / trace(RhoMix)) * RhoMix;
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
      //C = ExpandBasis1(R);
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
		 TriangularMPO const& Ham,
		 std::deque<StateComponent>& LeftBlockHam,
		 std::deque<StateComponent>& RightBlockHam,
		 StatesInfo const& SInfo, int MinIter, int NumIter,
		 moving_exponential<double>& FidelityAv, bool TwoSite, double MixFactor,
                 double RandomMixFactor,
                 double EvolveDelta)
{
   LinearWavefunction Result;

   LinearWavefunction::const_iterator I = Psi.begin();
   TriangularMPO::const_iterator H = Ham.begin();

   StateComponent L = prod(C_l, *I);
   MatrixOperator C = ExpandBasis2Used(L, *H);
   //MatrixOperator C = ExpandBasis2(L);
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
	 //C = ExpandBasis1(R);
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
                          Tol, MinIter, Verbose);
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
      if (RandomMixFactor > 0)
      {
	 MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
         RhoMix = herm(RhoMix) * RhoMix;
	 Rho += (RandomMixFactor / trace(RhoMix)) * RhoMix;
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
      //C = ExpandBasis2(L);
      LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), LeftBlockHam.back(), L));
      RightBlockHam.pop_front();

      ++I;
      ++H;
   }

   // cleanup
   Result.push_back(L);
   RightBlockHam = LeftBlockHam;
   Psi = Result;
   return C;
}

#if defined(ENABLE_ONE_SITE_SCHEME)
void OneSiteScheme(InfiniteWavefunction& Psi, LinearWavefunction& MyPsi, double& LastEnergy, MatrixOperator& C,
                   TriangularMPO const& HamMPO,
                   QuantumNumber const& QShift,
                   std::deque<StateComponent>& LeftBlock,
                   std::deque<StateComponent>& RightBlock,
                   StatesInfo const& SInfo, int MinIter, int NumIter,
                   double MixFactor, int NumSteps, bool Verbose)

{
   bool TwoSite = false;

   // initialization of the blocks

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
   int ReturnCode = 0; // return code becomes non-zero if we have a checkpoint
   int NumIterationsCompleted = 0;
   std::cout << "Starting iDMRG...\n";

   // initial energy.  If we started from the fixed point, this should be
   // the same as the energy eigenvalues
   LastEnergy = inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))).real();

   DEBUG_TRACE(LastEnergy);

   int UnitCellSize = Psi.size();
   moving_average<double> FidelityAv(UnitCellSize);
   FidelityAv.push(MaxTol); // initialization

   try
   {
      for (int i = 0; i < NumSteps; ++i)
      {
         C = DoDMRGSweepLeft(MyPsi, C, HamMPO, LeftBlock, RightBlock, SInfo, NumIter,
                             FidelityAv, TwoSite, MixFactor, RandomMixFactor);

         // now comes the slightly tricky part, where we turn around

         // retrieve the wrapped around left block from the last iteration
         LeftBlock = std::deque<StateComponent>(1, SaveLeftBlock);

         C = InvertDiagonal(DiagonalL, InverseTol) * C;
         C = herm(ExpanderL) * C;
         C = delta_shift(PsiR, QShift) * C;
         // C is now dm x dm

         DEBUG_CHECK_EQUAL(C.Basis1(), LeftBlock.back().Basis2());

         LeftBlock.back().back() -= LastEnergy * LeftBlock.back().front();

         // solve

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

            Energy = Lanczos(C, SuperblockMultiply(LeftBlock.back(), RightBlock.front()),
                             Iterations,
                             Tol, MinIter, Verbose);
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
               if (RandomMixFactor > 0)
               {
                  MatrixOperator RhoMix = MakeRandomMatrixOperator(RhoL.Basis1(), RhoL.Basis2());
                  RhoMix = herm(RhoMix) * RhoMix;
                  RhoL += (RandomMixFactor / trace(RhoMix)) * RhoMix;
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
                      << " Iter=" << Iterations
                      << " Tol=" << Tol
                      << '\n';

            C = TruncL * C;
            LeftBlock.back() = triple_prod(TruncL, LeftBlock.back(), herm(TruncL));

         }

         {
            // DiagonalL becomes the matrix of singular values in the m-dimensional truncated basis
            MatrixOperator U;
            SingularValueDecomposition(C, U, DiagonalL, ExpanderL);
            C = DiagonalL * ExpanderL;
            LeftBlock.back() = triple_prod(herm(U), LeftBlock.back(), U);
         }

         DEBUG_CHECK_EQUAL(C.Basis2(), RightBlock.front().Basis1());

         SaveRightBlock = RightBlock.front();  // the right block at the left-hand edge of the unit cell
         SaveRightBlock = delta_shift(SaveRightBlock, adjoint(QShift));

         // right-moving sweep

         C = DoDMRGSweepRight(C, MyPsi, HamMPO, LeftBlock, RightBlock, SInfo, NumIter,
                              FidelityAv, TwoSite, MixFactor, RandomMixFactor);

         // turn around at the right-hand side
         SaveLeftBlock = LeftBlock.back();
         SaveLeftBlock = delta_shift(SaveLeftBlock, QShift);

         // retrieve the wrapped-around block
         RightBlock = std::deque<StateComponent>(1, SaveRightBlock);

         C = C * InvertDiagonal(DiagonalR, InverseTol);
         C = C * herm(ExpanderR);
         C = C * delta_shift(PsiL, adjoint(QShift));

         DEBUG_CHECK_EQUAL(C.Basis2(), RightBlock.front().Basis1());

         // make the energy zero
         RightBlock.front().front() -= LastEnergy * RightBlock.front().back();

         // solve
         {
            Iterations = NumIter;
            Tol = std::min(std::sqrt(FidelityAv.value()) * FidelityScale, MaxTol);
            C *= 1.0 / norm_frob(C);
            MatrixOperator COld = C;
            Energy = Lanczos(C, SuperblockMultiply(LeftBlock.back(), RightBlock.front()),
                             Iterations,
                             Tol, MinIter, Verbose);
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
               if (RandomMixFactor > 0)
               {
                  MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
                  RhoMix = herm(RhoMix) * RhoMix;
                  Rho += (RandomMixFactor / trace(RhoMix)) * RhoMix;
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
         {
            // DiagonalR becomes the matrix of singular values in the m-dimensional truncated basis
            MatrixOperator U, Vt;
            SingularValueDecomposition(C, ExpanderR, DiagonalR, Vt);
            C = ExpanderR * DiagonalR;
            //MyPsi.set_back(prod(MyPsi.get_back(), U));
            RightBlock.front() = triple_prod(Vt, RightBlock.front(), herm(Vt));
            //LeftBlock.back() = triple_prod(herm(U), LeftBlock.back(), U);
         }

         //PsiR = C;

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
}
#endif // ENABLE_ONE_SITE_SCHEME

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
      //double MixFactor = 0.01;
      //bool TwoSite = false;
      int NumSteps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      //bool NoVariance = false;
      //bool UseDGKS = false;
      std::string FName;
      std::string HamStr;
      std::string CouplingFile;
      bool TwoSite = true;
      bool OneSite = false;
      int WavefuncUnitCellSize = 0;
      double MixFactor = 0.0;
      double RandomMixFactor = 0.0;
      bool NoFixedPoint = false;
      bool NoOrthogonalize = false;
      bool Create = false;
      bool ExactDiag = false;
#if defined(ENABLE_ONE_SITE_SCHEME)
      bool UseOneSiteScheme = false;
#endif
      bool DoRandom = false; // true if we want to start an iteration from a random centre matrix
      std::string TargetState;
      std::vector<std::string> BoundaryState;
      double EvolveDelta = 0.0;
      double InitialFidelity = 1E-7;
      //      bool TwoSiteTurn = true;  // we can control whether we want two sites at the turning points separately

      pvalue_ptr<InfiniteWavefunction> PsiPtr;

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

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
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || HamStr.empty())
      {
         print_copyright(std::cerr);
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
      TriangularMPO::const_iterator HI = HamMPO.begin();

      DEBUG_TRACE(inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))));
      DEBUG_TRACE(inner_prod(LeftBlock.back().back(), scalar_prod(C, herm(C))));
      DEBUG_TRACE(inner_prod(RightBlock.front().front(), scalar_prod(herm(C), C)));

      // set the initial energy to zero, so that we get the correct energy per unit cell
      // even from the first iteration
#if 1
      double LastEnergy = inner_prod(C, operator_prod(LeftBlock.back(), C, herm(RightBlock.front()))).real();
      DEBUG_TRACE(LastEnergy);
      // subtract the constant energy from the Hamiltonian component of the LeftBlock
      LeftBlock.back().back() -= LastEnergy * LeftBlock.back().front();
#else
      double LastEnergy = 0.0;
#endif

      DEBUG_TRACE(LeftBlock.back().front());

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
      DEBUG_TRACE(inner_prod(LeftBlock.back().back(), scalar_prod(C, herm(C))));
      DEBUG_TRACE(inner_prod(RightBlock.front().front(), scalar_prod(herm(C), C)));

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

      DEBUG_TRACE(HamMPO);

      try
      {
	 for (int i = 0; i < MyStates.size(); ++i)
	 {
	    SInfo.MaxStates = MyStates[i].NumStates;

	    C = DoDMRGSweepLeft(MyPsi, C, HamMPO, LeftBlock, RightBlock, SInfo, MinIter, NumIter,
				FidelityAv, TwoSite, MixFactor, RandomMixFactor, EvolveDelta);

	    // now comes the slightly tricky part, where we turn around

	    // retrieve the wrapped around left block from the last iteration
	    LeftBlock = std::deque<StateComponent>(1, SaveLeftBlock);

	    C = InvertDiagonal(DiagonalL, InverseTol) * C;
	    C = herm(ExpanderL) * C;
	    C = delta_shift(PsiR, QShift) * C;
	    // C is now dm x dm

            DEBUG_TRACE(trace(scalar_prod(C,herm(C))));

	    DEBUG_CHECK_EQUAL(C.Basis1(), LeftBlock.back().Basis2());

            // Subtract off the energy
	    LeftBlock.back().back() -= LastEnergy * LeftBlock.back().front();
	    //std::cout << "EShift=" << LastEnergy << '\n';

	    // solve

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
                  if (DoRandom)
                  {
                     if (Verbose)
                        std::cout << "Randomizing centre matrix.\n";
                     C = MakeRandomMatrixOperator(C.Basis1(), C.Basis2());
                     //DoRandom = false;
                  }
                  Energy = Lanczos(C, SuperblockMultiply(LeftBlock.back(), RightBlock.front()),
                                   Iterations,
                                   Tol, MinIter, Verbose);
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
               if (RandomMixFactor > 0)
               {
                  MatrixOperator RhoMix = MakeRandomMatrixOperator(RhoL.Basis1(), RhoL.Basis2());
                  RhoMix = herm(RhoMix) * RhoMix;
                  RhoL += (RandomMixFactor / trace(RhoMix)) * RhoMix;
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

	    {
	       // DiagonalL becomes the matrix of singular values in the m-dimensional truncated basis
	       MatrixOperator U;
	       SingularValueDecomposition(C, U, DiagonalL, ExpanderL);
	       C = DiagonalL * ExpanderL;
	       LeftBlock.back() = triple_prod(herm(U), LeftBlock.back(), U);
	    }

	    DEBUG_CHECK_EQUAL(C.Basis2(), RightBlock.front().Basis1());

	    SaveRightBlock = RightBlock.front();  // the right block at the left-hand edge of the unit cell
	    SaveRightBlock = delta_shift(SaveRightBlock, adjoint(QShift));

	    // right-moving sweep

	    C = DoDMRGSweepRight(C, MyPsi, HamMPO, LeftBlock, RightBlock, SInfo, MinIter, NumIter,
				 FidelityAv, TwoSite, MixFactor, RandomMixFactor, EvolveDelta);

	    // turn around at the right-hand side
	    SaveLeftBlock = LeftBlock.back();
	    SaveLeftBlock = delta_shift(SaveLeftBlock, QShift);

	    // retrieve the wrapped-around block
	    RightBlock = std::deque<StateComponent>(1, SaveRightBlock);

	    C = C * InvertDiagonal(DiagonalR, InverseTol);
	    C = C * herm(ExpanderR);
	    C = C * delta_shift(PsiL, adjoint(QShift));

            DEBUG_TRACE(trace(scalar_prod(C,herm(C))));

	    DEBUG_CHECK_EQUAL(C.Basis2(), RightBlock.front().Basis1());

	    // make the energy zero
	    RightBlock.front().front() -= LastEnergy * RightBlock.front().back();

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
                                   Tol, MinIter, Verbose);
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
               if (RandomMixFactor > 0)
               {
                  MatrixOperator RhoMix = MakeRandomMatrixOperator(RhoR.Basis1(), RhoR.Basis2());
                  RhoMix = herm(RhoMix) * RhoMix;
                  RhoR += (RandomMixFactor / trace(RhoMix)) * RhoMix;
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
	    {
	       // DiagonalR becomes the matrix of singular values in the m-dimensional truncated basis
	       MatrixOperator U, Vt;
	       SingularValueDecomposition(C, ExpanderR, DiagonalR, Vt);
	       C = ExpanderR * DiagonalR;
	       //MyPsi.set_back(prod(MyPsi.get_back(), U));
	       RightBlock.front() = triple_prod(Vt, RightBlock.front(), herm(Vt));
	       //LeftBlock.back() = triple_prod(herm(U), LeftBlock.back(), U);
	    }

	    //PsiR = C;

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

      PsiPtr = new InfiniteWavefunction(Psi);
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
