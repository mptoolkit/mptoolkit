// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-idmrg4.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/wavefunc-utils.h"
#include "matrixproduct/triangularoperator.h"
#include "matrixproduct/infinitewavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"
#include "mp-algorithms/arnoldi.h"
#include "mp-algorithms/gmres.h"

#include "models/spin-su2.h"
#include "models/spin-u1.h"
#include "models/spin.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

namespace prog_opt = boost::program_options;


struct ProductLeft
{
   typedef MPStateComponent result_type;
   typedef MPStateComponent argument_type;

   ProductLeft(LinearWavefunction const& Psi_, MpOpTriangular const& Op_)
      : Psi(Psi_), Op(Op_)
   {
   }

   MPStateComponent operator()(MPStateComponent const& In) const
   {
      MPStateComponent Guess = In;
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
      {
	 Guess = operator_prod(herm(Op.data()), herm(*I), Guess, *I);
      }
      return Guess;
   }

   LinearWavefunction const& Psi;
   MpOpTriangular const& Op;
};

struct ProductRight
{
   typedef MPStateComponent result_type;
   typedef MPStateComponent argument_type;

   ProductRight(LinearWavefunction const& Psi_, MpOpTriangular const& Op_)
      : Psi(Psi_), Op(Op_)
   {
   }

   MPStateComponent operator()(MPStateComponent const& In) const
   {
      MPStateComponent Guess = In;
      LinearWavefunction::const_iterator I = Psi.end();
      while (I != Psi.begin())
      {
	 --I;
	 Guess = operator_prod(Op.data(), *I, Guess, herm(*I));
      }
      return Guess;
   }

   LinearWavefunction const& Psi;
   MpOpTriangular const& Op;
};

struct FrontProductLeft
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   FrontProductLeft(LinearWavefunction const& Psi_, MpOpTriangular const& Op_,
		    MPStateComponent const& E_, double Energy_)
      : Psi(Psi_), Op(Op_), E(E_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MPStateComponent Guess = E;
      Guess.front() = In;
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
      {
	 Guess = operator_prod(herm(Op.data()), herm(*I), Guess, *I);
      }
      return Guess.front() - Energy * Guess.back();
   }

   LinearWavefunction const& Psi;
   MpOpTriangular const& Op;
   MPStateComponent const& E;
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

std::complex<double>
MPO_EigenvaluesLeft(MPStateComponent& Guess, LinearWavefunction const& Psi, 
		    QuantumNumber const& QShift, MpOpTriangular const& Op,
		    MatrixOperator const& Rho)
{
   ProductLeft Prod(Psi, Op);
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
   SubProductLeft ProdL(Psi, QShift);

   int m = 30;
   int max_iter = 1000;
   double tol = 1e-12;
   GmRes(Guess.front(), ProdL, H0, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());

   // remove the spurious constant term from the energy
   Guess.front() =  Guess.front() - inner_prod(Guess.front(), Rho) * Guess.back();

   return Energy;
}

std::complex<double>
MPO_EigenvaluesRight(MPStateComponent& Guess, LinearWavefunction const& Psi, 
		     QuantumNumber const& QShift, MpOpTriangular const& Op,
		     MatrixOperator const& Rho)
{
   ProductRight Prod(Psi, Op);
   Guess = Initial_F(Op, Psi.Basis1());
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

   SubProductRight ProdR(Psi, QShift);

   int m = 30;
   int max_iter = 1000;
   double tol = 1e-12;
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

   SuperblockMultiply(MPStateComponent const& Left_,
		      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }

   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(MPStateComponent const& Left_,
				       MPStateComponent const& Right_)
   : Left(Left_), Right(Right_)
{
}

bool ExpandL = true, ExpandR = true;

//
// Sweep to the left.  Given the wavefunction, the center matrix C_r that sits
// at the right hand edge, the queue of left block Hamiltonian operators
// and the right block Hamiltonian that sits on the right of C_r, sweep to the left.
// This removes Psi.size() components from LeftBlockHam and adds them to the RightBlockHam.
// On exit, the leftBlockHam is overwritten by the new RightBlockHam.
MatrixOperator
DoDMRGSweepLeft(LinearWavefunction& Psi, 
		MatrixOperator const& C_r, 
		SimpleMPOperator const& Ham,
		std::deque<MPStateComponent>& LeftBlockHam,
		MPStateComponent const& IncomingHam,
		StatesInfo const& SInfo, int NumIter)
{
   LinearWavefunction Result;
   std::deque<MPStateComponent> RightBlockHam;

   LinearWavefunction::const_iterator I = Psi.end();
   SimpleMPOperator::const_iterator H = Ham.end();
   --I; --H;

   MPStateComponent R = prod(*I, C_r);
   MatrixOperator C = ExpandBasis1(R);
   RightBlockHam.push_front(IncomingHam);
   RightBlockHam.push_front(operator_prod(*H, R, IncomingHam, herm(R)));
   LeftBlockHam.pop_back();
   while (I != Psi.begin())
   {
      // Expand the left matrix
      --I;
      --H;
      MPStateComponent L = prod(*I, C);
      C = ExpandBasis2(L);

      LeftBlockHam.pop_back();
      LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), LeftBlockHam.back(), L));

      // apply the solver
      int Iterations = NumIter;
      double Energy = Lanczos(C, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
			      Iterations);

      // truncate
      MatrixOperator Rho = scalar_prod(herm(C), C);
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
		<< " Entropy=" << Info.KeptEntropy() << '\n';

      C = C * herm(U);
      R = prod(U, R);
      RightBlockHam.front() = triple_prod(U, RightBlockHam.front(), herm(U));

      // shift left
      Result.push_front(R);
      R = prod(L, C);
      C = ExpandBasis1(R);
      RightBlockHam.push_front(operator_prod(*H, R, RightBlockHam.front(), herm(R)));
      LeftBlockHam.pop_back();
   }

   // cleanup
   Result.push_front(R);
   LeftBlockHam = RightBlockHam;
   Psi = Result;
   return C;
}

MatrixOperator
DoDMRGSweepRight(MatrixOperator const& C_l, 
		 LinearWavefunction& Psi, 
		 SimpleMPOperator const& Ham,
		 MPStateComponent const& IncomingHam,
		 std::deque<MPStateComponent>& RightBlockHam,
		 StatesInfo const& SInfo, int NumIter)
{
   LinearWavefunction Result;
   std::deque<MPStateComponent> LeftBlockHam;

   LinearWavefunction::const_iterator I = Psi.begin();
   SimpleMPOperator::const_iterator H = Ham.begin();

   MPStateComponent L = prod(C_l, *I);
   MatrixOperator C = ExpandBasis2(L);

   LeftBlockHam.push_back(IncomingHam);
   LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), IncomingHam, L));
   RightBlockHam.pop_front();

   ++I; ++H;

   while (I != Psi.end())
   {
      MPStateComponent R = prod(C, *I);
      C = ExpandBasis1(R);

      RightBlockHam.pop_front();
      RightBlockHam.push_front(operator_prod(*H, R, RightBlockHam.front(), herm(R)));

      // apply the solver
      int Iterations = NumIter;
      double Energy = Lanczos(C, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
			      Iterations);

      // truncate
      MatrixOperator Rho = scalar_prod(C, herm(C));
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
		<< " Entropy=" << Info.KeptEntropy() << '\n';

      C = U * C;
      L = prod(L, herm(U));
      LeftBlockHam.back() = triple_prod(U, LeftBlockHam.back(), herm(U));

      // shift right
      Result.push_back(L);
      L = prod(C, R);
      C = ExpandBasis2(L);
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

void DoIteration(LinearWavefunction& Left, MatrixOperator& C_LR, LinearWavefunction& Right,
		 SimpleMPOperator const& LeftHam, SimpleMPOperator const& RightHam,
		 MatrixOperator& C_RL, 
		 std::deque<MPStateComponent>& LeftBlockHam,
		 std::deque<MPStateComponent>& RightBlockHam,
		 StatesInfo SInfo, int NumIter)
{
   //TRACE(C_LR);

   // These two could be done in parallel
   MPStateComponent LBack = LeftBlockHam.back();
   MatrixOperator C_left = DoDMRGSweepLeft(Left, C_LR, LeftHam, 
					   LeftBlockHam, RightBlockHam.front(), SInfo, NumIter);

   MatrixOperator C_right = DoDMRGSweepRight(C_LR, Right, RightHam, 
					     LBack, RightBlockHam, SInfo, NumIter);
   // now we have swapped LeftBlockHam and RightBlockHam
   
   // adjust for the energy per site
   {
   MatrixOperator RhoR = scalar_prod(C_right, herm(C_right));
   double RightE = inner_prod(RhoR, RightBlockHam.back().front()).real();
   // and for the left
   MatrixOperator RhoL = scalar_prod(herm(C_left), C_left);
   double LeftE = inner_prod(RhoL, LeftBlockHam.front().back()).real();
   //   TRACE(LeftE)(RightE);
   RightBlockHam.back().front() -= RightE * MatrixOperator::make_identity(RightBlockHam.back().front().Basis1());
   LeftBlockHam.front().back() -= LeftE * MatrixOperator::make_identity(LeftBlockHam.front().back().Basis1());
   }

   //TRACE(C_RL)(InvertDiagonal(C_RL));

   // update C_RL
   //TRACE(C_RL)(InvertDiagonal(C_RL, 1E-7));
   C_RL = C_right * InvertDiagonal(C_RL, 1E-7) * C_left;
   //TRACE(C_RL)(C_right)(C_left);
   
   //TRACE(C_left)(C_LR)(C_RL)(C_right)(C_LR.Basis1())(C_LR.Basis2())(C_left.Basis2())(C_right.Basis1());

   // solve

   double Energy;
   {
      int Iterations = NumIter;
      Energy = Lanczos(C_RL, SuperblockMultiply(RightBlockHam.back(), LeftBlockHam.front()),
		       Iterations);
   }

   // truncate
   {
      MatrixOperator RhoL = scalar_prod(C_RL, herm(C_RL));
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
		<< " Entropy=" << Info.KeptEntropy() << '\n';
      //DML.DensityMatrixReport(std::cout);

      C_RL = TruncL * C_RL;
      Right.set_back(prod(Right.get_back(), herm(TruncL)));
      RightBlockHam.back() = triple_prod(TruncL, RightBlockHam.back(), herm(TruncL));

      MatrixOperator RhoR = scalar_prod(herm(C_RL), C_RL);
      DensityMatrix<MatrixOperator> DMR(RhoR);
      TruncationInfo InfoR;
      MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(), 
						     TruncateFixTruncationErrorAbsolute(DMR.begin(),
											DMR.end(),
											SInfo,
											InfoR));
      C_RL = C_RL * herm(TruncR);
      Left.set_front(prod(TruncR, Left.get_front()));
      LeftBlockHam.front() = triple_prod(TruncR, LeftBlockHam.front(), herm(TruncR));
   }

   //TRACE(C_RL);

   // sweep back
   MPStateComponent RBack = RightBlockHam.back();
   C_left = DoDMRGSweepLeft(Right, C_RL, RightHam, RightBlockHam, LeftBlockHam.front(), SInfo, NumIter);
   C_right = DoDMRGSweepRight(C_RL, Left, LeftHam, RBack, LeftBlockHam, SInfo, NumIter);
   
   // adjust for the energy per site
   {
   MatrixOperator RhoL = scalar_prod(herm(C_left), C_left);
   double LeftE = inner_prod(RhoL, RightBlockHam.front().back()).real();
   // and for the left
   MatrixOperator RhoR = scalar_prod(C_right, herm(C_right));
   double RightE = inner_prod(RhoR, LeftBlockHam.back().front()).real();
   //   TRACE(LeftE)(RightE);
   LeftBlockHam.back().front() -= RightE * MatrixOperator::make_identity(LeftBlockHam.back().front().Basis1());
   RightBlockHam.front().back() -= LeftE * MatrixOperator::make_identity(RightBlockHam.front().back().Basis1());
   }

   //TRACE(C_LR)(InvertDiagonal(C_LR, 1E-7));
   C_LR = C_right * InvertDiagonal(C_LR, 1E-7) * C_left;
   //TRACE(C_LR)(C_left)(C_right);
   // update C_LR

   //   TRACE(C_left)(C_LR)(C_right)(C_LR.Basis1())(C_LR.Basis2())(C_right.Basis2())(C_left.Basis1());

   // solve
   {
      int Iterations = NumIter;
      Energy = Lanczos(C_LR, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
		       Iterations);
   }

   // truncate
   {
      MatrixOperator RhoL = scalar_prod(C_LR, herm(C_LR));
      DensityMatrix<MatrixOperator> DML(RhoL);
      TruncationInfo Info;
      MatrixOperator TruncL = DML.ConstructTruncator(DML.begin(), 
						     TruncateFixTruncationErrorAbsolute(DML.begin(),
											DML.end(),
											SInfo,
											Info));
      std::cout << "B Energy=" << Energy 
		<< " States=" << Info.KeptStates()
		<< " TruncError=" << Info.TruncationError()
		<< " Entropy=" << Info.KeptEntropy() << '\n';
      //DML.DensityMatrixReport(std::cout);

      C_LR = TruncL * C_LR;
      Left.set_back(prod(Left.get_back(), herm(TruncL)));
      LeftBlockHam.back() = triple_prod(TruncL, LeftBlockHam.back(), herm(TruncL));

      MatrixOperator RhoR = scalar_prod(herm(C_LR), C_LR);
      DensityMatrix<MatrixOperator> DMR(RhoR);
      TruncationInfo InfoR;
      MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(), 
						     TruncateFixTruncationErrorAbsolute(DMR.begin(),
											DMR.end(),
											SInfo,
											InfoR));
      C_LR = C_LR * herm(TruncR);
      Right.set_front(prod(TruncR, Right.get_front()));
      RightBlockHam.front() = triple_prod(TruncR, RightBlockHam.front(), herm(TruncR));
   }
}

template <typename T>
std::string FormatDefault(std::string const& Text, T const& Value)
{
   return Text + " [default " + boost::lexical_cast<std::string>(Value) + "]";
}

std::vector<BasisList>
ExtractLocalBasis(MpOpTriangular const& H)
{
   std::vector<BasisList> Result(1, H.LocalBasis());
   return Result;
}

int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true);
   try
   {
      int NumIter = 10;
      int MinStates = 1;
      int MaxStates = 100000;
      double MixFactor = 0.01;
      bool TwoSite = false;
      int NumSteps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      bool NoVariance = false;
      bool UseDGKS = false;
      std::string FName;
      std::string HamStr;
      double Lambda = 1.0;
      double J2 = 0.0;
      double Theta = 0.0;
      half_int Spin = 0.5;
      bool NoFixedPoint = false;
      int Verbose = 0;
      bool NoOrthogonalize = false;
      bool Create = false;
      bool ExactDiag = false;
      int UnitCellSize;
      std::string TargetState;
      bool EarlyTermination = false;  // we set this to true if we get a checkpoint
      
      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr), 
          "model Hamiltonian.  Valid choices: itf, xxx-su2, xxx-u1")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
	 ("iter,i", prog_opt::value<int>(&NumIter), 
	  FormatDefault("Number of Lanczos iterations per step", NumIter).c_str())
	 ("max-states,m", prog_opt::value<int>(&MaxStates), 
	  FormatDefault("Maximum number of states to keep", MaxStates).c_str())
         ("min-states", prog_opt::value<int>(&MinStates), 
	  FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff), 
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff), 
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
	 ("random,a", prog_opt::bool_switch(&Create),
	  "Create a new wavefunction starting from a random state")
	 ("exactdiag,e", prog_opt::bool_switch(&ExactDiag),
	  "Start from an effective exact diagonalization of the unit cell")
	 ("unitcell,u", prog_opt::value(&UnitCellSize),
	  "Only if --create is specified, the size of the unit cell")
	 ("target,q", prog_opt::value(&TargetState),
	  "Only if --create is specified, the target quantum number per unit cell")
	 ("bootstrap,b", prog_opt::bool_switch(&NoFixedPoint),
	  "boostrap iterations by starting from a single unit cell, "
	  "instead of obtaining the fixed point Hamiltonian "
	  "('bootstrap' is necessary if the wavefunction is not orthonormal)")
	 ("steps,s", prog_opt::value<int>(&NumSteps), 
	  FormatDefault("Number of DMRG steps to perform", NumSteps).c_str())
	 ("no-orthogonalize", prog_opt::bool_switch(&NoOrthogonalize),
	  "Don't orthogonalize the wavefunction before saving")
	 ("spin", prog_opt::value(&Spin), 
	  FormatDefault("spin (for xxx,xxz,xyz hamiltonians)", Spin).c_str())
	 ("J2", prog_opt::value(&J2), 
	  FormatDefault("next-nearest-neighbor hopping J2 (for xxx)", J2).c_str())
	 ("theta", prog_opt::value(&Theta), 
	  FormatDefault("theta (for xxx)", Theta).c_str())
	 ("lambda", prog_opt::value(&Lambda), 
	  FormatDefault("transverse field strength (for itf hamiltonian)", Lambda).c_str())
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

      bool StartFromFixedPoint = !NoFixedPoint; // we've reversed the option

      // Hamiltonian
      MpOpTriangular Ham;
      if (HamStr == "itf")
      {
	 std::cout << "Hamiltonian is transverse-field Ising.\n";
	 SiteBlock Site = CreateSpinSite(0.5);
	 Ham = 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
	    + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
      }
      else if (HamStr == "xxx-su2")
      {
	 double J = cos(Theta * math_const::pi);
	 double Beta = sin(Theta * math_const::pi);
	 std::cout << "Hamiltonian is XXX model with spin S=" << Spin << ", theta="<<Theta
		   << ", J=" << J << ",beta=" << Beta << ", J2=" << J2 << '\n';
	 SiteBlock Site = CreateSU2SpinSite(Spin);
	 Ham = J*TriangularTwoSite(-sqrt(3.0)*Site["S"], Site["S"], Site["I"].TransformsAs());
	 // The Beta*1.2 here is an SU(2) factor, because we use Q.Q instead of (S.S)^2
	 if (Beta != 0.0)
	    Ham = Ham + (Beta*1.2) * TriangularTwoSite(sqrt(5.0)*Site["Q"], Site["Q"], Site["I"].TransformsAs());
	 if (J2 != 0.0)
	    Ham = Ham + J2 * TriangularThreeSite(-sqrt(3.0)*Site["S"], 
						 Site["I"], Site["S"]);
      }
      else if (HamStr == "xxx-u1")
      {
	 double J = cos(Theta * math_const::pi);
	 double Beta = sin(Theta * math_const::pi);
	 std::cout << "Hamiltonian is XXX model with spin S=" << Spin << ", theta="<<Theta
		   << ", J=" << J << ", J2=" << J2 << '\n';
	 SiteBlock Site = CreateU1SpinSite(Spin);
	 Ham = J*(TriangularTwoSite(Site["Sz"], Site["Sz"], Site["I"].TransformsAs())
                  + 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"], Site["I"].TransformsAs())
                           + TriangularTwoSite(Site["Sm"], Site["Sp"], Site["I"].TransformsAs())));
	 // The Beta*1.2 here is an SU(2) factor, because we use Q.Q instead of (S.S)^2
	 if (J2 != 0.0)
	    Ham = Ham + J2 * (TriangularThreeSite(Site["Sz"], Site["I"], Site["Sz"])
                              + 0.5 * (TriangularThreeSite(Site["Sp"], Site["I"], Site["Sm"])
                                       + TriangularThreeSite(Site["Sm"], Site["I"], Site["Sp"])));
      }
      else
      {
	 std::cerr << "mp-idmrg: error: Hamiltonian parameter must be one of itf, xxx-su2, xxx-u1.\n";
	 exit(1);
      }

      // load the wavefunction
      InfiniteWavefunction Psi;
      if (ExactDiag)
      {
	 std::cout << "Creating exact diagonalization basis.  Unit cell size = " << UnitCellSize << '\n';
	 pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize());
	 std::vector<BasisList> BL = ExtractLocalBasis(Ham);
	 std::vector<BasisList> FullBL = BL;
	 while (int(FullBL.size()) < UnitCellSize)
	    std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));

	 QuantumNumbers::QuantumNumber q(Ham.GetSymmetryList(), TargetState);
	 std::cout << "Target quantum number = " << q << '\n';

	 CenterWavefunction W;
	 int Sz = UnitCellSize / 2;
	 VectorBasis B1(Ham.GetSymmetryList());
	 B1.push_back(q, 1);
	 VectorBasis B2(Ham.GetSymmetryList());
	 B2.push_back(QuantumNumber(Ham.GetSymmetryList()), 1);
	 W.PushLeft(ConstructFromLeftBasis(FullBL.front(), B1));
	 for (int i = 1; i < Sz; ++i)
	 {
	    W.PushLeft(ConstructFromLeftBasis(FullBL[i], W.Left().Basis2()));
	 }
	 W.PushRight(ConstructFromRightBasis(FullBL.back(), B2));
	 for (int i = FullBL.size()-2; i >= Sz; --i)
	 {
	    W.PushRight(ConstructFromRightBasis(FullBL[i], W.Right().Basis1()));
	 }
	 W.Center() = MakeRandomMatrixOperator(W.Left().Basis2(), W.Right().Basis1());
	 while (W.RightSize() > 1)
	    W.RotateRight();
	 W.Center() = W.Center() * ExpandBasis1(W.Right());
	 while (W.LeftSize() > 1)
	    W.RotateLeft();
	 Psi.QShift = q;
	 Psi.C_old = MatrixOperator::make_identity(W.RightVacuumBasis());
	 MatrixOperator C = MatrixOperator::make_identity(B1);
	 Psi.Psi = inject_left_old_interface(C, W.AsLinearWavefunction());

	 Psi.C_right = Psi.C_old;
      }
      else if (Create)
      {
	 std::cout << "Creating wavefunction.  unit cell size = " << UnitCellSize << '\n';
	 pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize());
	 std::vector<BasisList> BL = ExtractLocalBasis(Ham);
	 std::vector<BasisList> FullBL = BL;
	 while (int(FullBL.size()) < UnitCellSize)
	    std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));

	 QuantumNumbers::QuantumNumber q(Ham.GetSymmetryList(), TargetState);
	 std::cout << "Target quantum number = " << q << '\n';
	 MPWavefunction W = CreateRandomWavefunction(FullBL, q, 3);
	 Psi.QShift = q;
	 Psi.C_old = MatrixOperator::make_identity(W.Basis2());
	 MatrixOperator C = MatrixOperator::make_identity(W.Basis1());
	 Psi.Psi = inject_left_old_interface(C, W);
	 Psi.C_right = Psi.C_old;
      }
      else
      {
	 long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
	 pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize);
	 Psi = *PsiPtr;
      }

      UnitCellSize = Psi.Psi.size();
      std::cout << "Unit cell size = " << UnitCellSize << '\n';

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      SimpleMPOperator LeftHam, RightHam;
      for (int i = 0; i < UnitCellSize; i += 2)
      {
	 LeftHam.push_back(Ham.data());
	 RightHam.push_front(Ham.data());
      }

      // Get the initial Hamiltonian matrix elements
      LinearWavefunction Lin = get_orthogonal_wavefunction(Psi);
      MPStateComponent BlockHamL = Initial_E(Ham, Lin.Basis1());
      if (StartFromFixedPoint)
      {
	 MatrixOperator Rho = scalar_prod(Psi.C_old, herm(Psi.C_old));
	 std::complex<double> Energy = MPO_EigenvaluesLeft(BlockHamL, Lin, Psi.QShift, Ham, Rho);
	 std::cout << "Starting energy (left eigenvalue) = " << Energy << '\n';
      }

      LinearWavefunction LinR = get_right_orthogonal_wavefunction(Psi);
      MPStateComponent BlockHamR = Initial_F(Ham, LinR.Basis2());
      if (StartFromFixedPoint)
      {
	 MatrixOperator Rho = scalar_prod(Psi.C_old, herm(Psi.C_old));
	 std::complex<double> Energy = MPO_EigenvaluesRight(BlockHamR, LinR, Psi.QShift, Ham, Rho);
	 std::cout << "Starting energy (right eigenvalue) = " << Energy << '\n';
      }

      // setup the wavefunction by splitting the unit cell
      LinearWavefunction Left, Right;
      MatrixOperator C_RL = Psi.C_old;
      MatrixOperator C_LR = Psi.C_right;
      LinearWavefunction::const_iterator I = Psi.Psi.end();
      for (unsigned i = 0; i < Psi.Psi.size()/2; ++i)
      {
	 --I;
	 MPStateComponent A = prod(*I, C_LR);
	 C_LR = TruncateBasis1(A);
	 Right.push_front(A);
      }
      while (I != Psi.Psi.begin())
      {
	 --I;
	 Left.push_front(*I);
      }

      // block hamiltonian
      std::deque<MPStateComponent> LeftBlockHam, RightBlockHam;
      LeftBlockHam.push_back(BlockHamL);
      I = Left.begin();
      SimpleMPOperator::const_iterator HI = LeftHam.begin();
      while (I != Left.end())
      {
	 LeftBlockHam.push_back(operator_prod(herm(*HI), herm(*I), LeftBlockHam.back(), *I));
	 ++HI;
	 ++I;
      }

      RightBlockHam.push_front(BlockHamR);
      I = Right.end();
      HI = RightHam.end();
      while (I != Right.begin())
      {
	 --HI;
	 --I;
	 RightBlockHam.push_front(operator_prod(*HI, *I, RightBlockHam.front(), herm(*I)));
      }

      // now do the DMRG
      int ReturnCode = 0; // return code becomes non-zero if we have a checkpoint
      int NumIterationsCompleted = 0;
      std::cout << "Starting iDMRG...\n";
      try
      {

	 {
	    // initial energy
	    int Iterations = NumIter;
	    double Energy = Lanczos(C_LR, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
				    Iterations);
	    std::cerr << "Initial energy = " << Energy << '\n';
	 }

	 for (int i = 0; i < NumSteps; ++i)
	 {
	    DoIteration(Left, C_LR, Right, LeftHam, RightHam, C_RL, 
			LeftBlockHam, RightBlockHam, SInfo, NumIter);
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

      if (Verbose >= 1)
	 std::cerr << "Saving wavefunction.\n";

      // convert back to an InfiniteWavefunction
      Psi.C_old = C_RL;
      Psi.Psi = LinearWavefunction();
      for (LinearWavefunction::const_iterator I = Left.begin(); I != Left.end(); ++I)
      {
	 Psi.Psi.push_back(*I);
      }

      MatrixOperator U = C_LR;
      for (LinearWavefunction::const_iterator I = Right.begin(); I != Right.end(); ++I)
      {
	 MPStateComponent x = prod(U, *I);
	 U = TruncateBasis2(x);
	 Psi.Psi.push_back(x);
      }
      Psi.C_right = U;
      Psi.QShift = QuantumNumbers::QuantumNumber(U.GetSymmetryList());

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
