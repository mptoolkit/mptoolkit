// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-itebd2.cpp
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

#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/prog_opt_accum.h"
#include "interface/inittemp.h"
#include "mps/state_component.h"
#include "mps/density.h"
#include "mp-algorithms/random_wavefunc.h"
#include "models/bosehubbard-spinless-u1.h"

#include "models/spin-su2.h"
#include "models/spin.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mps/infinitewavefunction.h"
#include "tensor/tensor_exponential.h"

namespace prog_opt = boost::program_options;

void RandomizeRowSigns(MatrixOperator& TruncL)
{
   for (unsigned j = 0; j < TruncL.Basis1().size(); ++j)
   {
      for (int k = 0; k < TruncL.Basis1().dim(j); ++k)
      {
	 int sign = ((rand() / 1000) % 2) * 2 - 1;
	 for (unsigned l = 0; l < TruncL.Basis2().size(); ++l)
	 {
	    if (iterate_at(TruncL.data(), j,l))
	    {
	       for (int m = 0; m < TruncL.Basis2().dim(l); ++m)
	       {
		  TruncL(j,l)(k,m) *= sign;
	       }
	    }
	 }
      }
   }
}   

TruncationInfo
DoIteration(StateComponent& A, RealDiagonalOperator& Center, StateComponent& B,
	    SimpleOperator const& EvolOperator,
	    StatesInfo const& SInfo, SimpleOperator const& EnergyOperator, double& Energy)
{
   // if the wavefunction has ~zero weight, things might get screwy
   if (norm_frob(Center) < 1e-10)
   {
      std::cerr << "warning: Center matrix has small norm!\n";
      //      Center = MakeRandomMatrixOperator(Center.Basis1(), Center.Basis2());
      //      Center *= 1.0 / norm_frob(Center);    // normalize
   }

   // coarse-grain the tensors
   StateComponent Pair = local_tensor_prod(prod(A, Center), B);


   // bond energy prior to evolution
   double Energy1 = inner_prod(Pair, local_prod(EnergyOperator, Pair)).real() / norm_frob(Pair);

   // apply the evolution operator
   Pair = local_prod(EvolOperator, Pair);

   // bond energy after to evolution
   double Energy2 = inner_prod(Pair, local_prod(EnergyOperator, Pair)).real() / norm_frob(Pair);

   Energy = (Energy1 + Energy2) / 2;  // this is rather hueristic, but a better guess then Energy2
   
   // truncate
   SingularDecomposition<StateComponent, StateComponent> 
      SL(Pair, Tensor::make_product_basis(A.LocalBasis(), B.LocalBasis()));

   TruncationInfo Info;
   RealDiagonalOperator C;
   SL.ConstructMatrices(SL.begin(), TruncateFixTruncationErrorRelative(SL.begin(),
								       SL.end(),
								       SInfo,
								       Info),
			A, C, B);

   // normalize
   C *= 1.0 / norm_frob(C);
   Center = C;

   return Info;
}

RealDiagonalOperator ProjectRealDiagonal(MatrixOperator const& M)
{
   RealDiagonalOperator Result(M.Basis1(), M.Basis2(), M.TransformsAs());

   for (unsigned i = 0; i < M.Basis1().size(); ++i)
   {
      if (iterate_at(M.data(), i, i))
      {
	 Result(i,i) = LinearAlgebra::DiagonalMatrix<double>(M.Basis2().dim(i), M.Basis1().dim(i), 0.0);
	 for (int j = 0; j < std::min(M.Basis1().dim(i), M.Basis2().dim(i)); ++j)
	 {
            Result(i,i)(j,j) = M(i,i)(j,j).real();
	 }
      }
   }
   return Result;
}

InfiniteWavefunction MakeWavefunction(StateComponent const& A, MatrixOperator const& Lambda2,
				      StateComponent const& B, MatrixOperator const& Lambda1,
                                      QuantumNumber const& QShift)
{
   // Convert to an InfiniteWavefunction 
   InfiniteWavefunction Psi;
   Psi.C_old = Lambda1;
   Psi.Psi.push_back(A);
   StateComponent BNew = prod(Lambda2, B);
   MatrixOperator Lambda2New = TruncateBasis2(BNew);
   Psi.Psi.push_back(BNew);
   Psi.C_right = Lambda2New;
   Psi.QShift = QShift;
   return Psi;
}

void LoadFromWavefunction(InfiniteWavefunction const& Psi,
                          StateComponent& A, RealDiagonalOperator& Lambda2,
                          StateComponent& B, RealDiagonalOperator& Lambda1, QuantumNumber& QShift)
{
   Lambda1 = ProjectRealDiagonal(Psi.C_old);

   A = Psi.Psi.get_front();
   B = Psi.Psi.get_back();
   B = prod(B, Psi.C_right);
   Lambda2 = RealDiagonalOperator::make_identity(B.Basis1());
   QShift = Psi.QShift;
   //prod(InvertDiagonal(Lambda2, Epsilon), Psi.Psi.get_back());
   //   Lambda2 = ProjectRealDiagonal(Psi.C_right);
}

int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true);
   try
   {
      double DeltaTau = 0.0;
      double DeltaT = 0.0;
      double Epsilon = 1e-7;      // epsilon for inverse of singular values
      int MinStates = 1;
      int MaxStates = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      std::string FName;
      std::string HamStr;
      bool CreateWavefunction = false;
      int NumIter = 10;
      QuantumNumber QShift;

      double J = 1.0;
      double U = 0;
      double Lambda = 1.0;
      int NMax = 3;

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian.  Valid choices: itf, xxx-su2, bh-u1")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction file to evolve (required)")
	 ("max-states,m", prog_opt::value<int>(&MaxStates),
	  FormatDefault("Maximum number of states to keep", MaxStates).c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
	  FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("create,c", prog_opt::bool_switch(&CreateWavefunction),
          "Create a new (random) wavefunction for the initial state, overwriting any existing file")
         ("numiter,n", prog_opt::value(&NumIter),
          FormatDefault("Number of timesteps to perform", NumIter).c_str())
         ("time,t", prog_opt::value(&DeltaT),
          FormatDefault("real part of the timestep per iteration", DeltaT).c_str())
         ("tau", prog_opt::value(&DeltaTau),
          FormatDefault("imaginary part of the timestep per iteration", DeltaTau).c_str())
         ("epsilon", prog_opt::value(&Epsilon),
          FormatDefault("Cutoff singular value for inversion", Epsilon).c_str())
	 ("lambda", prog_opt::value(&Lambda),
	  FormatDefault("transverse field strength (for itf hamiltonian)", Lambda).c_str())
	 ("J", prog_opt::value(&J),
	  FormatDefault("Hopping (for bh hamiltonian)", Lambda).c_str())
	 ("U", prog_opt::value(&U),
	  FormatDefault("Coulomb repulsion (for bh hamiltonian)", Lambda).c_str())
         ("nmax", prog_opt::value(&NMax),
          FormatDefault("Maximum number of particles (for bose-hubbard model)", NMax).c_str())
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
         std::cerr << "usage: mp-itebd [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // we require 1 or more iterations
      if (NumIter < 1)
      {
         std::cerr << "mp-itebd: error: numiter must be > 0\n";
         exit(1);
      }

      // Set up the Hamiltonian
      SiteBlock Site;
      SimpleOperator BondHamiltonian;
      if (HamStr == "itf")
      {
	 std::cout << "Hamiltonian is transverse-field Ising, J=" << J << ", Lambda=" << Lambda << "\n";
         Site = CreateSpinSite(0.5);
         BondHamiltonian = J * 4.0 * tensor_prod(Site["Sz"], Site["Sz"])
            + Lambda * (tensor_prod(Site["Sx"], Site["I"]) + tensor_prod(Site["I"], Site["Sx"]));
      }
      else if (HamStr == "xxx-su2")
      {
         Site = CreateSU2SpinSite(0.5);
         QuantumNumbers::QuantumNumber Ident(Site.GetSymmetryList());
         BondHamiltonian = -sqrt(3.0) * tensor_prod(Site["S"], Site["S"], Ident);
      }
      else if (HamStr == "bh-u1")
      {
	 std::cout << "Hamiltonian is Bose-Hubbard, J=" << J << ", U=" << U << "\n";
         Site = CreateBoseHubbardSpinlessU1Site(NMax);
         QuantumNumbers::QuantumNumber Ident(Site.GetSymmetryList());
         BondHamiltonian = -J * (tensor_prod(Site["BH"], Site["B"]) + tensor_prod(Site["B"], Site["BH"]))
            + (U / 4.0) * (tensor_prod(Site["N2"], Site["I"]) + tensor_prod(Site["I"], Site["N2"]));
      }
      else
      {
         std::cerr << "mp-itebd: error: unrecognized Hamiltonian\n";
         exit(1);
      }

      // Construct or load the MPS
      StateComponent A,B;
      RealDiagonalOperator Lambda1, Lambda2, Lambda1Inv, Lambda2Inv;

      if (CreateWavefunction)
      {
	 std::cout << "Creating wavefunction.\n";
	 pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize());

         VectorBasis B1 = VectorBasis(make_vacuum_basis(Site.GetSymmetryList()));

         InfiniteWavefunction Psi;
         Psi.C_old = MatrixOperator::make_identity(B1);
         Psi.C_right = MatrixOperator::make_identity(B1);
         Psi.QShift = QuantumNumbers::QuantumNumber(Site.GetSymmetryList());
         std::vector<BasisList> BasisVec(2, Site.Basis1().Basis());
         Psi.Psi = CreateRandomWavefunction(BasisVec, QuantumNumbers::QuantumNumber(Site.GetSymmetryList()), 3);
         LoadFromWavefunction(Psi, A, Lambda2, B, Lambda1, QShift);
      }
      else
      {
	 long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
	 pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize);
         LoadFromWavefunction(*PsiPtr, A, Lambda2, B, Lambda1, QShift);
      }

      // Initialize the inverse matrices
      double Eps = Epsilon;
      Lambda1Inv = InvertDiagonal(Lambda1, Eps);
      Lambda2Inv = InvertDiagonal(Lambda2, Eps);
      
      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      TruncationInfo Info;
      // Initial evolution operator for a half timestep
      SimpleOperator EvolOperator = Exponentiate(0.5 * std::complex<double>(-DeltaTau, -DeltaT) * BondHamiltonian);

      double Energy = 0.0;

      Lambda2 = Lambda2Inv;
      Info = DoIteration(A, Lambda2, B, EvolOperator, SInfo, BondHamiltonian, Energy);
      Eps = std::max(Epsilon, Info.LargestDiscardedEigenvalue());
      Lambda2Inv = InvertDiagonal(Lambda2, Eps);

      A = prod(A, Lambda2);
      B = delta_shift(prod(Lambda2, B), QShift);
      
      // main iterations
      EvolOperator = Exponentiate(std::complex<double>(-DeltaTau, -DeltaT) * BondHamiltonian);

      //Lambda1 = delta_shift(Lambda1Inv, adjoint(QShift));
      Lambda1 = Lambda1Inv;
      Info = DoIteration(B, Lambda1, A, EvolOperator, SInfo, BondHamiltonian, Energy);
      Eps = std::max(Epsilon, Info.LargestDiscardedEigenvalue());
      Lambda1Inv = InvertDiagonal(Lambda1, Eps);

      B = delta_shift(prod(B, Lambda1), adjoint(QShift));
      A = prod(Lambda1, A);

      std::cout << " BE=" << Energy
                << " Ent=" << Info.TotalEntropy()
                << " NS=" << Info.KeptStates() 
                << " TError=" << Info.TruncationError()
                << " KEigen=" << Info.SmallestKeptEigenvalue()
                << " DeltaT=" << DeltaTau
                << std::endl;

      for (int iter = 1; iter < NumIter; ++iter)
      {
         // do we want to play around with the timestep?
#if 0
	 DeltaTau *= 1.0 - 2e-5;
	 EvolOperator = Exponentiate(std::complex<double>(-DeltaTau, -DeltaT) * BondHamiltonian);
#endif

         Lambda2 = Lambda2Inv;
         //         Lambda2 = delta_shift(Lambda2Inv, adjoint(QShift));
         Info = DoIteration(A, Lambda2, B, EvolOperator, SInfo, BondHamiltonian, Energy);
         Eps = std::max(Epsilon, Info.LargestDiscardedEigenvalue());
         Lambda2Inv = InvertDiagonal(Lambda2, Eps);

	 std::cout << " BE=" << Energy
		   << " Ent=" << Info.TotalEntropy()
		   << " NS=" << Info.KeptStates() 
		   << " TError=" << Info.TruncationError()
		   << " KEigen=" << Info.SmallestKeptEigenvalue()
		   << " DeltaT=" << DeltaTau
		   << std::endl;

         A = prod(A, Lambda2);
         B = delta_shift(prod(Lambda2, B), QShift);
         
         //         Lambda1 = delta_shift(Lambda1Inv, adjoint(QShift));
         Lambda1 = Lambda1Inv;
         Info = DoIteration(B, Lambda1, A, EvolOperator, SInfo, BondHamiltonian, Energy);
         Eps = std::max(Epsilon, Info.LargestDiscardedEigenvalue());
         Lambda1Inv = InvertDiagonal(Lambda1, Eps);


	 std::cout << " BE=" << Energy
		   << " Ent=" << Info.TotalEntropy()
		   << " NS=" << Info.KeptStates() 
		   << " TError=" << Info.TruncationError()
		   << " KEigen=" << Info.SmallestKeptEigenvalue()
		   << " DeltaT=" << DeltaTau
		   << std::endl;

         B = delta_shift(prod(B, Lambda1), adjoint(QShift));
         A = prod(Lambda1, A);

      } // end for loop over main iterations



      // do a half-step, to finish off the 2nd order S-T decomposition
      EvolOperator = Exponentiate(0.5 * std::complex<double>(-DeltaTau, -DeltaT) * BondHamiltonian);

      Lambda2 = Lambda2Inv;
      Info = DoIteration(A, Lambda2, B, EvolOperator, SInfo, BondHamiltonian, Energy);

      pvalue_ptr<InfiniteWavefunction> PsiPtr = new InfiniteWavefunction(MakeWavefunction(A, Lambda2, B, Lambda1,
                                                                                          QShift));
      orthogonalize(*PsiPtr.lock());

      pheap::ShutdownPersistent(PsiPtr);

      ProcControl::Shutdown();
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
