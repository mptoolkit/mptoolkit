// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-idmrg-old.cpp
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
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"

namespace prog_opt = boost::program_options;

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

// function to determine the target states.  This might return
// multiple target states, in which case the map also gives the
// density matrix weights (NOT absolute normalized).
// The algorithm requires the target weight
std::map<QuantumNumbers::QuantumNumber, double>
FindTargetStates(VectorBasis const& B1, VectorBasis const& B2, double TargetWeight)
{
   using QuantumNumbers::QuantumNumber;
   using QuantumNumbers::Projection;

   std::map<QuantumNumber, double> Result;
   if (B1.is_empty() || B2.is_empty())  // our algorithm requires the basis is non-empty
      return Result;

   std::vector<Projection> B1Proj(B1.size()), B2Proj(B2.size());

   double ClosestWeight = 1E100;
   double const CloseEpsilon = 1E-10;  // this needs to be smaller than 1 / (largest possible lattice size)

   for (unsigned i = 0; i < B1.size(); ++i)
   {
      for (unsigned j = 0; j < B2.size(); ++j)
      {
         Projection p = QuantumNumbers::difference(B1[i], B2[j]);
         double ThisDiff = std::abs(TargetWeight - weight(p));
         if (ThisDiff < ClosestWeight - CloseEpsilon)
         {
            ClosestWeight = ThisDiff;
            Result.clear();
            Result[heighest_weight(p)] = 1.0;
            TRACE(p)(ThisDiff);
         }
         else if (ThisDiff < ClosestWeight + CloseEpsilon)
         {
            Result[heighest_weight(p)] = 1.0;
         }
      }
   }

   TRACE(TargetWeight)(ClosestWeight);

   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      int NumIter = 20;
      int MinStates = 50;
      int MaxStates = 100000;
      double MixFactor = 0.01;
      int NumRightSites = 2;
      double MinTrunc = 0;
      std::string TargetQStr;
      std::string HamiltonianStr;
      std::string OutPsiStr;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamiltonianStr),
          "Hamiltonian operator")
         ("out,o", prog_opt::value(&OutPsiStr),
          "output wavefunction filename")
         ("right-sites", prog_opt::value(&NumRightSites),
          "number of sites in the right block")
         ("iter,i", prog_opt::value<int>(&NumIter), ("Number of Lanczos iterations per step [default "
          +boost::lexical_cast<std::string>(NumIter)+"]").c_str())
         ("target,q", prog_opt::value(&TargetQStr), "target quantum number")
         ("max-states,m", prog_opt::value<int>(&MaxStates), "Maximum number of states to keep [default 100000]")
         ("min-states", prog_opt::value<int>(&MinStates), "Minimum number of states to keep [default 50]")
         ("min-trunc,t", prog_opt::value<double>(&MinTrunc),
          "Minimum desired truncation error (overriden by max-states) [default 0]")
         ("mix-factor,f", prog_opt::value<double>(&MixFactor), "Mixing coefficient for the density matrix [default 0.01]")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("out") == 0 || vm.count("Hamiltonian") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-idmrg [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Required options: out, target, Hamiltonian\n";
         return 1;
      }

      std::cout << "Creating wavefunction from infinite DMRG...\n";

      int PageSize = getenv_or_default("MP_PAGESIZE", DEFAULT_PAGE_SIZE);
      long CacheSize = getenv_or_default("MP_CACHESIZE", DEFAULT_PAGE_CACHE_SIZE);
      pheap::Initialize(OutPsiStr, 1, PageSize, CacheSize);

      OperatorList OpList;
      MPOperator H;
      std::tie(OpList, H) = ParseLatticeAndOperator(HamiltonianStr);
      Lattice const& Lat = OpList.GetLattice();

      // target quantum number
      QuantumNumbers::QuantumNumber TargetState(Lat.GetSymmetryList(), TargetQStr);

      // weight of the target state quantum number
      double TargetWeight = weight(difference(TargetState, QuantumNumber(Lat.GetSymmetryList())));
      TRACE(TargetWeight);

      // left-most basis
      BasisList Vacuum = make_vacuum_basis(Lat.GetSymmetryList());
      VectorBasis LeftVacuumBasis(Lat.GetSymmetryList());
      LeftVacuumBasis.push_back(TargetState, 1);
      TRACE(LeftVacuumBasis);

      TRACE(Lat.GetSymmetryList());
      CenterWavefunction Psi;

      Lattice::const_iterator LeftLatIter = Lat.begin();
      LinearOperator::const_iterator HamIter = H.begin();

      // Hamiltonian E matrices
      MPStateComponent Elem = MPStateComponent(Vacuum, LeftVacuumBasis, LeftVacuumBasis);

      Elem[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1);
      SuperblockOperator HamMatrices;
      HamMatrices.PushLeft(Elem); // add the vacuum operator, for convenience

      // initial block
      Psi.PushLeft(ConstructFromLeftBasis(LeftLatIter->Basis1().Basis(), LeftVacuumBasis));
      Elem = operator_prod(herm(*HamIter),
                           herm(Psi.Left()),
                           Elem,
                           Psi.Left());
      HamMatrices.PushLeft(Elem);

      ++LeftLatIter;  // 2 sites in initial left block
      ++HamIter;
      Psi.PushLeft(ConstructFromLeftBasis(LeftLatIter->Basis1().Basis(), Psi.Left().Basis2()));
      Elem = operator_prod(herm(*HamIter),
                           herm(Psi.Left()),
                           Elem,
                           Psi.Left());
      HamMatrices.PushLeft(Elem);

      // construct the right block
      Lattice::const_iterator RightLatIter = LeftLatIter;
      LinearOperator::const_iterator RightHamIter = HamIter;
      std::advance(RightLatIter, NumRightSites);
      std::advance(RightHamIter, NumRightSites);

      // what quantum number to choose for the right vacuum basis?
      VectorBasis RightVacuumBasis(Lat.GetSymmetryList());
      RightVacuumBasis.push_back(QuantumNumbers::QuantumNumber(Lat.GetSymmetryList()), 1);

      LinearOperator::const_iterator RHI = H.end();
      Lattice::const_iterator RLI = Lat.end();

      MPStateComponent RElem(Vacuum, RightVacuumBasis, RightVacuumBasis);
      RElem[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1);
      HamMatrices.PushRight(RElem);
      --RHI;
      --RLI;
      while (RHI != RightHamIter)
      {
         RElem = local_prod((1.0 / RHI->SiteBasis().total_degree()) * local_trace(*RHI),
                            RElem);
         HamMatrices.PushRight(RElem);
         --RHI;
         --RLI;
      }

      Psi.PushRight(ConstructFromRightBasis(RightLatIter->Basis1().Basis(), RightVacuumBasis));
      RElem = operator_prod(*RightHamIter,
                            Psi.Right(),
                            RElem,
                            herm(Psi.Right()));
      HamMatrices.PushRight(RElem);

      --RightLatIter;
      --RightHamIter;

      // Keep building up the initial right block
      while (RightLatIter != LeftLatIter)
      {
         Psi.PushRight(ConstructFromRightBasis(RightLatIter->Basis1().Basis(), Psi.Right().Basis1()));
         RElem = operator_prod(*RightHamIter,
                               Psi.Right(),
                               RElem,
                               herm(Psi.Right()));
         HamMatrices.PushRight(RElem);

         --RightLatIter;
         --RightHamIter;
      }

      TRACE(Psi.Left().Basis2())(Psi.Right().Basis1());

      double ThisWeight = (1.0 - 4.0 / Lat.size()) * TargetWeight;
      TRACE(ThisWeight);
      std::map<QuantumNumbers::QuantumNumber, double> Targets = FindTargetStates(Psi.Left().Basis2(),
                                                                                 Psi.Right().Basis1(),
                                                                                 ThisWeight);
      for (std::map<QuantumNumbers::QuantumNumber, double>::const_iterator I = Targets.begin(); I != Targets.end(); ++I)
      {
         TRACE(I->first);
      }

      Psi.Center() = MakeRandomMatrixOperator(Psi.Left().Basis2(), Psi.Right().Basis1(), Targets.begin()->first);

      int Iterations = NumIter;
      double Energy = Lanczos(Psi.Center(),
                              SuperblockMultiply(HamMatrices.Left(),
                                                 HamMatrices.Right()),
                              Iterations);

      TRACE(Energy);

      // do the next iteration


      for (int n = 0; n < 4; ++n)
      {
         // extend the left block by one site
         ++LeftLatIter;
         ++HamIter;
         Psi.PushLeft(ConstructFromLeftBasis(LeftLatIter->Basis1().Basis(), Psi.Left().Basis2()));
         Elem = operator_prod(herm(*HamIter),
                              herm(Psi.Left()),
                              Elem,
                              Psi.Left());
         HamMatrices.PushLeft(Elem);
         HamMatrices.PopRight();
         Psi.PopRight();

         // new right blocks
         while (Psi.RightSize() > 0)
         {
            Psi.PopRight();
            HamMatrices.PopRight();
         }

         RightLatIter = LeftLatIter;
         RightHamIter = HamIter;
         std::advance(RightLatIter, NumRightSites);
         std::advance(RightHamIter, NumRightSites);

         Psi.PushRight(ConstructFromRightBasis(RightLatIter->Basis1().Basis(), RightVacuumBasis));
         RElem = HamMatrices.Right();
         RElem = operator_prod(*RightHamIter,
                               Psi.Right(),
                               RElem,
                               herm(Psi.Right()));
         HamMatrices.PushRight(RElem);
         --RightLatIter;
         --RightHamIter;

         RElem = HamMatrices.Right();
         while (RightLatIter != LeftLatIter)
         {
            Psi.PushRight(ConstructFromRightBasis(RightLatIter->Basis1().Basis(), Psi.Right().Basis1()));
            RElem = operator_prod(*RightHamIter,
                                  Psi.Right(),
                                  RElem,
                                  herm(Psi.Right()));
            HamMatrices.PushRight(RElem);
            --RightLatIter;
            --RightHamIter;
         }

         RightLatIter = LeftLatIter;

      ThisWeight = (1.0 - 5.0 / Lat.size()) * TargetWeight;
      TRACE(ThisWeight);
      Targets = FindTargetStates(Psi.Left().Basis2(), Psi.Right().Basis1(), ThisWeight);
      for (std::map<QuantumNumbers::QuantumNumber, double>::const_iterator I = Targets.begin(); I != Targets.end(); ++I)
      {
         TRACE(I->first);
      }


      Psi.Center() = MakeRandomMatrixOperator(Psi.Left().Basis2(), Psi.Right().Basis1(), Targets.begin()->first);
      Iterations = NumIter;
      Energy = Lanczos(Psi.Center(),
                       SuperblockMultiply(HamMatrices.Left(),
                                          HamMatrices.Right()),
                       Iterations);
      TRACE(Energy);

      }

      pvalue_ptr<MPWavefunction> P = pvalue_ptr<MPWavefunction>(new MPWavefunction(Psi.AsLinearWavefunction()));
      pheap::ShutdownPersistent(P);
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
