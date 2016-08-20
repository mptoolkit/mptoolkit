// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-dmrg-infinite-klm.cpp
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
#include "matrixproduct/operatoratsite.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/dmrg.h"
#include "mp-algorithms/random_wavefunc.h"
#include "matrixproduct/wavefunc-utils.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"

// the Hamiltonian
#include "models/hubbard-so4.h"
#include "models/spin-su2.h"

namespace prog_opt = boost::program_options;

// The local basis for a unit cell.  This must agree with the choice used in klm-so4.cpp,
// otherwise we cannot use the wavefunction with a klm lattice file.
// Note that the unit cell is 2 real lattice sites, 4 DMRG sites.
std::vector<SiteBasis> ConstructUnitCell()
{
   SiteBlock bSiteA = CreateSO4HubbardSiteA("Q", "S");
   SiteBlock bSiteB = CreateSO4HubbardSiteB("Q", "S");
   SiteBlock bSiteS = CreateSU2SpinSite(0.5, "S");
   // The spin site doesn't contain the Q pseudospin quantum number.
   // Force it so that all sites have the same symmetry list
   // (the pseudospin automatically ends up as zero).
   CoerceSymmetryList(bSiteS, SymmetryList("Q:SU(2),S:SU(2)"));

   std::vector<SiteBasis> Cell;
   Cell.push_back(bSiteS.Basis1());
   Cell.push_back(bSiteA.Basis1());
   Cell.push_back(bSiteS.Basis1());
   Cell.push_back(bSiteB.Basis1());
   return Cell;
}

// The Hamiltonian for L sites.  It is very inefficient to reconstruct the Hamiltonian
// every iteration, but until we have recursive operators it is the only choice.
MPOperator ConstructHamiltonian(int L, double Jk, double Jh)
{
   CHECK(L%2 == 0); // make sure we have an integer number of unit cells

   // This part is copied from klm-so4.cpp
   SiteBlock bSiteA = CreateSO4HubbardSiteA("Q", "S");
   SiteBlock bSiteB = CreateSO4HubbardSiteB("Q", "S");
   SiteBlock bSiteS = CreateSU2SpinSite(0.5, "S");

   Lattice MyLattice(SymmetryList("Q:SU(2),S:SU(2)"));
   for (int i = 1; i <= L; ++i)
   {
      SiteBlock const& bSiteC = (i % 2 == 1) ? bSiteA : bSiteB;  // bipartite

      std::string Site = boost::lexical_cast<std::string>(i);
      MyLattice.Append("0," + Site, bSiteS);
      MyLattice.Append("1," + Site, bSiteC);
   }

   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int, int> S(OpList, "S");
   MPOperator Hamiltonian;

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   // these are the SU(2) coefficients for the dot products
   double HoppingValue = -2.0;
   double SpinValue = -sqrt(3.0);

   for (int i = 1; i < L; ++i)
   {
      // Hopping
      Hamiltonian += HoppingValue * prod(C(1,i), CH(1,i%L+1), Ident);
      // Direct Heisenberg term
      if (Jh != 0) Hamiltonian += (Jh * SpinValue) * prod(S(0,i), S(0,i%L+1), Ident);
   }
   // Kondo coupling, and set the convenience operators
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian += (Jk * SpinValue) * prod(S(0,i), S(1,i), Ident);
   }

   return Hamiltonian;
}

QuantumNumber GetTarget(QuantumNumber const& FinalTarget, int FinalL, int L)
{
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2,QuantumNumbers::SU2>
      QN(FinalTarget.GetSymmetryList());

   half_int s_target = FinalTarget.Get<QuantumNumbers::SU2>("S").j;
   half_int q_target = FinalTarget.Get<QuantumNumbers::SU2>("Q").j;

   TRACE(s_target)(q_target);

   double Desired_s = s_target.to_double() * L / FinalL;
   double Desired_q = q_target.to_double() * L / FinalL;

   // s+q must be an integer
   int SpQ = int(Desired_s + Desired_q + 0.5);
   int Q2 = int(Desired_q * 2 + 0.5);  // twice q
   if (Q2 > SpQ*2) Q2 = SpQ*2;
   int S2 = 2*SpQ - Q2;

   TRACE(Desired_s)(Desired_q)(Q2)(S2);

   return QN(half_int(Q2 / 2.0), half_int(S2 / 2.0));
}

int main(int argc, char** argv)
{
   try
   {
      int NumIter = 4;
      int MinStates = 1;
      int MaxStates = 20;
      double MixFactor = 0.01;
      bool TwoSite = false;
      int NumSweeps = 2;
      double MinTrunc = 0;
      std::string TargetStateStr;
      int TargetL;
      std::string OutputPsi;
      double Jk;
      double Jh = 0;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("max-states,m", prog_opt::value<int>(&MaxStates),
          ("Maximum number of states to keep [default "
           + boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
          ("Minimum number of states to keep [default "
           +boost::lexical_cast<std::string>(MinStates)+"]").c_str())
         ("min-trunc,t", prog_opt::value<double>(&MinTrunc),
          "Minimum desired truncation error (overriden by max-states) [default 0]")
         ("mix-factor,f", prog_opt::value<double>(&MixFactor),
          "Mixing coefficient for the density matrix [default 0.01]")
         ("size,s", prog_opt::value(&TargetL),
          "Target lattice size (must be multiple of unit cell size) [required]")
         ("target,r", prog_opt::value<std::string>(&TargetStateStr),
          "Target quantum number")
         ("out,o", prog_opt::value(&OutputPsi),
          "Filename for the output wavefunction [required]")
         ("Jk", prog_opt::value(&Jk), "Kondo coupling")
         ("Jh", prog_opt::value(&Jh),
          ("Direct Heisenberg coupling [default "
           + boost::lexical_cast<std::string>(Jh) + "]").c_str())
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || !vm.count("size") || !vm.count("out") || !vm.count("Jk"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-dmrg-infinite [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      pheap::Initialize(OutputPsi, 1, PageSize, CacheSize);

      std::cout << "Starting infinite-size DMRG\n";
      std::cout << "Target lattice size is " << TargetL << '\n';

      // Get the unit cell size
      std::vector<SiteBasis> UnitCell = ConstructUnitCell();
      SymmetryList MySList = UnitCell[0].GetSymmetryList();

      QuantumNumber Target(MySList, TargetStateStr);
      std::cout << "Target quantum number is " << Target << '\n';

      // Split the unit cell in half
      int SplitLoc = (UnitCell.size()+1) / 2;
      std::vector<SiteBasis> UnitCellLeft(UnitCell.begin(), UnitCell.begin()+SplitLoc);
      std::vector<SiteBasis> UnitCellRight(UnitCell.begin()+SplitLoc, UnitCell.end());

      std::cout << "Splitting the unit cell into partition ("
                << UnitCellLeft.size() << "," << UnitCellRight.size() << ")\n";

      // Initial lattice
      int L = UnitCell.size();
      std::vector<SiteBasis> LeftStartLattice = UnitCellLeft;
      std::vector<SiteBasis> RightStartLattice = UnitCellRight;
      std::swap(UnitCellLeft, UnitCellRight);
      while (L < 6)
      {
         L += UnitCell.size();
         LeftStartLattice.insert(LeftStartLattice.end(),
                                 UnitCellLeft.begin(), UnitCellLeft.end());
         RightStartLattice.insert(RightStartLattice.begin(),
                                  UnitCellRight.begin(), UnitCellRight.end());
         std::swap(UnitCellLeft, UnitCellRight);
      }
      QuantumNumber CurrentTarget = GetTarget(Target, TargetL, L);

      // initial wavefunction
      MPWavefunction Psi;
      VectorBasis LeftBasis(MySList);
      //LeftBasis.push_back(Target, 1);  with no delta shift
      LeftBasis.push_back(CurrentTarget, 1);
      VectorBasis RightBasis(MySList);
      RightBasis.push_back(QuantumNumber(MySList), 1);
      for (unsigned i = 0; i < LeftStartLattice.size(); ++i)
      {
         Psi.PushLeft(MPStateComponent::ConstructFullBasis2(LeftBasis,
                                                            LeftStartLattice[i].Basis()));
         TRACE(Psi.Left())(scalar_prod(herm(Psi.Left()),Psi.Left()));
         LeftBasis = Psi.Left().Basis2();
      }
      for (int i = RightStartLattice.size()-1; i >= 0; --i)
      {
         Psi.PushRight(MPStateComponent::ConstructFullBasis1(RightStartLattice[i].Basis(),
                                                             RightBasis));
         RightBasis = Psi.Right().Basis1();
      }
      // Initial Center matrix
      Psi.Center() = MakeRandomMatrixOperator(Psi.Left().Basis2(), Psi.Right().Basis1());

      // Make the DMRG object
      DMRG dmrg(Psi, ConstructHamiltonian(L, Jk, Jh));

      // first iteration
      dmrg.Solve(NumIter);
      dmrg.TruncateLeft(MinStates, MaxStates, MinTrunc, MixFactor);
      dmrg.TruncateRight(MinStates, MaxStates, MinTrunc, MixFactor);

      TRACE(Psi.Left())(scalar_prod(herm(Psi.Left()),Psi.Left()));

      TRACE(Psi.LookupLeft(0).Basis1())(Psi.LookupLeft(0).Basis2());
      TRACE(Psi.LookupLeft(0)[0].Basis1())(Psi.LookupLeft(0)[0].Basis2());

      while (L < TargetL)
      {
         // insert sites
         L += UnitCell.size();
         CurrentTarget = GetTarget(Target, TargetL, L);

         std::cout << "Constructing Hamiltonian" << std::endl;
         MPOperator Ham = ConstructHamiltonian(L, Jk, Jh);
         std::cout << "Inserting sites" << std::endl;
         dmrg.InsertSitesDeltaShift(UnitCellLeft, UnitCellRight, CurrentTarget, Ham);

         std::swap(UnitCellLeft, UnitCellRight);

         std::cout << "Solving" << std::endl;
         double Energy = dmrg.Solve(NumIter);
         std::cout << "Truncating" << std::endl;
         dmrg.TruncateLeft(MinStates, MaxStates, MinTrunc, MixFactor);
         dmrg.TruncateRight(MinStates, MaxStates, MinTrunc, MixFactor);

         std::cout << "L=" << L << ", Target=" << CurrentTarget
                   << ", Energy=" << Energy << std::endl;

         DEBUG_TRACE(dmrg.Wavefunction().LeftVacuumBasis());
      }

      std::cout << "Finished." << std::endl;
      pvalue_ptr<MPWavefunction> P = new MPWavefunction(dmrg.Wavefunction());
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
