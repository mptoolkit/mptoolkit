// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-construct.cpp
//
// Copyright (C) 2012-2017 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "common/environment.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "lattice/infinitelattice.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

int
GetNextBasisState(std::string States, int& Pos, SiteBasis LocalBasis)
{
   // Read the local state label.
   int End = States.find_last_not_of(":( ", Pos);
   int Start = States.find_last_of(":) ", End);
   std::string Label = States.substr(Start+1, End-Start);

   int Index;

   // If there is only one local basis state, we do not require it to be specified.
   if (Label == "" && LocalBasis.size() == 1)
      Index = 0;
   else
   {
      Index = LocalBasis.LookupOrNeg(Label);
      // If Label is invalid.
      if (Index < 0)
      {
         std::string ErrorMessage = "fatal: Invalid basis state:\n" + States + '\n';
         for (int i = 0; i <= Start; ++i)
            ErrorMessage += ' ';
         for (int i = 0; i < End-Start; ++i)
            ErrorMessage += '^';
         ErrorMessage += "\nThe valid basis states are:";
         for (int i = 0; i < LocalBasis.size(); ++i)
            ErrorMessage += ' ' + LocalBasis.Label(i);

         throw std::runtime_error(ErrorMessage);
      }
   }

   Pos = Start;

   return Index;
}

QuantumNumber
GetNextBoundaryQN(std::string States, int& Pos, SymmetryList SL, QuantumNumberList QL)
{
   // Read the quantum number label.
   int End = States.find_last_not_of(") ", Pos);
   int Start = States.find_last_of("( ", End);
   std::string Label = States.substr(Start+1, End-Start);

   Pos = States.find_last_not_of("( ", Start) + 1;

   QuantumNumber Q;
   bool Valid = false;
   try
   {
      // If the QN is unspecified and there is only one allowed option, use that.
      if (Label == "" && QL.size() == 1)
         Q = QL[0];
      else
         Q = QuantumNumber(SL, Label);

      // Check whether this is a valid QN by checking against the allowed QNs in QL.
      for (auto const& I : QL)
         if (Q == I)
            Valid = true;

      // If we do not provide a QL (e.g. for the first bond), treat the QN as valid.
      if (QL.size() == 0)
         Valid = true;
   }
   catch (invalid_string_conversion& e)
   {
      Valid = false;
   }

   if (!Valid)
   {
      std::string ErrorMessage = "fatal: Invalid quantum number:\n" + States + '\n';
      for (int i = 0; i <= Start; ++i)
         ErrorMessage += ' ';
      for (int i = 0; i < End-Start; ++i)
         ErrorMessage += '^';
      if (QL.size() > 0)
      {
         ErrorMessage += "\nThe allowed quantum numbers at this bond (given the state and quantum number to the right) are:";
         for (auto const& I : QL)
            ErrorMessage += ' ' + I.ToString();
      }

      throw std::runtime_error(ErrorMessage);
   }

   return Q;
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string LatticeFilename;
      std::string OutputFilename;
      std::string States;
      bool Infinite = false;
      bool Finite = false;
      bool Force = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("lattice,l", prog_opt::value(&LatticeFilename), "Lattice filename [required]")
         ("output,o", prog_opt::value(&OutputFilename), "Output filename [required]")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("infinite", prog_opt::bool_switch(&Infinite), "Create an infinite wavefunction [default]")
         ("finite", prog_opt::bool_switch(&Finite), "Create a finite wavefunction")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("states", prog_opt::value(&States), "states")
         ;

      prog_opt::positional_options_description p;
      p.add("states", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("lattice") == 0 || vm.count("output") == 0 || vm.count("states") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <states>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      if (Infinite && Finite)
      {
         std::cerr << "fatal: Cannot use --infinite and --finite simultaneously!" << std::endl;
         return 1;
      }
      else if (!Infinite && !Finite)
         Infinite = true; // This is the current default behavior.

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFilename);
      SymmetryList SL = Lattice->GetSymmetryList();

      LinearWavefunction Psi;

      if (States.find(':') != -1 && States.find_first_of("()") != -1)
      {
         std::cerr << "fatal: Invalid string: you cannot use both of the \":\" and \"()\" separators." << std::endl;
         return 1;
      }

      // This flag tells us which type of string we are parsing,
      // i.e. one with () separators if true, or : separators if false.
      bool ExplicitQNs = States.find_first_of("()") != -1;

      int Pos = -1;

      QuantumNumber QLast(SL);
      QuantumNumber QPrev = QLast;

      // Get the first quantum number if it is specified.
      if (ExplicitQNs)
      {
         QLast = GetNextBoundaryQN(States, Pos, SL, QuantumNumberList());
         QPrev = QLast;
      }

      VectorBasis BasisPrev(SL);
      BasisPrev.push_back(QPrev, 1);

      auto Site = Lattice->GetUnitCell().end();

      // Iterate backwards from the end.
      do
      {
         --Site;

         // Read the local state label.
         SiteBasis LocalBasis = Site->Basis2();
         int Index = GetNextBasisState(States, Pos, LocalBasis);

         // Get bond quantum number.
         QuantumNumber Q;
         if (ExplicitQNs)
         {
            // The allowed quantum numbers for this bond.
            QuantumNumberList QL = transform_targets(QPrev, LocalBasis.qn(Index));
            Q = GetNextBoundaryQN(States, Pos, SL, QL);
         }
         else
            Q = delta_shift(QPrev, LocalBasis.qn(Index));

         VectorBasis Basis(SL);
         Basis.push_back(Q, 1);

         // Form the local A-matrix.
         StateComponent A(LocalBasis, Basis, BasisPrev);
         A[Index](0,0) = LinearAlgebra::Matrix<double>(1,1,1.0);
         Psi.push_front(A);

         QPrev = Q;
         BasisPrev = Basis;
         if (Site == Lattice->GetUnitCell().begin())
            Site = Lattice->GetUnitCell().end();
      }
      while (Pos > 0);

      if (Site != Lattice->GetUnitCell().end())
      {
         std::cerr << "fatal: The wavefunction unit cell size must be a multiple of the lattice unit cell size ("
                   << Lattice->GetUnitCell().size() << ")." << std::endl;
         return 1;
      }

      MPWavefunction Wavefunction;

      if (Infinite)
      {
         VectorBasis Basis(SL);
         Basis.push_back(QLast, 1);
         RealDiagonalOperator Lambda = RealDiagonalOperator::make_identity(Basis);

         // Find the QShift if the first and last quantum numbers are Abelian.
         // If they are non-Abelian, then we require that the first and last quantum number are the same.
         QuantumNumber QShift(SL);
         QuantumNumberList QL = transform_targets(QPrev, adjoint(QLast));
         if (QL.size() == 1 && QL[0].degree() == 1)
            QShift = QL[0];
         else if (QPrev != QLast)
         {
            std::cerr << "fatal: Cannot have a nontrivial qshift in a non-Abelian quantum number." << std::endl;
            return 1;
         }

         InfiniteWavefunctionLeft PsiOut = InfiniteWavefunctionLeft::ConstructFromOrthogonal(Psi, QShift, Lambda);
         PsiOut.check_structure();

         Wavefunction.Wavefunction() = std::move(PsiOut);
      }
      else if (Finite)
      {
         FiniteWavefunctionLeft PsiOut = FiniteWavefunctionLeft::Construct(Psi);
         PsiOut.check_structure();

         Wavefunction.Wavefunction() = std::move(PsiOut);
      }

      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
      Wavefunction.SetDefaultAttributes();

      pvalue_ptr<MPWavefunction> PsiPtr(new MPWavefunction(Wavefunction));
      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      if (e.Why == "File exists")
         std::cerr << "Note: Use --force (-f) option to overwrite." << std::endl;
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      pheap::Cleanup();
      return 1;
   }
}
