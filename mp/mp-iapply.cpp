// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-iapply.cpp
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcell-parser.h"
#include "lattice/infinite-parser.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OpStr;
      std::string InputFile;
      std::string OutputFile;
      bool Force = false;
      bool AssumeOrthogonal = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("assume-orthogonal", prog_opt::bool_switch(&AssumeOrthogonal),
          "skip the orthogonalization step")
         ("force,f", prog_opt::bool_switch(&Force),
          "allow overwriting the output file, if it already exists")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("op", prog_opt::value(&OpStr), "op")
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("op", 1);
      p.add("psi1", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") < 1)
      {
         print_copyright(std::cerr, "tools", "mp-iapply");
         std::cerr << "usage: " << basename(argv[0]) << " [options] <operator> <input-psi> <output-psi>\n";
         std::cerr << desc << '\n';
         std::cerr << "This tool calculates the action of an operator on an iMPS.\n";
         std::cerr << "The operator must be of the ProductMPO form.\n";
         std::cerr << "If the operator is unitary and the boundary of the unit cell is invariant,\n"
                   << "that is, the operator is a ProductMPO with 1-dimensional boundaries\n"
                   << "(use mp-show-operator to check), then the final wavefunction will already be\n"
                   << "orthogonal already, and the orthonormalization step can be ommitted by using\n"
                   << "the --assume-orthogonal option.\n";

         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFile == OutputFile)
         PsiPtr = pheap::OpenPersistent(InputFile.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFile);
      }

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

      InfiniteLattice Lattice;
      ProductMPO StringOp;
      std::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(OpStr);

      if (Psi.size() != StringOp.size())
      {
         int Size = statistics::lcm(Psi.size(), StringOp.size());
         if (Psi.size() < Size)
         {
            std::cerr << "mp-iapply: warning: extending wavefunction size to lowest common multiple, which is "
                      << Size << " sites.\n";
         }
         Psi = repeat(Psi, Size / Psi.size());
         StringOp = repeat(StringOp, Size / StringOp.size());
      }

      LinearWavefunction PsiL = get_left_canonical(Psi).first;

      if (Verbose > 0)
         std::cout << "Applying operator..." << std::endl;

      // apply the string operator
      ProductMPO::const_iterator MI = StringOp.begin();
      for (LinearWavefunction::iterator I = PsiL.begin(); I != PsiL.end(); ++I, ++MI)
      {
         (*I) = aux_tensor_prod(*MI, *I);
      }

      if (AssumeOrthogonal)
      {
         RealDiagonalOperator Lambda = Psi.lambda_r();
         PsiPtr.mutate()->Wavefunction() = InfiniteWavefunctionLeft::ConstructFromOrthogonal(PsiL, Psi.qshift(), Lambda, Psi.log_amplitude(), Verbose);
      }
      else
      {
         RealDiagonalOperator Lambda = Psi.lambda_r();
         MatrixOperator Rho = Lambda*Lambda;
         PsiPtr.mutate()->Wavefunction() = InfiniteWavefunctionLeft::Construct(PsiL, Psi.qshift(), Rho, Psi.log_amplitude(), Verbose);
      }

      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      if (Verbose > 0)
         std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
