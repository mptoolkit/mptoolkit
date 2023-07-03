// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-aux-matrix.cpp
//
// Copyright (C) 2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mp-algorithms/write-matrix.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/formatting.h"
#include "tensor/regularize.h"
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace prog_opt = boost::program_options;

// returns true if Name exists and is a regular file
bool FileExists(std::string const& Name)
{
   struct stat buf;
   return stat(Name.c_str(), &buf) != -1 && S_ISREG(buf.st_mode);
}

MatrixOperator
Regularize(MatrixOperator const& M)
{
   Regularizer R1(M.Basis1());
   Regularizer R2(M.Basis2());
   return RegularizeBasis12(R1, M, R2);
}

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      std::string OutFile;
      bool Force = false;
      int Verbose = 0;
      bool Quiet = false;
      std::string Format = "python";

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("wavefunction,w", prog_opt::value(&PsiStr), "Wavefunction [required]")
         ("output,o", prog_opt::value(&OutFile), "output file [required]")
         ("quiet,q", prog_opt::bool_switch(&Quiet), "suppress comments in the output")
         ("format,f", prog_opt::value(&Format), FormatDefault("output format", Format).c_str())
         ("force", prog_opt::bool_switch(&Force), "overwrite the output file, if it exists")
         ;

      prog_opt::options_description hidden("Hidden options");

      prog_opt::positional_options_description p;

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1 || vm.count("output") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] -w <psi> -o <file> \n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      Quiet = true;  // TODO: the matlab format doesn't work without this

      // Load the wavefunction
      pvalue_ptr<MPWavefunction> Psi
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      InfiniteWavefunctionLeft InfPsi = Psi->get<InfiniteWavefunctionLeft>();

      // orthogonalize the wavefunction
      LinearWavefunction PsiL;
      RealDiagonalOperator D;
      std::tie(PsiL, D) = get_left_canonical(InfPsi);
      MatrixOperator Rho = D;
      Rho = Regularize(scalar_prod(Rho, herm(Rho)));
      if (Rho.size1() != 1 || Rho.size2() != 1)
      {
         std::cerr << "mp-matrix: error: mps has non-trivial symmetries.\n";
         TRACE(Rho.Basis1());
         return 1;
      }

      if (!Force && FileExists(OutFile))
      {
         std::cerr << "mp-matrix: error: output file " << OutFile << " already exists, "
            "use --force to overwrite.\n";
         return 1;
      }

      std::ofstream Out(OutFile, std::ios::out | std::ios::trunc);
      if (!Out.good())
      {
         std::cerr << "mp-matrix: failed to open file " << OutFile << "\n";
         return 1;
      }

      Out.precision(getenv_or_default("MP_PRECISION", 14));

      if (Format == "matlab")
      {
         WriteMPS_MATLAB(Out, PsiL, Rho, Quiet);
      }
      else if (Format == "python")
      {
         WriteMPS_Python(Out, PsiL, Rho, Quiet);
      }
      else
      {
         std::cerr << "mp-matrix: error: unknown output format.\n";
         return 1;
      }

      pheap::Shutdown();
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
