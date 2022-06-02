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

#include "wavefunction/mpwavefunction.h"
#include "mps/packunpack.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include "linearalgebra/arpack_wrapper.h"
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>
#include "common/formatting.h"
#include <tuple>

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
   MatrixOperator U = Regularize(M.Basis1());
   MatrixOperator V = Regularize(M.Basis2());
   return U * M * herm(V);
}

StateComponent
Regularize(StateComponent const& M)
{
   MatrixOperator U = Regularize(M.Basis1());
   MatrixOperator V = Regularize(M.Basis2());
   return prod(U, prod(M, herm(V)));
}

//
// MATLAB output format
//

void WriteMatrixFormatMATLAB(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                             bool Quiet = false)
{
   out << "complex( [\n";
   if (!Quiet)
      "% real part\n ";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << ";\n    ";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << " ";
         out << M(i,j).real();
      }
   }
   out << "  ],";
   out << "[\n";
   if (!Quiet)
      out << "\n%imaginary part\n  ";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << ";\n    ";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << " ";
         out << M(i,j).imag();
      }
   }
   out << "  ] )";
}

void WriteRealMatrixFormatMATLAB(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                                 bool quiet = false)
{
   out << "[ ";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << ";\n    ";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << " ";
         out << M(i,j).real();
      }
   }
   out << "  ]";
}

void WriteMPS_MATLAB(std::ostream& out, LinearWavefunction const& Psi, MatrixOperator const& Rho, bool Quiet = false)
{
   LinearWavefunction::const_iterator I = Psi.begin();
   int Site = 0;
   out << "MPS = cell(1," << Psi.size() << ");\n";
   while (I != Psi.end())
   {
      if (!Quiet)
         out << "% site " << Site << "\n";
      //         out << "[\n";
      out << "MPS{1," << (Site+1) << "} = cat(3, ... \n";
      StateComponent A = Regularize(*I);
      for (int i = 0; i < A.size(); ++i)
      {
         if (!Quiet)
            out << "% site " << Site << "  basis state " << i << "\n";
         if (i != 0)
            out << ", ...\n";
         if (A.Basis1().size() != 1 || A.Basis2().size() != 1)
         {
            throw std::runtime_error("mp-matrix: error: mps has non-trivial symmetries");
         }
         WriteMatrixFormatMATLAB(out, A[i](0,0), Quiet);
         //out << "\n";
      }
      if (!Quiet)
         out << "% end of site " << Site << "\n";
      out << ");\n";

      ++I;
      ++Site;
   }
   out << "\n";
   if (!Quiet)
      out << "% end of MPS\n";

   if (!Quiet)
      out << "\n% density matrix\n";
   out << "RHO = ";
   WriteRealMatrixFormatMATLAB(out, Rho(0,0), Quiet);
   out << ";\n";
   if (!Quiet)
      out << "% end of density matrix\n";
}

//
// Python output format
//

void WriteMatrixFormatPython(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                             std::string Prefix = "", bool Quiet = false)
{
   out << "np.array([[";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << "],\n" << Prefix << " [";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << ", ";
         out << M(i,j).real() << '+' << M(i,j).imag() << 'j';
      }
   }
   out << "]\n" << Prefix << "])";
}

void WriteRealMatrixFormatPython(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                                 std::string Prefix = "", bool quiet = false)
{
   out << "np.array([[";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << "],\n" << Prefix << " [";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << ", ";
         out << M(i,j).real();
      }
   }
   out << "]\n" << Prefix << "])";
}

void WriteMPS_Python(std::ostream& out, LinearWavefunction const& Psi, MatrixOperator const& Rho, bool Quiet = false)
{
   LinearWavefunction::const_iterator I = Psi.begin();
   int Site = 0;
   out <<"import numpy as np\n";
   out << "MPS = [\n";
   while (I != Psi.end())
   {
      if (!Quiet)
         out << "# site " << Site << "\n";
      out << " np.array([";
      StateComponent A = Regularize(*I);
      for (int i = 0; i < A.size(); ++i)
      {
         if (!Quiet)
            out << "  # site " << Site << "  basis state " << i << "\n";
         if (i != 0)
            out << " ,\n";
         if (A.Basis1().size() != 1 || A.Basis2().size() != 1)
         {
            throw std::runtime_error("mp-matrix: error: mps has non-trivial symmetries");
         }
         out << "  ";
         WriteMatrixFormatPython(out, A[i](0,0), "  ", Quiet);
         out << "\n";
      }
      if (!Quiet)
         out << "# end of site " << Site << "\n";
      out << "])\n";

      ++I;
      ++Site;
      if (I != Psi.end())
      {
         out << ",\n";
      }
   }
   out << "]\n";
   if (!Quiet)
      out << "# end of MPS\n";

   if (!Quiet)
      out << "\n# density matrix\n";
   out << "RHO = ";
   WriteRealMatrixFormatPython(out, Rho(0,0), "", Quiet);
   out << "\n";
   if (!Quiet)
      out << "# end of density matrix\n";
}

//
//
//

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
      MatrixOperator Rho = Regularize(D);
      Rho = scalar_prod(Rho, herm(Rho));
      if (Rho.size1() != 1 || Rho.size2() != 1)
      {
         std::cerr << "mp-matrix: error: mps has non-trivial symmetries.\n";
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
