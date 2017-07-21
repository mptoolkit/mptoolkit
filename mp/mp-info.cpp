// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-info.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int MaxEigenvalues = 1000000;
bool ShowBasic = false;
bool Base2 = false;
bool ShowEntropy = false;
bool ShowStates = false;
bool ShowLocalBasis = false;
bool ShowBasis = false;
bool ShowCasimir = false;
bool ShowDensity = false;
bool ShowDegen = false;
bool Quiet = false;
std::vector<int> Partition;

void ShowBasicInfo(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   out << "Wavefunction is an InfiniteWavefunction in the left canonical basis.\n";
   out << "Symmetry list = " << Psi.GetSymmetryList() << '\n';
   out << "Unit cell size = " << Psi.size() << '\n';
   out << "Quantum number per unit cell = " << Psi.qshift() << '\n';
   out << "Number of states = " << Psi.Basis1().total_dimension() << '\n';
   out << std::endl;
}

void ShowBasicInfo(IBCWavefunction const& Psi, std::ostream& out)
{
   out << "Wavefunction is an Infinite Boundary Condition wavefunction in the left/left/right canonical basis.\n";
   out << "Symmetry list = " << Psi.Window.GetSymmetryList() << '\n';
   out << "Left semi-infinite strip unit cell size = " << Psi.Left.size() << '\n';
   if (!Psi.Left.empty())
   {
      out << "Quantum number per unit cell (left) = " << Psi.Left.qshift() << '\n';
   }
   out << "Right semi-infinite strip unit cell size = " << Psi.Right.size() << '\n';
   if (!Psi.Right.empty())
   {
      out << "Quantum number per unit cell (right) = " << Psi.Right.qshift() << '\n';
   }

   out << "Window size = " << Psi.Window.size() << '\n';
   out << "Offset of first site of the window = " << Psi.window_offset() << '\n';

   out << "Number of states (left edge of window) = " << Psi.Window.Basis1().total_dimension() << '\n';
   out << std::endl;
}

void ShowBasicInfo(FiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   out << "Wavefunction is a FiniteWavefunction in the left canonical basis.\n";
   out << "Symmetry list = " << Psi.GetSymmetryList() << '\n';
   out << "Length = " << Psi.size() << '\n';
   out << "Quantum number = " << Psi.TransformsAs() << '\n';
   out << "Norm = " << norm_2(Psi) << '\n';
   out << std::endl;
}

void ShowStateInfo(CanonicalWavefunctionBase const& Psi, std::ostream& out)
{
   out << "#partition  #dimension  #degree\n";

   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      VectorBasis B = Psi.lambda(*I).Basis1();
      out << std::setw(10) << (*I) << "  "
          << std::setw(10) << B.total_dimension() << "  "
          << std::setw(7) <<  B.total_degree()
          << '\n';
   }
   out << std::endl;
}

void ShowBasisInfo(CanonicalWavefunctionBase const& Psi, std::ostream& out)
{
   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      VectorBasis B = Psi.lambda(*I).Basis1();
      out << "Basis at partition " << (*I) << ":\n" << B << '\n';
   }
   out << std::endl;
}

void ShowEntropyInfo(CanonicalWavefunctionBase const& Psi, std::ostream& out)
{
   if (!Quiet)
      out << "#left-size #right-size   #entropy" << (Base2 ? "(base-2)" : "(base-e)") << '\n';
   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      RealDiagonalOperator Rho = Psi.lambda(*I);
      Rho = scalar_prod(Rho, herm(Rho));
      DensityMatrix<MatrixOperator> DM(Rho);

      double Entropy = DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2);
      out << std::setw(10) << (*I) << ' ' << std::setw(11) << (Psi.size() - (*I)) << ' '
          << std::setw(18) << Entropy << '\n';
   }
   out << std::endl;
}

void ShowDM(CanonicalWavefunctionBase const& Psi, std::ostream& out)
{
   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      MatrixOperator Rho = Psi.lambda(*I);
      Rho = scalar_prod(Rho, herm(Rho));

      DensityMatrix<MatrixOperator> DM(Rho);

      if (!Quiet)
	 out << "#Reduced density matrix at partition ("
	     << (*I) << "," << (Psi.size() - *I) << ") :\n";
      DM.DensityMatrixReport(out, MaxEigenvalues, Base2, ShowDegen, Quiet);
      out << std::endl;
   }
}

void ShowCasimirInfo(CanonicalWavefunctionBase const& Psi, std::ostream& out)
{
   out << "#left-size #right-size   ";
   QuantumNumbers::SymmetryList SList = Psi.GetSymmetryList();
   int NumCasimir = SList.NumCasimirOperators();
   for (int i = 0; i < NumCasimir; ++i)
   {
      if (i != 0)
         out << " ";
      std::string Name = "#" + SList.CasimirName(i);
      out << std::setw(20) << Name << ' '
          << std::setw(20) << std::string(Name+"^2");
   }
   out << '\n';

   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      MatrixOperator Rho = Psi.lambda(*I);
      Rho = scalar_prod(Rho, herm(Rho));

      DensityMatrix<MatrixOperator> DM(Rho);

      out << std::setw(10) << (*I) << ' ' << std::setw(11) << (Psi.size() - *I) << "   ";
      for (int i = 0; i < NumCasimir; ++i)
      {
         if (i != 0)
            out << ' ';
         out << std::setw(20) << DM.EvaluateCasimir(i);
         out << std::setw(20) << DM.EvaluateCasimirMoment(i);
      }
      out << '\n';
   }
   out << std::endl;
}

void ShowLocalBasisInfo(CanonicalWavefunctionBase const& Psi, std::ostream& out)
{
   for (int i = 0; i < Psi.size(); ++i)
   {
      out << "Basis at site " << i << '\n';
      out << Psi[i].LocalBasis() << std::endl;
   }
   out << std::endl;
}

struct ShowWavefunctionBasicInfo : public boost::static_visitor<void>
{
   template <typename T>
   void operator()(T const& Psi) const
   {
      ShowBasicInfo(Psi, std::cout);
   }
};


struct ShowWavefunction : public boost::static_visitor<void>
{
   void operator()(InfiniteWavefunctionLeft const& Psi) const;
   void operator()(IBCWavefunction const& Psi) const;
   void operator()(FiniteWavefunctionLeft const& Psi) const;
};

void
ShowWavefunction::operator()(InfiniteWavefunctionLeft const& Psi) const
{
   std::sort(Partition.begin(), Partition.end());
   if (Partition.empty())
   {
      // all partitions
      for (int i = 0; i <= Psi.size(); ++i)
         Partition.push_back(i);
   }

   if (ShowStates)
      ShowStateInfo(Psi, std::cout);

   if (ShowBasis)
      ShowBasisInfo(Psi, std::cout);

   if (ShowEntropy)
      ShowEntropyInfo(Psi, std::cout);

   if (ShowDensity)
      ShowDM(Psi, std::cout);

   if (ShowCasimir)
      ShowCasimirInfo(Psi, std::cout);

   if (ShowLocalBasis)
      ShowLocalBasisInfo(Psi, std::cout);
}

void
ShowWavefunction::operator()(IBCWavefunction const& Psi) const
{
   std::sort(Partition.begin(), Partition.end());
   if (Partition.empty())
   {
      // all partitions in the window
      for (int i = 0; i <= Psi.window_size(); ++i)
         Partition.push_back(i);
   }

   if (ShowStates)
      ShowStateInfo(Psi.Window, std::cout);

   if (ShowBasis)
      ShowBasisInfo(Psi.Window, std::cout);

   if (ShowEntropy)
      ShowEntropyInfo(Psi.Window, std::cout);

   if (ShowDensity)
      ShowDM(Psi.Window, std::cout);

   if (ShowCasimir)
      ShowCasimirInfo(Psi.Window, std::cout);

   if (ShowLocalBasis)
      ShowLocalBasisInfo(Psi.Window, std::cout);
}

void
ShowWavefunction::operator()(FiniteWavefunctionLeft const& Psi) const
{
   std::sort(Partition.begin(), Partition.end());
   if (Partition.empty())
   {
      // all partitions
      for (int i = 0; i <= Psi.size(); ++i)
         Partition.push_back(i);
   }

   if (ShowStates)
      ShowStateInfo(Psi, std::cout);

   if (ShowBasis)
      ShowBasisInfo(Psi, std::cout);

   if (ShowEntropy)
      ShowEntropyInfo(Psi, std::cout);

   if (ShowDensity)
      ShowDM(Psi, std::cout);

   if (ShowCasimir)
      ShowCasimirInfo(Psi, std::cout);

   if (ShowLocalBasis)
      ShowLocalBasisInfo(Psi, std::cout);
}

int main(int argc, char** argv)
{
   try
   {
      std::string WavefuncFile;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("entropy,e", prog_opt::bool_switch(&ShowEntropy), "show the entropy at each partition")
         ("states,s", prog_opt::bool_switch(&ShowStates), "show the number of states at each partition")
         ("basis,a", prog_opt::bool_switch(&ShowBasis), "show the complete basis at each partition")
         ("density-matrix,d", prog_opt::bool_switch(&ShowDensity), "show the density matrix eigenvalues")
         ("degen", prog_opt::bool_switch(&ShowDegen),
          "Show degeneracies in the density matrix as repeated eigenvalues (implies -d)")
         ("limit,l", prog_opt::value<int>(&MaxEigenvalues),
          "limit the density matrix display to N eigenvalues (implies -d)")
         ("casimir,c", prog_opt::bool_switch(&ShowCasimir),
          "show the values of the casimir invariant operators and 2nd moments at each partition")
         ("localbasis,b", prog_opt::bool_switch(&ShowLocalBasis),
          "Show the local basis at each site")
         ("base2,2", prog_opt::bool_switch(&Base2), "show the entropy using base 2 instead of base e")
         ("partition,p", prog_opt::value(&Partition),
          "show quantities only for this parition (zero-based, can be used more than once; use --partition 0 to show only quantities at the edge of the unit cell")
	 ("quiet", prog_opt::bool_switch(&Quiet), "don't show column headings")
         ("warranty", "show the complete explanation as to why there is NO WARRANTY, to the extent permitted by law")
         ("copying", "This program is free software, show the conditions under which it may be copied or modified")
         ("citations", "show information about the citations that publications using this software should reference")
         ;
      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("input-wavefunction", prog_opt::value(&WavefuncFile), "input wavefunction (required)")
         ;

      prog_opt::positional_options_description p;
      p.add("input-wavefunction", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("warranty"))
         print_warranty(std::cout);

      if (vm.count("copying"))
         print_copying(std::cout);

      if (vm.count("citations"))
         print_citations(std::cout);

      if (vm.count("warranty") || vm.count("copying") || vm.count("citations"))
         return 0;

      if (vm.count("help") || vm.count("input-wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <wavefunction>\n";
         std::cerr << desc << "\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (vm.count("limit"))
         ShowDensity = true;

      if (ShowDegen)
         ShowDensity = true;

      if (!ShowEntropy && !ShowStates && !ShowLocalBasis && !ShowBasis && !ShowCasimir && !ShowDensity)
         ShowBasic = true;

      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(WavefuncFile, mp_pheap::CacheSize(), true);

      if (ShowBasic)
      {
         boost::apply_visitor(ShowWavefunctionBasicInfo(), Psi->Wavefunction());
	 std::cout << "File format version " << Psi->version() << '\n';

         std::cout << "Attributes:\n" << Psi->Attributes();

         std::cout << "\nLast history entry:\n" << Psi->History().back() << "\n\n";
      }

      boost::apply_visitor(ShowWavefunction(), Psi->Wavefunction());

   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
   }
   catch(std::exception& e)
   {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
   }
   catch(...)
   {
      std::cerr << "Exception of unknown type!\n";
   }
   return 0;
}
