// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-info.cpp
//
// Copyright (C) 2004-2022 Ian McCulloch <ian@qusim.net>
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
   out << "Infinite wavefunction in the left canonical basis.\n";
   out << "Symmetry list = " << Psi.GetSymmetryList() << '\n';
   out << "Unit cell size = " << Psi.size() << '\n';
   out << "Quantum number per unit cell = " << Psi.qshift() << '\n';
   out << "Log amplitude per unit cell = " << formatting::format_complex(Psi.log_amplitude()) << '\n';
   out << "Number of states = " << Psi.Basis1().total_dimension() << '\n';
}

void ShowBasicInfo(InfiniteWavefunctionRight const& Psi, std::ostream& out)
{
   out << "Infinite wavefunction in the right canonical basis.\n";
   out << "Symmetry list = " << Psi.GetSymmetryList() << '\n';
   out << "Unit cell size = " << Psi.size() << '\n';
   out << "Quantum number per unit cell = " << Psi.qshift() << '\n';
   out << "Log amplitude per unit cell = " << formatting::format_complex(Psi.log_amplitude()) << '\n';
   out << "Number of states = " << Psi.Basis1().total_dimension() << '\n';
}

void ShowBasicInfo(IBCWavefunction const& Psi, std::ostream& out)
{
   out << "Infinite Boundary Condition wavefunction in the left/left/right orthogonal basis.\n";
   out << "Symmetry list = " << Psi.window().GetSymmetryList() << '\n';

   if (Psi.left().empty())
   {
      out << "No left strip.\n";
   }
   else
   {
      std::string lw = Psi.get_left_filename();
      if (lw.empty())
         out << "Left semi-infinite strip is stored directly.\n";
      else
         out << "Left semi-infinite strip is stored in the file \"" << lw << "\"\n";
      out << "Left semi-infinite strip unit cell size = " << Psi.left().size() << '\n';

      out << "Sites incorporated into the window from the left unit cell = " << Psi.window_left_sites() << '\n';
      out << "Quantum number per unit cell (left) = " << Psi.left().qshift() << '\n';
      out << "Left quantum number shift = " << Psi.left_qshift() << '\n';
   }

   if (Psi.right().empty())
   {
      out << "No right strip.\n";
   }
   else
   {
      std::string rw = Psi.get_right_filename();
      if (rw.empty())
         out << "Right semi-infinite strip is stored directly.\n";
      else
         out << "Right semi-infinite strip is stored in the file \"" << rw << "\"\n";
      out << "Right semi-infinite strip unit cell size = " << Psi.right().size() << '\n';

      out << "Sites incorporated into the window from the right unit cell = " << Psi.window_right_sites() << '\n';
      out << "Quantum number per unit cell (right) = " << Psi.right().qshift() << '\n';
      out << "Right quantum number shift = " << Psi.right_qshift() << '\n';
   }

   out << "Window size = " << Psi.window_size() << '\n';
   out << "Offset of first site of the window = " << Psi.window_offset() << '\n';

   out << "Number of states (left edge of window) = " << Psi.window().Basis1().total_dimension() << '\n';
}

void ShowBasicInfo(EAWavefunction const& Psi, std::ostream& out)
{
   out << "Excitation Ansatz wavefunction in the left/left/right orthogonal basis.\n";
   out << "Symmetry list = " << Psi.left().GetSymmetryList() << '\n';

   if (Psi.left().empty())
   {
      out << "No left strip.\n";
   }
   else
   {
      std::string lw = Psi.get_left_filename();
      if (lw.empty())
         out << "Left semi-infinite strip is stored directly.\n";
      else
         out << "Left semi-infinite strip is stored in the file \"" << lw << "\"\n";
      out << "Left semi-infinite strip unit cell size = " << Psi.left().size() << '\n';
      out << "Left unit cell starting index = " << Psi.left_index() << '\n';
      out << "Quantum number per unit cell (left) = " << Psi.left().qshift() << '\n';
      out << "Left quantum number shift = " << Psi.left_qshift() << '\n';
   }

   if (Psi.right().empty())
   {
      out << "No right strip.\n";
   }
   else
   {
      std::string rw = Psi.get_right_filename();
      if (rw.empty())
         out << "Right semi-infinite strip is stored directly.\n";
      else
         out << "Right semi-infinite strip is stored in the file \"" << rw << "\"\n";
      out << "Right semi-infinite strip unit cell size = " << Psi.right().size() << '\n';

      out << "Right unit cell starting index = " << Psi.right_index() << '\n';
      out << "Quantum number per unit cell (right) = " << Psi.right().qshift() << '\n';
      out << "Right quantum number shift = " << Psi.right_qshift() << '\n';
   }

   out << "Window size = " << Psi.window_size() << '\n';
   out << "ExpIK = " << Psi.exp_ik() << '\n';

   if (Psi.gs_overlap() != 0.0)
      out << "GSOverlap = " << Psi.gs_overlap() << '\n';

   out << std::endl;
}

void ShowBasicInfo(FiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   out << "Finite wavefunction in the left canonical basis.\n";
   out << "Symmetry list = " << Psi.GetSymmetryList() << '\n';
   out << "Length = " << Psi.size() << '\n';
   out << "Quantum number = " << Psi.TransformsAs() << '\n';
   out << "Norm = " << norm_2(Psi) << '\n';
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
          << std::setw(20) << std::string("central-"+Name+"^2");
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
         out << ' ' << std::setw(20) << DM.EvaluateCasimirMoment(i);
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
   void operator()(InfiniteWavefunctionRight const& Psi) const;
   void operator()(IBCWavefunction const& Psi) const;
   void operator()(EAWavefunction const& Psi) const;
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
ShowWavefunction::operator()(InfiniteWavefunctionRight const& Psi) const
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
      ShowStateInfo(Psi.window(), std::cout);

   if (ShowBasis)
      ShowBasisInfo(Psi.window(), std::cout);

   if (ShowEntropy)
      ShowEntropyInfo(Psi.window(), std::cout);

   if (ShowDensity)
      ShowDM(Psi.window(), std::cout);

   if (ShowCasimir)
      ShowCasimirInfo(Psi.window(), std::cout);

   if (ShowLocalBasis)
      ShowLocalBasisInfo(Psi.window(), std::cout);
}

// TODO: We will eventually want to be able to look at each of the windows, but
// it is not currently obvious how to do this.
void
ShowWavefunction::operator()(EAWavefunction const& Psi) const
{
   std::sort(Partition.begin(), Partition.end());
   if (Partition.empty())
   {
      // all partitions in the left boundary
      for (int i = 0; i <= Psi.window_size(); ++i)
         Partition.push_back(i);
   }

   if (ShowStates)
      ShowStateInfo(Psi.window_vec().front(), std::cout);

   if (ShowBasis)
      ShowBasisInfo(Psi.window_vec().front(), std::cout);

   if (ShowEntropy)
      ShowEntropyInfo(Psi.window_vec().front(), std::cout);

   if (ShowDensity)
      ShowDM(Psi.window_vec().front(), std::cout);

   if (ShowCasimir)
      ShowCasimirInfo(Psi.window_vec().front(), std::cout);

   if (ShowLocalBasis)
      ShowLocalBasisInfo(Psi.window_vec().front(), std::cout);
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
         std::cout << "File format version " << Psi->version() << "\n\n";

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
