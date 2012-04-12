// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

void ShowBasicInfo(MPWavefunction const& Psi, std::ostream& out)
{
   out << "Symmetry list is " << Psi.GetSymmetryList() << '\n';
   out << "State transforms as " << Psi.TransformsAs() << '\n';
   out << "Number of sites = " << Psi.size() << '\n';
}

void ShowStateInfo(LinearWavefunction const& Psi, std::ostream& out)
{
   LinearWavefunction::const_iterator I = Psi.begin();
   int Bond = 0;
   out << "#bond    #dimension  #degree\n";
   out << std::setw(5) << Bond << "    "
       << std::setw(10) << I->Basis1().total_dimension() << "  "
       << std::setw(7) <<  I->Basis1().total_degree()
       << '\n';
   ++Bond;
   while (I != Psi.end())
   {
      out << std::setw(5) << Bond << "    "
          << std::setw(10) << I->Basis2().total_dimension() << "  "
          << std::setw(7) <<  I->Basis2().total_degree()
          << '\n';
      ++Bond; ++I;
   }
}

void ShowBasisInfo(LinearWavefunction const& Psi, std::ostream& out)
{
   LinearWavefunction::const_iterator I = Psi.begin();
   int Bond = 0;
   out << "Basis at bond " << Bond << ":\n";
   out << I->Basis1() << '\n';
   ++Bond;
   while (I != Psi.end())
   {
      out << "Basis at bond " << Bond << ":\n";
      out << I->Basis1() << '\n';
      ++Bond; ++I;
   }
}
void ShowEntropyInfo(CenterWavefunction Psi, std::ostream& out, bool Base2)
{
   while (Psi.LeftSize() > 1) Psi.RotateLeft();

   out << "#left-size #right-size   #entropy" << (Base2 ? "(base-2)" : "(base-e)") << '\n';

   DensityMatrix<MatrixOperator> DM(scalar_prod(Psi.Center(), herm(Psi.Center())));
   out.precision(16);
   double Entropy = DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2);
   out << std::setw(10) << Psi.LeftSize() << ' ' << std::setw(11) << Psi.RightSize() << ' ' 
       << std::setw(18) << Entropy << '\n';

   while (Psi.RightSize() > 1)
   {
      Psi.RotateRight();
      DM = scalar_prod(Psi.Center(), herm(Psi.Center()));
      Entropy = DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2);
      out << std::setw(10) << Psi.LeftSize() << ' ' << std::setw(11) << Psi.RightSize() << ' ' 
	  << std::setw(18) << Entropy << '\n';
   }
   out.flush();
}

void ShowCasimirInfo(CenterWavefunction Psi, std::ostream& out)
{
   while (Psi.LeftSize() > 1) Psi.RotateLeft();

   SymmetryList const SList = Psi.GetSymmetryList();
   int const NumCasimir = SList.NumCasimirOperators();
   out.precision(16);
   out.setf(std::ios::fixed);

   out << "#left-size #right-size   ";
   for (int i = 0; i < NumCasimir; ++i)
   {
      if (i != 0)
	 out << " ";
      std::string Name = "#" + SList.CasimirName(i);
      out << std::setw(20) << Name;
   }
   out << '\n';

   DensityMatrix<MatrixOperator> DM(scalar_prod(Psi.Center(), herm(Psi.Center())));
   out << std::setw(10) << Psi.LeftSize() << ' ' << std::setw(11) << Psi.RightSize() << "   ";
   for (int i = 0; i < NumCasimir; ++i)
   {
      if (i != 0) 
	 out << ' ';
      out << std::setw(20) << DM.EvaluateCasimir(i);
   }
   out << '\n';

   while (Psi.RightSize() > 1)
   {
      Psi.RotateRight();
      DM = scalar_prod(Psi.Center(), herm(Psi.Center()));
      out << std::setw(10) << Psi.LeftSize() << ' ' << std::setw(11) << Psi.RightSize() << "   ";
      for (int i = 0; i < NumCasimir; ++i)
      {
	 if (i != 0) 
	    out << ' ';
	 out << std::setw(20) << DM.EvaluateCasimir(i);
      }
      out << '\n';
   }
   out.flush();
}

void ShowLocalBasis(CenterWavefunction Psi, std::ostream& out)
{
   for (int i = 0; i < Psi.size(); ++i)
   {
      out << "Basis at site " << i << '\n';
      out << Psi.Lookup(i).SiteBasis() << std::endl;
   }
}

void ShowDM(CenterWavefunction Psi, std::ostream& out, int Max, bool Base2)
{
   DensityMatrix<MatrixOperator> DM(scalar_prod(Psi.Center(), herm(Psi.Center())));

   if (Max > 0)
   {
      out << "\nReduced density matrix at partition (" 
          << Psi.LeftSize() << "," << Psi.RightSize() << ") :\n";
      DM.DensityMatrixReport(out, Max, Base2);
   }
   else
   {
      out << "Entropy at partition ("
          << Psi.LeftSize() << "," << Psi.RightSize() << ") = "
          << DensityEntropy(DM.begin(), DM.end(), DM.EigenSum()) << '\n';
   }
}

void ShowWavefunc(CenterWavefunction Psi, std::ostream& out, bool ShowEntropy, bool ShowStates, int Max,
		  bool Base2)
{
   out << "Symmetry list is " << Psi.GetSymmetryList() << '\n';
   out << "State transforms as " << Psi.TransformsAs() << '\n';
   out << "Number of sites = " << Psi.size() << '\n';

   if (Max == -1 && !ShowEntropy && !ShowStates) return;

   if (Max == 0)
      Max = 100000;

   //   out << "Left-most basis:\n" << Psi.Left().Basis1();
   while (Psi.LeftSize() > 1) Psi.RotateLeft();

   if (ShowStates)
   {
      std::cout << "Number of states at partition ("
                << Psi.LeftSize() << "," << Psi.RightSize()
                << ") = " << Psi.Center().Basis1().total_dimension() 
                << " , degree = " << Psi.Center().Basis1().total_degree() << '\n';
   }
   else
      ShowDM(Psi, out, Max, Base2);
   while (Psi.RightSize() > 1)
   {
      Psi.RotateRight();
      if (ShowStates)
      {
         std::cout << "Number of states at partition ("
                   << Psi.LeftSize() << "," << Psi.RightSize()
                   << ") = " << Psi.Center().Basis1().total_dimension() 
                   << " , degree = " << Psi.Center().Basis1().total_degree() << '\n';
      }
      else
         ShowDM(Psi, out, Max, Base2);
   }

#if 0
   DensityMatrix<MatrixOperator> DM(scalar_prod(Psi.Center(), herm(Psi.Center())));
   out << "\nReduced density matrix at partition (" 
       << Psi.LeftSize() << "," << Psi.RightSize() << ") :\n";
   DM.DensityMatrixReport(out, Max, Base2);

   while (Psi.LeftSize() > 1)
   {
      DensityMatrix<MatrixOperator> DM(scalar_prod(herm(Psi.Center()), Psi.Center()));
      out << "\nReduced density matrix at partition (" 
	  << Psi.LeftSize() << "," << Psi.RightSize() << ") :\n";
      DM.DensityMatrixReport(out, Max, Base2);
      Psi.RotateLeft();
   }
#endif
}

int main(int argc, char** argv)
{
   try 
   {
      int MaxEigenvalues = -1;
      bool ShowLocal = false;
      bool Base2 = false;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("entropy,e", "show the entropy at each partition")
         ("states,s", "show the number of states at each partition")
	 ("basis,a", "show the complete basis at each partition")
         ("density-matrix,d", "show the density matrix eigenvalues")
	 ("casimir,c", "show the values of the casimir invariant operators at each partition")
         ("limit,l", prog_opt::value<int>(&MaxEigenvalues), 
          "limit the density matrix display to N eigenvalues (implies --density-matrix)")
         ("localbasis,b", prog_opt::bool_switch(&ShowLocal),
          "Show the local basis at each site")
	 ("base2,2", prog_opt::bool_switch(&Base2), "show the entropy using base 2 instead of base e")
         ;
      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("input-wavefunction", prog_opt::value<std::string>(), "input wavefunction (required)")
         ;

      prog_opt::positional_options_description p;
      p.add("input-wavefunction", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || vm.count("input-wavefunction") == 0) 
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-info [options] input-wavefunction\n";
         std::cerr << desc << "\n";
         return 1;
      }
      
      bool ShowEntropy = false;
      if (vm.count("entropy"))
         ShowEntropy = true;

      if (vm.count("density-matrix"))
         MaxEigenvalues = std::max(MaxEigenvalues, 0);

      bool ShowStates = vm.count("states");
      bool ShowCasimir = vm.count("casimir");
      bool ShowBasis = vm.count("basis");

      std::string Wavefunc = vm["input-wavefunction"].as<std::string>();

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(Wavefunc, CacheSize, true);

      if (ShowEntropy)
      {
	 ShowEntropyInfo(CenterWavefunction(*Psi), std::cout, Base2);
      }
      else if (ShowStates)
      {
         ShowStateInfo(*Psi, std::cout);
      }
      else if (ShowBasis)
      {
	 ShowBasisInfo(*Psi, std::cout);
      }
      else if (ShowCasimir)
      {
	 ShowCasimirInfo(CenterWavefunction(*Psi), std::cout);
      }
      else if (ShowLocal)
         ShowLocalBasis(*Psi, std::cout);
      else
         ShowWavefunc(CenterWavefunction(*Psi), std::cout, ShowEntropy, ShowStates, MaxEigenvalues, Base2);
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
