// -*- C++ -*- $Id$

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
std::vector<int> Partition;

void ShowBasicInfo(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   out << "Wavefunction is an InfiniteWavefunction in the left canonical basis.\n";
   out << "Symmetry list is " << Psi.GetSymmetryList() << '\n';
   out << "Unit cell size = " << Psi.size() << '\n';
   out << "Quantum number per unit cell is " << Psi.qshift() << '\n';
   out << "Number of states = " << Psi.Basis1().total_dimension() << '\n';
   out << std::endl;
}

void ShowStateInfo(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   out << "#partition  #dimension  #degree\n";

   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      VectorBasis B;
      if (*I == 0)
	 B = Psi.Basis1();
      else
	 B = Psi.lambda(*I-1).Basis1();

      out << std::setw(10) << (*I) << "  "
	  << std::setw(10) << B.total_dimension() << "  "
	  << std::setw(7) <<  B.total_degree()
	  << '\n';
   }
   out << std::endl;
}

void ShowBasisInfo(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      VectorBasis B;
      if (*I == 0)
	 B = Psi.Basis1();
      else
	 B = Psi.lambda(*I-1).Basis1();

      out << "Basis at partition " << (*I) << ":\n" << B << '\n';
   }
   out << std::endl;
}

void ShowEntropyInfo(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   out << "#left-size #right-size   #entropy" << (Base2 ? "(base-2)" : "(base-e)") << '\n';
   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      RealDiagonalOperator Lambda;
      if (*I == 0)
	 Lambda = Psi.lambda_r();
      else
	 Lambda = Psi.lambda(*I-1);

      MatrixOperator Rho = Lambda;
      Rho = scalar_prod(Rho, herm(Rho));

      DensityMatrix<MatrixOperator> DM(Rho);

      double Entropy = DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2);
      out << std::setw(10) << (*I) << ' ' << std::setw(11) << (Psi.size() - (*I)) << ' ' 
	  << std::setw(18) << Entropy << '\n';
   }
   out << std::endl;
}

void ShowDM(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      RealDiagonalOperator Lambda;
      if (*I == 0)
	 Lambda = Psi.lambda_r();
      else
	 Lambda = Psi.lambda(*I-1);

      MatrixOperator Rho = Lambda;
      Rho = scalar_prod(Rho, herm(Rho));

      DensityMatrix<MatrixOperator> DM(Rho);

      out << "Reduced density matrix at partition (" 
          << (*I) << "," << (Psi.size() - *I) << ") :\n";
      DM.DensityMatrixReport(out, MaxEigenvalues, Base2);
      out << std::endl;
   }
}

void ShowCasimirInfo(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   out << "#left-size #right-size   ";
   QuantumNumbers::SymmetryList SList = Psi.GetSymmetryList();
   int NumCasimir = SList.NumCasimirOperators();
   for (int i = 0; i < NumCasimir; ++i)
   {
      if (i != 0)
	 out << " ";
      std::string Name = "#" + SList.CasimirName(i);
      out << std::setw(20) << Name;
   }
   out << '\n';

   for (std::vector<int>::const_iterator I = Partition.begin(); I != Partition.end(); ++I)
   {
      RealDiagonalOperator Lambda;
      if (*I == 0)
	 Lambda = delta_shift(Psi.lambda_r(), Psi.qshift());
      else
	 Lambda = Psi.lambda(*I-1);

      MatrixOperator Rho = Lambda;
      Rho = scalar_prod(Rho, herm(Rho));

      DensityMatrix<MatrixOperator> DM(Rho);

      out << std::setw(10) << (*I) << ' ' << std::setw(11) << (Psi.size() - *I) << "   ";
      for (int i = 0; i < NumCasimir; ++i)
      {
	 if (i != 0) 
	    out << ' ';
	 out << std::setw(20) << DM.EvaluateCasimir(i);
      }
      out << '\n';
   }
   out << std::endl;
}

void ShowLocalBasisInfo(InfiniteWavefunctionLeft const& Psi, std::ostream& out)
{
   for (int i = 0; i < Psi.size(); ++i)
   {
      out << "Basis at site " << i << '\n';
      out << Psi[i].LocalBasis() << std::endl;
   }
   out << std::endl;
}

struct ShowWavefunction : public boost::static_visitor<void>
{
   void operator()(InfiniteWavefunctionLeft const& Psi) const;

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

   if (ShowBasic)
      ShowBasicInfo(Psi, std::cout);

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
	 ("casimir,c", prog_opt::bool_switch(&ShowCasimir), 
	  "show the values of the casimir invariant operators at each partition")
         ("limit,l", prog_opt::value<int>(&MaxEigenvalues), 
          "limit the density matrix display to N eigenvalues (implies --density-matrix)")
         ("localbasis,b", prog_opt::bool_switch(&ShowLocalBasis),
          "Show the local basis at each site")
	 ("base2,2", prog_opt::bool_switch(&Base2), "show the entropy using base 2 instead of base e")
         ("partition,p", prog_opt::value(&Partition), 
          "show quantities only for this parition (zero-based, can be used more than once; use --partition 0 to show only quantities at the edge of the unit cell")
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
      
      if (vm.count("help") || vm.count("input-wavefunction") == 0) 
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] input-wavefunction\n";
         std::cerr << desc << "\n";
         return 1;
      }
      
      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (vm.count("limit"))
	 ShowDensity = true;

      if (!ShowEntropy && !ShowStates && !ShowLocalBasis && !ShowBasis && !ShowCasimir && !ShowDensity)
	 ShowBasic = true;

      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(WavefuncFile, mp_pheap::CacheSize(), true);

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
