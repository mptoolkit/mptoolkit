// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/solver-gmres.h"
#include "mp/copyright.h"
#include "common/conflist.h"
#include "mp-algorithms/stateslist.h"
#include "mp-algorithms/dmrgloop.h"
#include "common/terminal.h"
#include "interface/operator-parser.h"
#include <boost/program_options.hpp>
#include <iostream>

namespace prog_opt = boost::program_options;
typedef std::complex<double> complex;

// helper function to load an attribute from the command line or wavefunction
template <typename T>
void LoadAttribute(prog_opt::variables_map& Options, MPWavefunction& Psi, 
                   std::string const& Name, T& Value)
{
   if (Options.count(Name) == 1)
      {
         Value = Options[Name]. template as<T>();
         Psi.Attributes()[Name] = Value;
      }
      else
      {
         if (Psi.Attributes().count(Name) == 0)
         {
            std::ostringstream ostr;
            ostr << "fatal: missing parameter: " << Name << ".";
            throw std::runtime_error(ostr.str());
         }
         Value = Psi.Attributes()[Name].template as<T>();
      }
   std::cout << Name << ": " << Value << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      double GroundstateEnergy = 0;
      double Frequency = 0;
      double Broadening = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(), 
          "operator to use for the Hamiltonian"
          " (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "initial correction vector wavefunction (required)")
         ("lanczos,l", prog_opt::value<std::string>(),
          "fixed Lanczos vector for the right hand side (required)")
         ("config,c", prog_opt::value<std::string>(), "configuration file (required)")
         ("out,o", prog_opt::value<std::string>(), 
          "initial part of filename to use for output files (required)")
         ("GroundstateEnergy,G", prog_opt::value(&GroundstateEnergy),
          "groundstate energy of the Hamiltonian"
          " (wavefunction attribute \"GroundstateEnergy\")")
         ("Frequency,F", prog_opt::value(&Frequency),
          "frequency of the correction vector (wavefunction attribute \"Frequency\")")
         ("Broadening,B", prog_opt::value(&Broadening),
          "broadening constant eta (wavefunction attribute \"Broadening\")")
         ("broadening-propto-energy", prog_opt::bool_switch(),
          "use eta*(H-E) as the broadening term")
         ;

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || !vm.count("wavefunction") || !vm.count("out")
          || !vm.count("config") || !vm.count("lanczos")) 
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-gmres-init [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout << "Starting correction vector GMRES.\n";

      std::string ConfigFile = vm["config"].as<std::string>();
      std::cout << "Configuration file: " << ConfigFile << '\n';
      ConfList Conf(ConfigFile);

      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << '\n';

      std::string LanczosVector = vm["lanczos"].as<std::string>();
      std::cout << "Lanczos vector: " << LanczosVector << '\n';

      std::string BasePathFull = vm["out"].as<std::string>();
      std::cout << "Base filename: " << BasePathFull << '\n';
      std::string BasePath, FileName;
      boost::tie(BasePath, FileName) = pheap::SplitPathFile(BasePathFull);
      if (BasePath.empty()) BasePath = "./";

      long PageSize = Conf.GetBytes("PageSize", 64*1024);
      long PageCacheSize = Conf.GetBytes("PageCacheSize", 16*1024*1024);
      int NumPageFiles = Conf.Get("NumPageFiles", 1);
      std::string BinPath = Conf.Get("BinPath", std::string("."));
      if (BinPath[BinPath.size()-1] != '/') BinPath += '/';
      std::cout << "Storing checkpoint files in directory: " << BinPath << '\n';

      MessageLogger::Logger("PHEAPFS").SetThreshold(Conf.Get("PHeapLogLevel", 0));

      // copy the configuration file
      std::string sstr = "cp -f " + ConfigFile + " " + BasePath + FileName + ".conf";
      system(sstr.c_str());

      // Make the params file
      std::ofstream Params((BasePath + FileName + ".params").c_str());
      std::copy(argv, argv+argc, std::ostream_iterator<char*>(Params, " "));
      Params << '\n';
      Params.close();

      pheap::Initialize(BinPath+FileName+".bin", NumPageFiles, PageSize, PageCacheSize);

      // load the wavefunctions
      MPWavefunction Psi, Lanczos;
      {
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(InputWavefunction);
         Psi = *P;
         P = pheap::ImportHeap(LanczosVector);
         Lanczos = *P;
      }
      
      // load the Hamiltonian
      std::string HamString;
      LoadAttribute(vm, Psi, "Hamiltonian", HamString);
      MPOperator Hamiltonian;
      OperatorList Lattice;
      boost::tie(Lattice, Hamiltonian) = ParseLatticeAndOperator(HamString);

      // Set up the ground state energy, frequency and broadening
      LoadAttribute(vm, Psi, "GroundstateEnergy", GroundstateEnergy);
      LoadAttribute(vm, Psi, "Frequency", Frequency);
      LoadAttribute(vm, Psi, "Broadening", Broadening);
 
      std::string NumStatesStr = Conf.Get("NumStates", "");
      StatesList States(NumStatesStr);
      States.ShowOptions(std::cout);

      MPOperator Identity = Lattice["I"];
      MPOperator A, A2;
      if (vm["broadening-propto-energy"].as<bool>())
      {
         std::cout << "Using broadening proportional to the energy.\n";
	 PANIC("This option no longer works, as the broadening is now handled differently.");
         // For the case of broadening proportional to energy, we have
         // A = E+w-H + i*eta*(H-E)
         // *** To compensate for the bug in GMRES, we take the conjugate here
         complex ScalarPart(GroundstateEnergy+Frequency, GroundstateEnergy*Broadening);
         complex HamPart(-1.0, -Broadening);
         A = ScalarPart*Identity + HamPart*Hamiltonian;

         MPOperator Part = Hamiltonian - GroundstateEnergy*Identity;
         A2 = (1.0 + Broadening*Broadening)*prod(Part, Part, Part.TransformsAs())
            + 2.0*Frequency*Part
            + Frequency*Frequency*Identity;
      }
      else
      {
         // For this case, we use the usual formula
         // A = E+w-H + i*eta
         // *** To compensate for the bug in GMRES, we take the conjugate here
	 //         complex ScalarPart(GroundstateEnergy+Frequency, -Broadening);
         A = (GroundstateEnergy+Frequency)*Identity - Hamiltonian;
         //MPOperator Part = (GroundstateEnergy + Frequency)*Identity - Hamiltonian;
         //A2 = prod(Part, Part, Part.TransformsAs()) 
	 //+ (Broadening*Broadening)*Identity;
	 A2 = A*A;
      }

      SolverGmres solver((CenterWavefunction(Psi)), SplitOperator(A), SplitOperator(A2), 
                         CenterWavefunction(Lanczos), Frequency, Broadening);
      solver.CreateLogFiles(BasePathFull, Conf);   
      pvalue_ptr<DMRGLoop<SolverGmres> > 
         Solver(new DMRGLoop<SolverGmres>(solver, States));

      pheap::ShutdownPersistent(Solver);
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
