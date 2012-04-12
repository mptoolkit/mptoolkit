// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/dmrg.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <common/prog_options.h>
#include <iostream>
#include "common/environment.h"

namespace prog_opt = boost::program_options;

void SweepRight(DMRG& dmrg, bool TwoSite, int NumIter, StatesInfo const& SInfo,
                double MixFactor)
{
   double SweepTruncation = 0;
   dmrg.StartSweep();
   dmrg.ExpandLeft();
   if (TwoSite) dmrg.ExpandRight();
   double E = dmrg.Solve(NumIter);
   TruncationInfo States = dmrg.TruncateLeft(SInfo, MixFactor);
   std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize()
	     << ") Energy:" << E 
	     << " States:" << States.KeptStates() 
	     << " Truncation-error:" << States.TruncationError() << '\n';
   SweepTruncation += States.TruncationError();
   // sweep right
   while (dmrg.RightSize() > 1)
   {
      dmrg.ShiftRightAndExpand();
      if (TwoSite) dmrg.ExpandRight();
      E = dmrg.Solve(NumIter);
      States = dmrg.TruncateLeft(SInfo, MixFactor);
      std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize()
		<< ") Energy:" << E
		<< " States:" << States.KeptStates() 
		<< " Truncation-error:" << States.TruncationError() << '\n';
      SweepTruncation += States.TruncationError();
   }
   std::cout << "Cumumative truncation error for sweep: " << SweepTruncation << '\n';
}

void SweepLeft(DMRG& dmrg, bool TwoSite, int NumIter, StatesInfo const& SInfo, double MixFactor)
{
   double SweepTruncation = 0;
   dmrg.StartSweep();
   dmrg.ExpandRight();
   if (TwoSite) dmrg.ExpandLeft();
   double E = dmrg.Solve(NumIter);
   TruncationInfo States = dmrg.TruncateRight(SInfo, MixFactor);
   std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize()
	     << ") Energy:" << E
	     << " States:" << States.KeptStates() 
	     << " Truncation-error:" << States.TruncationError() << '\n';
   SweepTruncation += States.TruncationError();

   // sweep left
   while (dmrg.LeftSize() > 1)
   {
      dmrg.ShiftLeftAndExpand();
      if (TwoSite) dmrg.ExpandLeft();
      E = dmrg.Solve(NumIter);
      States = dmrg.TruncateRight(SInfo, MixFactor);
      std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize()
		<< ") Energy:" << E
		<< " States:" << States.KeptStates() 
		<< " Truncation-error:" << States.TruncationError() << '\n';
      SweepTruncation += States.TruncationError();
   }
   std::cout << "Cumumative truncation error for sweep: " << SweepTruncation << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      int NumIter = 4;
      int MinStates = 1;
      int MaxStates = 100000;
      double MixFactor = 0.01;
      bool TwoSite = false;
      int NumSweeps = 2;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      bool NoVariance = false;
      bool UseDGKS = false;
      std::string Solver = "lanczos";
      
      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(), 
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "wavefunction to apply DMRG (required)")
	 ("two-site,2", "modify 2 neighboring sites at once (traditional DMRG)")
	 ("iter,i", prog_opt::value(&NumIter), 
          FormatDefault("Number of Lanczos iterations per step", NumIter).c_str())
	 ("max-states,m", prog_opt::value(&MaxStates), 
          FormatDefault("Maximum number of states to keep", MaxStates).c_str())
         ("min-states", prog_opt::value(&MinStates), 
          FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("trunc,r", prog_opt::value(&TruncCutoff), 
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff), 
          FormatDefault("Cutoff threshold for density matrix eigenvalues (alternative to truncation error)",
			EigenCutoff).c_str())
	 ("mix-factor,f", prog_opt::value(&MixFactor), 
          FormatDefault("Mixing coefficient for the density matrix", MixFactor).c_str())
	 ("sweeps,s", prog_opt::value(&NumSweeps), 
          FormatDefault("Number of half-sweeps to perform", NumSweeps).c_str())
         ("Solver,S", prog_opt::value(&Solver), 
          FormatDefault("Eigensoler to use ('lanczos', 'arnoldi', 'davidson')", Solver).c_str())
         ("orthogonal", prog_opt::value<std::vector<std::string> >(), 
          "force the wavefunction to be orthogonal to this state")
	 ("no-variance", prog_opt::bool_switch(&NoVariance), "Don't calculate the variance")
	 ("dgks", prog_opt::bool_switch(&UseDGKS), "Use DGKS correction for the orthogonality vectors")
	  ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") || vm.count("wavefunction") == 0) 
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-dmrg [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout << "Starting DMRG...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

      // Open the wavefunction
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> P = pheap::OpenPersistent(InputWavefunction, CacheSize);
      CenterWavefunction Psi = *P;
      P = pvalue_ptr<MPWavefunction>();
      // make sure the wavefunction is normalized
      Psi.normalize();

      // Make sure the center matrix is at one edge
      if (Psi.LeftSize() != 1 && Psi.RightSize() != 1)
      {
	 TRACE(Psi.LeftSize())(Psi.RightSize());
	 std::cout << "The center matrix is not located at an edge.  Rotating..." << std::flush;
	 if (Psi.LeftSize() > Psi.RightSize())
	 {
	    while (Psi.RightSize() > 1)
	       Psi.RotateRight();
	 }
	 else
	 {
	    while (Psi.LeftSize() > 1)
	       Psi.RotateLeft();
	 }
	 std::cout << "done" << std::endl;
      }

      // Set up the Hamiltonian
      std::string HamString;
      if (vm.count("Hamiltonian") == 1)
      {
         HamString = vm["Hamiltonian"].as<std::string>();
         Psi.Attributes()["Hamiltonian"] = HamString;
      }
      else
      {
         if (Psi.Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: No Hamiltonian specified.\n"
               "Specify a Hamiltonian either with --Hamiltonian parameter, "
               "or as an attribute of the initial wavefunction.\n";
            return 1;
         }
         HamString = Psi.Attributes()["Hamiltonian"].as<std::string>();
      }
      std::cout << "Hamiltonian: " << HamString << std::endl;
      SplitOperator Hamiltonian = ParseOperator(HamString);

      // Now we can construct the actual DMRG object
      DMRG dmrg(Psi, Hamiltonian);
      Psi = CenterWavefunction(); // we don't need Psi anymore, it will take up space on disk

      // Add the orthogonal states
      dmrg.Wavefunction().Attributes()["OrthogonalSet"] = "";
      if (vm.count("orthogonal") != 0)
      {
         std::vector<std::string> OrthoSet = vm["orthogonal"].as<std::vector<std::string> >();
         for (std::size_t j = 0; j < OrthoSet.size(); ++j)
         {
            pvalue_ptr<MPWavefunction> Ortho = pheap::ImportHeap(OrthoSet[j]);
	    std::cout << "Adding orthogonality constraint for " << OrthoSet[j] << std::endl;
            dmrg.AddOrthogonalState(*Ortho);
            if (!dmrg.Wavefunction().Attributes()["OrthogonalSet"].empty())
               dmrg.Wavefunction().Attributes()["OrthogonalSet"] += " ";
            dmrg.Wavefunction().Attributes()["OrthogonalSet"] += "--orthogonal " 
               + quote_shell(OrthoSet[j]);
         }
      }

      dmrg.UseDGKS = UseDGKS;
      dmrg.Solver = Solver;

      TwoSite = vm.count("two-site");
      if (TwoSite)
	 std::cout << "Optimizing two sites at a time." << std::endl;

      std::cout << "Density matrix mixing coefficient: " << MixFactor << std::endl;
      std::cout << "Number of Lanczos iterations: " << NumIter << std::endl;
      std::cout << "Number of half-sweeps: " << NumSweeps << std::endl;
      if (UseDGKS) std::cout << "Using DGKS correction." << std::endl;
      std::cout << "Using solver: " << Solver << std::endl;

      bool CalculateH2 = !NoVariance;
      std::cout << "Will ";
      if (!CalculateH2)
	 std::cout << "not ";
      std::cout << "calculate the variance at the end of each sweep.\n";

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      CenterWavefunction OldPsi = dmrg.Wavefunction();

      SplitOperator Ham2 = prod(dmrg.Ham, dmrg.Ham, dmrg.Ident);

      for (int Sweeps = 0; Sweeps < NumSweeps; ++Sweeps)
      {
	 if (dmrg.LeftSize() == 1)
	    SweepRight(dmrg, TwoSite, NumIter, SInfo, MixFactor);
	 else
	    SweepLeft(dmrg, TwoSite, NumIter, SInfo, MixFactor);

	 // the dmrg.Wavefunction() is not normalized anymore
	 double Norm2 = norm_2_sq(dmrg.Wavefunction());

	 // We need to re-calculate the energy, since it will have changed slightly after the truncation 
	 double E = dmrg.Energy()/Norm2;
	 std::cout << "E = " << E << '\n';

	 if (CalculateH2)
	 {
	    double h2 = std::abs(expectation(dmrg.Wavefunction(), Ham2, dmrg.Wavefunction()))/Norm2;
	    double nh2 = h2 - E*E;
	    std::cout << "(H-E)^2 = " << nh2 << '\n';
	 }

	 Psi = dmrg.Wavefunction();
	 double Overlap = dmrg.FidelityLoss();
	 
	 std::cout << "Wavefunction difference from last half-sweep = " << Overlap << '\n';
	 OldPsi = Psi;
      }

      std::cout << "Finished." << std::endl;
      P = pvalue_ptr<MPWavefunction>(new MPWavefunction(dmrg.Wavefunction().AsLinearWavefunction()));
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
