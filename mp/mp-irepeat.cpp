// -*- C++ -*- $Id$

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
      int Times = 2;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("times", prog_opt::value(&Times),
	  "repeat the unit cell this many times [default 2]")
	 ("force,f", prog_opt::bool_switch(&Force),
	  "allow overwriting the output file, if it already exists")
	 ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
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
         print_copyright(std::cerr, "tools", "mp-irepeat");
         std::cerr << "usage: " << basename(argv[0]) << " [options] <operator> <input-psi> <output-psi>\n";
         std::cerr << desc << '\n';

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

      Psi = repeat(Psi, Times);

      PsiPtr.mutate()->Wavefunction() = Psi;
      PsiPtr.mutate()->AppendHistory(EscapeCommandline(argc, argv));
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
