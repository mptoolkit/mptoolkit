// -*- C++ -*-

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OutputFile;
      std::string InputFile;
      int Offset = 0;
      int WindowSize = -1;
      bool Force = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("wavefunction,w", prog_opt::value(&InputFile), "input iMPS wavefunction [required]")
	 ("output,o", prog_opt::value(&OutputFile), "output IBC wavefuction [required]")
	 ("force,f", prog_opt::bool_switch(&Force), "overwrite output files")
	 ("size,s", prog_opt::value(&WindowSize), 
	  "initial size (in sites) of the finite window [default is the unit cell size of the input wavefunction]")
	 ("offset", prog_opt::value(&Offset), "site index of the first site of the window [default 0]")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
	 ;

      prog_opt::options_description hidden("Hidden options");

      prog_opt::positional_options_description p;

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || InputFile.empty() || OutputFile.empty())
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] -w <input_psi> -o <output_psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFile);
      InfiniteWavefunctionLeft Psi = InPsi->get<InfiniteWavefunctionLeft>();

      if (WindowSize == -1)
	 WindowSize = Psi.size();

      if (WindowSize % Psi.size() != 0)
      {
	 PANIC("FIXME: currently the window size must be a multiple of the wavefunction unit cell.");
      }

      WavefunctionSectionLeft Window(repeat(Psi, WindowSize / Psi.size()));
      InfiniteWavefunctionRight Right(Psi);
      IBCWavefunction ResultPsi(Psi, Window, Right, Offset);
      ResultPsi.check_structure();

      MPWavefunction Result(ResultPsi);

      // Attributes
      // TODO: add some

      // History log
      Result.AppendHistory(EscapeCommandline(argc, argv));

      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Result));
      pheap::ShutdownPersistent(P);
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
