// -*- C++ -*- $Id$

#include "quantumnumbers/all_symmetries.h"
#include "wavefunction/mpwavefunction.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "common/terminal.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      bool Reverse = false;
      std::string Filename;
      std::string Message;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("message,m", prog_opt::value(&Message), "add a new history entry");
	 ("reverse,r", prog_opt::bool_switch(&Reverse), "reverse order, newest first")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&Filename), "psi")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);
      
      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("psi") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] <wavefunction>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (vm.count("message"))
      {
	 pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(Filename, mp_pheap::CacheSize());
	 Psi.mutate()->AppendHistory(Message);
	 pheap::ShutdownPersistent(Psi);
      }
      else
      {
	 pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(Filename, mp_pheap::CacheSize(), true);
	 if (Reverse)
	 {
	    Psi->History().print_newest_first(std::cout);
	 }
	 else
	 {
	    Psi->History().print(std::cout);
	 }
      }
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
