// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      bool ShowRealPart = false, ShowImagPart = false;
      std::string LhsStr, RhsStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("real,r", prog_opt::bool_switch(&ShowRealPart),
	  "display only the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImagPart),
	  "display only the imaginary part of the result")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&LhsStr), "psi1")
         ("rhs", prog_opt::value<std::string>(&RhsStr), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);
      p.add("rhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("rhs") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-overlap [options] <psi1> <psi2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> Psi1 = pheap::OpenPersistent(LhsStr, CacheSize, true);
      pvalue_ptr<MPWavefunction> Psi2 = pheap::ImportHeap(RhsStr);

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::complex<double> x = overlap(*Psi1, *Psi2);
      
      if (ShowRealPart || ShowImagPart)
      {
	 if (ShowRealPart)
	 {
	    std::cout << x.real();
	    if (ShowImagPart)
	       std::cout << "   " << x.imag();
	 }
	 else // if we get here then ShowImagPart is true and ShowRealPart is false
	    std::cout << x.imag(); 
      }
      else // default to C++ complex output
	 std::cout << x;

      std::cout << '\n';

      pheap::Shutdown();

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
