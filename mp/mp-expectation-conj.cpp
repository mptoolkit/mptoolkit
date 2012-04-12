// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "interface/inittemp.h"
#include "interface/operator-parser.h"
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      bool Verbose = false;
      bool NoTempFile = false;
      bool ShowRealPart = false, ShowImagPart = false;
      std::string LhsStr, OpStr, RhsStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("r,real", prog_opt::bool_switch(&ShowRealPart),
	  "display only the real part of the result")
	 ("i,imag", prog_opt::bool_switch(&ShowImagPart),
	  "display only the imaginary part of the result")
         ("notempfile", prog_opt::bool_switch(&NoTempFile),
          "don't use a temporary data file, keep everything in RAM "
          "(faster, but needs enough RAM)")
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "extra debug output")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&LhsStr), "psi1")
         ("operator", prog_opt::value<std::string>(&OpStr), "operator")
         ("rhs", prog_opt::value<std::string>(&RhsStr), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);
      p.add("operator", 1);
      p.add("rhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("operator") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-expectation-conj [options] <psi1> <operator> [<psi2>]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (Verbose)
         std::cout << "Loading LHS wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi1;
      if (NoTempFile)
         Psi1 = pheap::OpenPersistent(LhsStr, mp_pheap::CacheSize(), true);
      else
      {
         mp_pheap::InitializeTempPHeap(Verbose);
         Psi1 = pheap::ImportHeap(LhsStr);
      }

      if (Verbose)
         std::cout << "Parsing operator...\n";

      MPOperator Op = ParseOperator(OpStr);

      if (Verbose)
         std::cout << "Loading RHS wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi2 = RhsStr.empty() ? Psi1 :
         pheap::ImportHeap(RhsStr);

      if (Verbose)
         std::cout << "Calculating expectation value...\n";

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::complex<double> x = expectation_conj(*Psi1, Op, *Psi2);
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
