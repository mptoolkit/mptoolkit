// -*- C++ -*-

#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "tensor/wigner_eckart.h"
#include "tensor/regularize.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

// functor to use the visitor pattern with wavefunction types
struct ApplyConj : public boost::static_visitor<void>
{
   template <typename T>
   void operator()(T& Psi) const
   {
      inplace_conj(Psi);
   }
};

int main(int argc, char** argv)
{
   try
   {
      bool Force = false;
      std::string SList;
      std::string InputFile;
      std::string OutputFile;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("force,f", prog_opt::bool_switch(&Force),
	  "overwrite the output file, if it exists")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("inpsi", prog_opt::value(&InputFile), "psi1")
         ("outpsi", prog_opt::value(&OutputFile), "psi1")
         ;

      prog_opt::positional_options_description p;
      p.add("inpsi", 1);
      p.add("outpsi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("inpsi") < 1)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] <input-psi> [output-psi]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Complex conjugation of a wavefunction.\n";
         return 1;
      }

      pvalue_ptr<MPWavefunction> Psi;

      if (OutputFile.empty())
      {
	 // re-use the input file as the output file
	 Psi = pheap::OpenPersistent(InputFile, mp_pheap::CacheSize());
      }
      else
      {
	 // create a new file for output
	 pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
	 // and load the input wavefunction
	 Psi = pheap::ImportHeap(InputFile);
      }

      boost::apply_visitor(ApplyConj(), Psi.mutate()->Wavefunction());

      Psi.mutate()->AppendHistory(EscapeCommandline(argc, argv));
      Psi.mutate()->SetDefaultAttributes();

      pheap::ShutdownPersistent(Psi);
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
