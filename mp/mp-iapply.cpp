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
      std::string InPsi;
      std::string OutPsi;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("op", prog_opt::value(&OpStr), "op")
         ("psi1", prog_opt::value(&InPsi), "psi1")
         ("psi2", prog_opt::value(&OutPsi), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("op", 1);
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
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] <operator> <input-psi> <output-psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (Verbose > 0)
	 std::cout << "Loading wavefunction..." << std::endl;

      pvalue_ptr<MPWavefunction> PsiPtr;
      long CacheSize = getenv_or_default("MP_CACHESIZE", DEFAULT_PAGE_CACHE_SIZE);
      if (InPsi == OutPsi)
	 PsiPtr = pheap::OpenPersistent(InPsi.c_str(), CacheSize);
      else
      {
	 int PageSize = getenv_or_default("MP_PAGESIZE", DEFAULT_PAGE_SIZE);
	 pheap::Initialize(OutPsi, 1, PageSize, CacheSize);
	 PsiPtr = pheap::ImportHeap(InPsi.c_str());
      }

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

      InfiniteLattice Lattice;
      ProductMPO StringOp;
      boost::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(OpStr);
      
      if (Psi.size() != StringOp.size())
      {
	 int Size = statistics::lcm(Psi.size(), StringOp.size());
	 if (Psi.size() < Size)
	 {
	    std::cerr << "mp-iapply: warning: extending wavefunction size to lowest common multiple, which is " << Size << " sites.\n";
	 }
	 Psi = repeat(Psi, Size / Psi.size());
	 StringOp = repeat(StringOp, Size / StringOp.size());
      }

      LinearWavefunction PsiL = get_left_canonical(Psi).first;

      if (Verbose > 0)
	 std::cout << "Applying operator..." << std::endl;

      // apply the string operator
      ProductMPO::const_iterator MI = StringOp.begin();
      for (LinearWavefunction::iterator I = PsiL.begin(); I != PsiL.end(); ++I, ++MI)
      {
	 (*I) = aux_tensor_prod(*MI, *I);
      }

      if (Verbose > 0)
	 std::cout << "Constructing canonical form wavefunction..." << std::endl;

      PsiPtr.mutate()->Wavefunction() = InfiniteWavefunctionLeft(PsiL, Psi.qshift(), Verbose);
      PsiPtr.mutate()->AppendHistory(EscapeCommandline(argc, argv));

      if (Verbose > 0)
	 std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
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
