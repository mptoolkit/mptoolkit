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
#include "lattice/product-parser.h"

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
      LinearWavefunction PsiL = get_left_canonical(Psi).first;

      InfiniteLattice Lattice;
      ProductMPO StringOp;
      boost::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(OpStr);
      
      CHECK(Psi.size() % StringOp.size() == 0)
	 ("Wavefunction size must be a multiple of the string operator size")
	 (Psi.size())(StringOp.size());
      StringOp = repeat(StringOp, Psi.size() / StringOp.size());

      // apply the string operator
      ProductMPO::const_iterator MI = StringOp.begin();
      for (LinearWavefunction::iterator I = PsiL.begin(); I != PsiL.end(); ++I, ++MI)
      {
	 (*I) = aux_tensor_prod(*MI, *I);
      }

      PsiPtr.mutate()->Wavefunction() = InfiniteWavefunctionLeft(PsiL, Psi.qshift());
      PsiPtr.mutate()->AppendHistory(EscapeCommandline(argc, argv));

      pheap::ShutdownPersistent(PsiPtr);
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
