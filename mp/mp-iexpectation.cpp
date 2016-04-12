// -*- C++ -*- $Id$

#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/unitcell-parser.h"
#include "wavefunction/operator_actions.h"

namespace prog_opt = boost::program_options;

void DisplayHeading(bool ShowReal, bool ShowImag)
{
   // Output the heading
   std::cout << "#i #j";
   if (ShowReal)
      std::cout << " #real";
   if (ShowImag)
      std::cout << " #imag";
   std::cout << '\n';
}

void
Display(std::complex<double> x, int s1, int s2, bool ShowReal, bool ShowImag)
{
   std::cout << s1 << "    " << s2 << "   ";
   if (ShowReal)
      std::cout << x.real() << "   ";
   if (ShowImag)
      std::cout << x.imag();
   std::cout << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      bool ShowReal = false, ShowImag = false;
      std::string PsiStr;
      std::string OpStr;
      int Verbose = 0;
      bool Print = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("real,r", prog_opt::bool_switch(&ShowReal),
	  "display only the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImag),
	  "display only the imaginary part of the result")
         ("print,p", prog_opt::bool_switch(&Print), "Print the MPO to standard output (use --verbose to see more detail)")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "Verbose output (use multiple times for more output)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&PsiStr), "psi")
         ("op", prog_opt::value(&OpStr), "op")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("op", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("op") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // If no real or imag specifiation is used, show both parts
      if (!ShowReal && !ShowImag)
         ShowReal = ShowImag = true;

      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(PsiStr);

      if (!PsiPtr->is<InfiniteWavefunctionLeft>())
      {
	 std::cerr << "fatal: wavefunction is not an iMPS!\n";
	 exit(1);
      }

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

      UnitCellMPO Op;
      InfiniteLattice Lattice;
      boost::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(OpStr);     

      CHECK(Op.GetSiteList() == Lattice.GetUnitCell().GetSiteList());

      if (Print)
      {
	 print_structure(Op.MPO(), std::cout);
	 if (Verbose > 0)
	 {
	    std::cout << Op.MPO() << '\n';
	 }
	 if (Verbose > 1)
	 {
	    SimpleRedOperator x = coarse_grain(Op.MPO());
	    std::cout << x << "\n";
	 }
	 //	 std::cout << Op << '\n';
	 //std::cout << "\nTransfer matrix:" << construct_transfer_matrix(herm(GenericMPO(Op.MPO())),
	 //							GenericMPO(Op.MPO())) << '\n';
      };

      // Check that Op is bosonic, otherwise it is not defined
      CHECK(Op.Commute() == LatticeCommute::Bosonic)("Cannot evaluate non-bosonic operator")(Op.Commute());

      // extend Op1 to a multiple of the wavefunction size
      Op.ExtendToCoverUnitCell(Psi.size());

      std::complex<double> x = expectation(Psi, Op.MPO());

      if (ShowReal)
	 std::cout << x.real() << "   ";
      if (ShowImag)
	 std::cout << x.imag();
      std::cout << '\n';

      pheap::Shutdown();

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
