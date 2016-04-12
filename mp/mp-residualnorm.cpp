// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/environment.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      std::string OperatorStr, RhsStr, PsiStr;
      double ExpectationA2 = -1;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("a2", prog_opt::value(&ExpectationA2),
          "Use this pre-computed value of <rhs|adjoint(A)*A|rhs> to speed up the residual calculation")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("operator", prog_opt::value<std::string>(&OperatorStr), "input operator")
         ("rhs", prog_opt::value<std::string>(&RhsStr), "input wavefunction")
         ("psi", prog_opt::value<std::string>(&PsiStr), "output wavefunction")
         ;

      prog_opt::positional_options_description p;
      p.add("operator", 1);
      p.add("rhs", 1);
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-residualnorm [options] <operator> <rhs> <psi>\n";
         std::cerr << "calculates the residual norm ||r|| where r = psi - operator*rhs\n";
         std::cerr << desc << '\n';
         return 1;
      }

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(PsiStr, CacheSize, true);
      MPOperator Op = ParseOperator(OperatorStr);
      pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(RhsStr);

      double PR_CosTheta = expectation(*Psi, Op, *Rhs).real();
      double PsiNorm2 = norm_2_sq(*Psi);
      if (ExpectationA2 == -1)
      {
         ExpectationA2 = expectation(*Rhs, prod(adjoint(Op), 
                                                Op, 
                                                QuantumNumber(Op.GetSymmetryList())), 
                                     *Rhs).real();
         std::cout << "<rhs|adjoint(A)*A|rhs> = " << ExpectationA2 << '\n';
      }
      double OpRhsNorm2 = ExpectationA2;
      double Residual = PsiNorm2 + OpRhsNorm2 - 2.0 * PR_CosTheta;
      double Theta2 = 2.0 * (1.0 - PR_CosTheta / std::sqrt(PsiNorm2 * OpRhsNorm2));
      std::cout.precision(16);
      std::cout << "|Residual|^2 = " << Residual << '\n';
      std::cout << "Theta^2 = " << Theta2 << '\n';
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
   pheap::Shutdown();
}
