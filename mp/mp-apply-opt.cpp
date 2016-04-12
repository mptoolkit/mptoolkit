// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "matrixproduct/centerwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/prodoptimizer.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int MinStates = 1;
      int MaxStates = 10000;
      double MinTrunc = 0;
      double EigenCutoff = -1;
      double MixFactor = 0.01;
      bool TwoSite = false;
      int NumSweeps = 2;
      std::string OperatorStr, RhsStr, PsiStr;
      bool ShouldCalculateVariance = true;
      bool NoCalculateVariance = false;
      double ExpectationA2 = -1;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("two-site,2", "modify 2 neighboring sites at once")
	 ("min-states", prog_opt::value<int>(&MinStates), "Minimum number of states to keep [default 1]")
	 ("max-states,m", prog_opt::value<int>(&MaxStates), "Maximum number of states to keep [default 10000]")
         ("trunc,r", prog_opt::value<double>(&MinTrunc), 
          "Cutoff truncation error per site [default 0]")
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff), 
          ("Cutoff threshold for density matrix eigenvalues [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
	 ("mix-factor,f", prog_opt::value<double>(&MixFactor), 
          "Mixing coefficient for the density matrix [default 0.01]")
	 ("sweeps,s", prog_opt::value<int>(&NumSweeps), "Number of half-sweeps to perform [default 2]")
         ("a2", prog_opt::value(&ExpectationA2),
          "Use this pre-computed value of <rhs|adjoint(A)*A|rhs> to speed up the residual calculation")
         ("no-resid", prog_opt::bool_switch(&NoCalculateVariance), 
          "Don't calculate the final residual norm");
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
         std::cerr << "usage: mp-apply-opt [options] <operator> <input> <output>\n";
         std::cerr << "optimizes |output> = operator * |input>\n";
         std::cerr << "where |output> is an existing guess vector\n";
         std::cerr << desc << '\n';
         return 1;
      }

      TwoSite = vm.count("two-site");
      ShouldCalculateVariance = !NoCalculateVariance;
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(PsiStr, CacheSize);
      MPOperator Op = ParseOperator(OperatorStr);
      SplitOperator Operator = Op;
      pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(RhsStr);

      ProductOptimizer Optimizer((CenterWavefunction(*Psi)), Operator, CenterWavefunction(*Rhs));
      Psi = pvalue_ptr<MPWavefunction>();

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = MinTrunc;
      SInfo.EigenvalueCutoff = EigenCutoff;

      double SweepTruncation = 0;

      for (int Sweeps = 0; Sweeps < NumSweeps; ++Sweeps)
      {
         SweepTruncation = 0;
         if (Optimizer.Wavefunction().LeftSize() < Optimizer.Wavefunction().RightSize())
         {
            // sweep right
            Optimizer.ExpandLeft();
            if (TwoSite) Optimizer.ExpandRight();
            Optimizer.Solve();
            TruncationInfo States = Optimizer.TruncateLeft(SInfo, MixFactor);

            // sweep right
            while (Optimizer.RightSize() > 1)
            {
               Optimizer.ShiftRightAndExpand();
               if (TwoSite) Optimizer.ExpandRight();
               double Norm = Optimizer.Solve();
               TruncationInfo States = Optimizer.TruncateLeft(SInfo, MixFactor);
               SweepTruncation += States.TruncationError();

               std::cout << "partition=(" << Optimizer.LeftSize() << ',' << Optimizer.RightSize() << ")"
                         << ", states=" << States.KeptStates()
                         << ", trunc=" << States.TruncationError()
                         << ", largest_discarded_evalue=" << States.LargestDiscardedEigenvalue()
                         << ", norm=" << Norm
                         << '\n';
            }
         }
         else
         {
            Optimizer.ExpandRight();
            if (TwoSite) Optimizer.ExpandLeft();
            Optimizer.Solve();
            Optimizer.TruncateRight(SInfo, MixFactor);

            // sweep left
            while (Optimizer.LeftSize() > 1)
            {
               Optimizer.ShiftLeftAndExpand();
               if (TwoSite) Optimizer.ExpandLeft();
               double Norm = Optimizer.Solve();
               TruncationInfo States = Optimizer.TruncateRight(SInfo, MixFactor);
               SweepTruncation += States.TruncationError();

               std::cout << "partition=(" << Optimizer.LeftSize() << ',' << Optimizer.RightSize() << ")"
                         << ", states=" << States.KeptStates()
                         << ", trunc=" << States.TruncationError()
                         << ", largest_discarded_evalue=" << States.LargestDiscardedEigenvalue()
                         << ", norm=" << Norm
                         << '\n';
            }
         }
         std::cout << "Sweep finished, total truncation = " << SweepTruncation << '\n';
      }
      Psi = new MPWavefunction(Optimizer.Wavefunction().AsLinearWavefunction());

      if (ShouldCalculateVariance)
      {
         // calculate the residual norm 
         // double PR_CosTheta = expectation(*Psi, Op, *Rhs).real();
         double PR_CosTheta = Optimizer.Energy();  // this is faster than expectation()
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
         std::cout << "|Residual|^2 = " << Residual << '\n';
         std::cout << "Theta^2 = " << Theta2 << '\n';
      }
      
      Rhs = pvalue_ptr<MPWavefunction>();  // this might end up in the saved file if we don't do this
      pheap::ShutdownPersistent(Psi);

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
