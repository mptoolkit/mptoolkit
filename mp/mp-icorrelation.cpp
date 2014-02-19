
#include "mps/infinitewavefunction.h"
#include "siteoperator/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"

#include "tensor/tensor_eigen.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/spin-z2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/kondo-u1u1.h"
#include "models/hubbard-so4.h"
#include "models/hubbard-u1u1.h"
#include "models/bosehubbard-spinless-u1.h"
#include "models/bosehubbard-2bosons-u1z2.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string StringOpStr = "I";
      bool ShowReal = false, ShowImag = false;
      std::string Model, PsiStr;
      std::string Op1Str, Op2Str;
      int Length = 100;
      bool IncludeFirst = false, Fermionic = false;
      int Verbose = 0;
      int NMax = 5;
      double Spin = 0.5;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("real,r", prog_opt::bool_switch(&ShowReal),
	  "display only the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImag),
	  "display only the imaginary part of the result")
         ("string,s", prog_opt::value(&StringOpStr),
          "calculate a string correlation with the given operator S, as operator1 "
	  "\\otimes S \\otimes S .... S \\otimes operator2")
         ("includefirst,d", prog_opt::bool_switch(&IncludeFirst),
          "when calculating a string correlation, apply the string operator also "
          "to the first site, as operator1*S")
         ("spin", prog_opt::value(&Spin),
          "spin (for su(2) models)")
         ("fermionic,f", prog_opt::bool_switch(&Fermionic),
          "take the correlator to be fermionic; equivalent to \"--string P --includefirst\"")
         ("length,n", prog_opt::value(&Length),
          "calculate correlator to this length [default 100]")
	 ("nmax", prog_opt::value(&NMax),
	  FormatDefault("Maximum number of bosons per site [for bosonic models]", NMax).c_str())
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "Verbose output (use multiple times for more output)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("model", prog_opt::value(&Model), "model")
         ("psi", prog_opt::value(&PsiStr), "psi")
         ("op1", prog_opt::value(&Op1Str), "op1")
         ("op2", prog_opt::value(&Op2Str), "op2")
         ;

      prog_opt::positional_options_description p;
      p.add("model", 1);
      p.add("psi", 1);
      p.add("op1", 1);
      p.add("op2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

     if (Fermionic)
      {
         if (IncludeFirst || StringOpStr != "I")
         {
            std::cerr << "mp-icorrelation: error: cannot combine string correlators with --fermionic.\n\n";
            return 1;
         }
         else
         {
            // a fermionic operator is equivalent to the string correlator
            // operator1*P \otimes P \otimes .... \otimes P \otimes operator2
            // which is equivalent to "--includefirst --string P" options.
            IncludeFirst = true;
            StringOpStr = "P";  // TODO: this should actually be the SignOperator() of the site operator at Op1
         }
      }

      if (vm.count("help") > 0 || vm.count("op2") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-icorrelation [options] <model> <psi> <operator1> <operator2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // If no real or imag specifiation is used, show both parts
      if (!ShowReal && !ShowImag)
         ShowReal = ShowImag = true;

      bool FermionicWarning = false;  // set to true if we have already warned the user about 
      // a possible fermionic problem

      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<InfiniteWavefunction> Psi = pheap::ImportHeap(PsiStr);

      LatticeSite Site;
      if (Model == "sf-u1")
      {
	 Site = CreateU1SpinlessFermion();
      }
      else if (Model == "spin")
      {
	 Site = CreateSpinSite(Spin);
      }
      else if (Model == "spin1")
      {
	 Site = CreateSpinSite(1.0);
      }
      else if (Model == "spin-su2")
      {
	 Site = CreateSU2SpinSite(Spin);
      }
      else if (Model == "spin-u1")
      {
	 Site = CreateU1SpinSite(Spin);
      }
      else if (Model == "spin-z2")
      {
	 Site = CreateZ2SpinSite(0.5);
      }
      else if (Model == "tj-u1")
      {
	 Site = CreateU1tJSite();
      }
      else if (Model == "tj-u1su2")
      {
	 Site = CreateU1SU2tJSite();
      }
      else if (Model == "klm-u1su2")
      {
	 Site = CreateU1SU2KondoSite();
      }
      else if (Model == "klm-u1")
      {
	 Site = CreateU1KondoSite();
      }
      else if (Model == "klm-u1u1")
      {
	 Site = CreateU1U1KondoSite();
      }
      else if (Model == "hubbard-so4")
      {
	 Site = CreateSO4HubbardSiteA();
      }
      else if (Model == "hubbard-u1u1")
      {
	 Site = CreateU1U1HubbardSite();
      }
      else if (Model == "bh-u1")
      {
	 Site = CreateBoseHubbardSpinlessU1Site(NMax);
      }
      else if (Model == "bh2-u1z2")
      {
	 Site = CreateBoseHubbard2BosonsU1Z2Site(NMax);
      }
      else
      {
	 PANIC("mp-icorrelation: fatal: model parameter should be one of tj-u1, tj-u1su2, sf-u1, klm-u1su2.");
      }
 
      QuantumNumber Ident(Psi->GetSymmetryList());
 
      SimpleOperator Op1 = Site[Op1Str];
      SimpleOperator Op2 = Site[Op2Str];
      SimpleOperator Op1Op2 = prod(Op1, Op2, Ident);  // do this now so that we don't include the string term
      SimpleOperator StringOp = Site[StringOpStr];
      if (IncludeFirst)
      {
	 Op1 = Op1 * StringOp;
      }

      LinearWavefunction PsiOrtho = get_orthogonal_wavefunction(*Psi);
      MatrixOperator Rho = scalar_prod(Psi->C_right, herm(Psi->C_right));

      std::cout.precision(14);
      
      // for each site in the unit cell, make the right hand operator
      std::list<MatrixOperator> F;
      LinearWavefunction::const_iterator I = PsiOrtho.end();
      while (I != PsiOrtho.begin())
      {
	 --I;
	 // add string operator to each existing operator
	 for (std::list<MatrixOperator>::iterator J = F.begin(); J != F.end(); ++J)
	 {
	    *J = operator_prod(StringOp, *I, *J, herm(*I));
	 }

      // next site
      F.push_back(operator_prod(Op2, *I, herm(*I)));
   }

      // now the left hand operators.   We also do the correlators for the first unit cell
      typedef std::list<std::pair<int, MatrixOperator> > EType;
      typedef std::list<std::pair<std::pair<int, int>, MatrixOperator> > CorrType;
      std::list<std::pair<int, MatrixOperator> > E;
      std::list<std::pair<std::pair<int, int>, MatrixOperator> > Corr;  // first unit cell correlators
      I = PsiOrtho.begin();
      int pLoc = 1;
      while (I != PsiOrtho.end())
      {
	 // add identity operator to each existing first unit cell operator
	 for (CorrType::iterator J = Corr.begin(); J != Corr.end(); ++J)
	 {
	    J->second = operator_prod(herm(*I), J->second, *I);
	 }
	 // Add the zero-point correlator to the unit cell correlations
	 Corr.push_back(std::make_pair(std::make_pair(pLoc, pLoc),
				       operator_prod(herm(Op1Op2), herm(*I), *I)));
	 // add the next site to the first unit cell correlator
	 for (EType::const_iterator J = E.begin(); J != E.end(); ++J)
	 {
	    Corr.push_back(std::make_pair(std::make_pair(J->first, pLoc),
					  operator_prod(herm(Op2), herm(*I), J->second, *I, Ident)));
	 }
	 // add the string operator to each existing operator
	 for (EType::iterator J = E.begin(); J != E.end(); ++J)
	 {
	    J->second = operator_prod(herm(StringOp), herm(*I), J->second, *I);
	 }
	 // next site
	 E.push_back(std::make_pair(pLoc, operator_prod(herm(Op1), herm(*I), *I)));
	 ++I;
	 ++pLoc;
      }

      // Output the heading
      std::cout << "#i #j";
      if (ShowReal)
	 std::cout << " #real";
      if (ShowImag)
	 std::cout << " #imag";
      std::cout << '\n';
      // output the first unit cell
      for (CorrType::const_iterator J = Corr.begin(); J != Corr.end(); ++J)
      {
	 std::complex<double> x = inner_prod(J->second, Rho);
	 std::cout << J->first.first << "    " << J->first.second << "   ";
	 if (ShowReal)
	    std::cout << x.real() << "   ";
	 if (ShowImag)
	    std::cout << x.imag();
	 std::cout << '\n';
      }

      int NCell = 1; // number of cells
      while (NCell * int(PsiOrtho.size()) < Length)
      {
	 for (EType::iterator J = E.begin(); J != E.end(); ++J)
	 {
	    J->second = delta_shift(J->second, Psi->QShift);
	 }

	 Corr.clear();
	 for (LinearWavefunction::const_iterator I = PsiOrtho.begin(); I != PsiOrtho.end(); ++I)
	 {
	    // add identity site to operators already in Corr
	    for (CorrType::iterator J = Corr.begin(); J != Corr.end(); ++J)
	    {
	       J->second = operator_prod(herm(*I), J->second, *I);
	    }
	    // add next site
	    for (EType::const_iterator J = E.begin(); J != E.end(); ++J)
	    {
	       Corr.push_back(std::make_pair(std::make_pair(J->first, pLoc), 
					     operator_prod(herm(Op2), herm(*I), J->second, *I, Ident)));
	    }
	    // add identity operator to E
	    for (EType::iterator J = E.begin(); J != E.end(); ++J)
	    {
	       J->second = operator_prod(herm(*I), J->second, *I);
	    }
	    ++pLoc;
	 }

	 // output
	 for (CorrType::const_iterator J = Corr.begin(); J != Corr.end(); ++J)
	 {
	    std::complex<double> x = inner_prod(J->second, Rho);
	    std::cout << J->first.first << "    " << J->first.second << "   ";
	    if (ShowReal)
	       std::cout << x.real() << "   ";
	    if (ShowImag)
	       std::cout << x.imag();
	    std::cout << '\n';
	 }
	 ++NCell;
      }

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
