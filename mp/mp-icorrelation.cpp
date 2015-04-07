// -*- C++ -*- $Id$

#include "mps/infinitewavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/siteoperator-parser.h"
#include "mpo/generic_mpo.h"
#include "mps/operator_actions.h"

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
#include "models/hubbard-u1su2.h"
#include "models/hubbard-u1u1-old.h"
#include "models/boson-u1.h"
#include "models/boson.h"
#include "models/boson-2component-u1z2.h"

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
      std::string StringOpStr = "I";
      bool ShowReal = false, ShowImag = false;
      std::string Model, PsiStr;
      std::string Op1Str, Op2Str;
      int Length = 100;
      bool IncludeFirst = false, Fermionic = false;
      int Verbose = 0;
      int NMax = 5;
      double Spin = 0.5;

      std::cout.precision(14);

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
      else if (Model == "hubbard-u1su2")
      {
	 Site = CreateU1SU2HubbardSite();
      }
      else if (Model == "hubbard-u1u1-old")
      {
	 Site = CreateU1U1HubbardOldOrderingSite();
      }
      else if (Model == "bh-u1")
      {
	 Site = BosonU1(NMax);
      }
      else if (Model == "bh")
      {
	 Site = Boson(NMax);
      }
      else if (Model == "bh2-u1z2")
      {
	 Site = CreateBoseHubbard2BosonsU1Z2Site(NMax);
      }
      else
      {
	 PANIC("mp-icorrelation: fatal: model parameter not recognised.");
      }

      QuantumNumber Ident(Psi->GetSymmetryList());
      
      SimpleOperator Op1 = ParseSiteOperator(Site, Op1Str);
      SimpleOperator Op2 = ParseSiteOperator(Site, Op2Str);
      SimpleOperator Op1Op2 = prod(Op1, Op2, Ident);  // do this now so that we don't include the string term
      SimpleOperator StringOp = ParseSiteOperator(Site, StringOpStr);
      if (IncludeFirst)
      {
	 Op1 = Op1 * StringOp;
      }

      LinearWavefunction PsiOrtho = get_orthogonal_wavefunction(*Psi);
      MatrixOperator Rho = scalar_prod(Psi->C_right, herm(Psi->C_right));
      MatrixOperator Identity = MatrixOperator::make_identity(PsiOrtho.Basis1());
      QuantumNumber QShift = Psi->shift();

      int PsiSize = PsiOrtho.size();

      // paranoid check the orthogonalization of the wavefunction
      DEBUG_CHECK(norm_frob(delta_shift(inject_left(Identity, PsiOrtho, PsiOrtho), QShift) - Identity) < 1E-10);
      DEBUG_CHECK(norm_frob(inject_right(Rho, PsiOrtho, PsiOrtho) - delta_shift(Rho, QShift)) < 1E-10);

      // Make a GenericMPO for our string operator.  The 'vacuum' basis here is
      // the quantum number of our operator
      BasisList StringMPOBasis(StringOp.GetSymmetryList());
      StringMPOBasis.push_back(adjoint(Op1.TransformsAs()));
      GenericMPO StringMPO(PsiSize, OperatorComponent(StringOp.Basis1(), StringMPOBasis, StringMPOBasis));
      for (int i = 0; i < PsiSize; ++i)
      {
	 StringMPO[i](0,0) = StringOp;
      }

      local_basis_compatible_or_abort(Psi->Psi, StringMPO);

      // for each site in the unit cell, make the right hand operator.
      // We start with Rho (the right-identity-vector), and construct the operator
      // at the i'th site of the unit cell.
      // F[i] is Op2(i)
      std::vector<MatrixOperator> F(PsiSize, Rho);
      LinearWavefunction::const_iterator I = PsiOrtho.end();
      int s2 = PsiSize;
      while (I != PsiOrtho.begin())
      {
	 --I; --s2;
	 for (int i = 0; i < PsiSize; ++i)
	 {
	    if (i < s2)
	       F[i] = operator_prod(*I, F[i], herm(*I));
	    else if (i == s2)
	       F[i] = operator_prod(Op2, *I, F[i], herm(*I));
	    else if (i > s2)
	       F[i] = operator_prod(StringOp, *I, F[i], herm(*I));
	 }
      }
      // Shift the F-matrices back to the Basis2()
      for (int i = 0; i < PsiSize; ++i)
      {
	 F[i] = delta_shift(F[i], adjoint(QShift));
      }

      // now the left hand operators.  We now start from the Identity (left-identity-vector)
      // and similarly construct the Op1 for each site of the unit cell
      // E[i] is Op1(i)
      std::vector<MatrixOperator> E(PsiSize, Identity);
      int s1 = 0;
      I = PsiOrtho.begin();
      while (I != PsiOrtho.end())
      {
	 for (int i = 0; i < PsiSize; ++i)
	 {
	    if (i < s1)
	       E[i] = operator_prod(herm(StringOp), herm(*I), E[i], *I);
	    else if (i == s1)
	       E[i] = operator_prod(herm(Op1), herm(*I), E[i], *I);
	    else if (i > s1)
	       E[i] = operator_prod(herm(*I), E[i], *I);
	 }
	 ++I; ++s1;
      }

      // Start of output
      DisplayHeading(ShowReal, ShowImag);

      // Terms internal to the unit cell
      MatrixOperator Initial = Identity;  // track the identity operator as it passes through the unit cell
      s1 = 0;
      I = PsiOrtho.begin();
      while (I != PsiOrtho.end())
      {
	 // The (s1,s1) component
	 MatrixOperator E_ss = operator_prod(herm(prod(Op1,Op2,Ident)), herm(*I), Initial, *I);
	 // E-matrix for Op1(s1)
	 MatrixOperator ThisE = operator_prod(herm(Op1), herm(*I), Initial, *I);

	 int const s1sRemain = PsiSize-s1-1;
	 std::vector<MatrixOperator> E_s_2(s1sRemain);  // operators Op1(s1) * Op2(s2)
	 // Now loop over remaining sites
	 LinearWavefunction::const_iterator I2 = I;
	 int s2 = s1;
	 ++I2; ++s2;
	 while (I2 != PsiOrtho.end())
	 {
	    // firstly do the (s1,s1) component
	    E_ss = operator_prod(herm(*I2), E_ss, *I2);

	    // now the (s1,s2) terms
	    for (int j = 0; j < s2-s1-1; ++j)
	    {
	       // add identity to the terms we've already constructed
	       E_s_2[j] = operator_prod(herm(*I2), E_s_2[j], *I2);
	    }
	    // construct the next term
	    E_s_2[s2-s1-1] = operator_prod(herm(Op2), herm(*I2), ThisE, *I2, Ident);
	    ThisE = operator_prod(herm(StringOp), herm(*I2), ThisE, *I2); 	    // add a site to ThisE
	    ++I2; ++s2;
	 }

	 // print out the results
	 std::complex<double> x = inner_prod(E_ss, Rho);
	 Display(x, s1, s1, ShowReal, ShowImag);

	 for (int j = 0; j < s1sRemain; ++j)
	 {
	    // s2 = s1+j+1
	    x = inner_prod(E_s_2[j], Rho);
	    Display(x, s1, s1+j+1, ShowReal, ShowImag);
	 }

	 // Update the Initializer identity matrix
	 Initial = operator_prod(herm(*I), Initial, *I);
	 ++I; ++s1;
      }

      // terms between unit cells
      int NCell = 1; // number of cells
      while (NCell*PsiSize < Length)
      {
	 for (int ei = 0; ei < PsiSize; ++ei)
	 {
	    for (int fi = 0; fi < PsiSize; ++fi)
	    {
	       std::complex<double> x = inner_prod(E[ei], F[fi]);
	       Display(x, ei, NCell*PsiSize + fi, ShowReal, ShowImag);
	    }
	 }

	 // next unit cell
	 for (int i = 0; i < PsiSize; ++i)
	 {
	    E[i] = inject_left(delta_shift(E[i], Psi->shift()), PsiOrtho, StringMPO, PsiOrtho);
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
