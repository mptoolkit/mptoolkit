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
#include "lattice/unitcell-parser.h"
#include "mpo/finite_mpo.h"
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

// inject_left for a FiniteMPO.  This can have support on multiple wavefunction unit cells
MatrixOperator
inject_left(MatrixOperator const& m, 
            LinearWavefunction const& Psi1,
	    QuantumNumbers::QuantumNumber const& QShift,
            FiniteMPO const& Op, 
            LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   DEBUG_CHECK_EQUAL(m.Basis1(), Psi1.Basis1());
   DEBUG_CHECK_EQUAL(m.Basis2(), Psi2.Basis1());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   MatrixOperator Result = m;
   StateComponent E(Op.Basis1(), m.Basis1(), m.Basis2());
   E[0] = m;
   E.debug_check_structure();
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   FiniteMPO::const_iterator OpIter = Op.begin();
   while (OpIter != Op.end())
   {
      if (I1 == Psi1.end())
      {
	 I1 = Psi1.begin();
	 I2 = Psi2.begin();
	 E = delta_shift(E, QShift);
      }
      E = operator_prod(herm(*OpIter), herm(*I1), E, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return E[0];
}

int main(int argc, char** argv)
{
   try
   {
      bool ShowReal = false, ShowImag = false;
      std::string Model, PsiStr;
      std::string OpStr;
      int UnitCellSize = 1;
      int Verbose = 0;
      int NMax = 5;
      double Spin = 0.5;
      bool Print = false;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("real,r", prog_opt::bool_switch(&ShowReal),
	  "display only the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImag),
	  "display only the imaginary part of the result")
	 ("unitcell,u", prog_opt::value(&UnitCellSize),
	  "Size of the Hamiltonian unit cell")
         ("spin", prog_opt::value(&Spin),
          "spin (for su(2) models)")
         ("print,p", prog_opt::bool_switch(&Print), "Print the MPO to standard output")
	 ("nmax", prog_opt::value(&NMax),
	  FormatDefault("Maximum number of bosons per site [for bosonic models]", NMax).c_str())
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "Verbose output (use multiple times for more output)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("model", prog_opt::value(&Model), "model")
         ("psi", prog_opt::value(&PsiStr), "psi")
         ("op", prog_opt::value(&OpStr), "op")
         ;

      prog_opt::positional_options_description p;
      p.add("model", 1);
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
         print_copyright(std::cerr);
         std::cerr << "usage: mp-iexpectation [options] <model> <psi> <operator>\n";
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
      else if (Model == "hubbard-u1su2")
      {
	 Site = CreateU1SU2HubbardSite();
      }
      else if (Model == "hubbard-u1u1")
      {
	 Site = CreateU1U1HubbardSite();
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
	 PANIC("mp-iexpectation: fatal: model parameter unrecognised.");
      }

      QuantumNumber Ident(Psi->GetSymmetryList());
      
      UnitCell Cell = repeat(Site, UnitCellSize);

      FiniteMPO Op = ParseUnitCellOperator(Cell, 0, OpStr);

      if (Print)
      {
	 std::cout << Op << '\n';
      };

      // Check that Op is bosonic, otherwise it is not defined
      CHECK(Op.Commute() == LatticeCommute::Bosonic)("Cannot evaluate non-bosonic operator")(Op.Commute());

      LinearWavefunction PsiOrtho = get_orthogonal_wavefunction(*Psi);
      MatrixOperator Rho = scalar_prod(Psi->C_right, herm(Psi->C_right));
      MatrixOperator Identity = MatrixOperator::make_identity(PsiOrtho.Basis1());
      QuantumNumber QShift = Psi->shift();

      int PsiSize = PsiOrtho.size();

      // paranoid check the orthogonalization of the wavefunction
      DEBUG_CHECK(norm_frob(delta_shift(inject_left(Identity, PsiOrtho, PsiOrtho), QShift) - Identity) < 1E-10);
      DEBUG_CHECK(norm_frob(inject_right(Rho, PsiOrtho, PsiOrtho) - delta_shift(Rho, QShift)) < 1E-10);

      // extend Op1 to a multiple of the wavefunction size
      while (Op.size() % PsiOrtho.size() != 0)
      {
	 Op = join(Op, identity_mpo(Cell));
      }

      // now calculate the actual expectation value
      MatrixOperator X = Identity;
      X = inject_left(X, PsiOrtho, QShift, Op, PsiOrtho);

      std::complex<double> x = inner_prod(X, Rho);

      if (ShowReal)
	 std::cout << x.real() << "   ";
      if (ShowImag)
	 std::cout << x.imag();
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
