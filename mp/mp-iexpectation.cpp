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
#include "mps/operator_actions.h"

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
   return delta_shift(E[0], QShift);
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

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("real,r", prog_opt::bool_switch(&ShowReal),
	  "display only the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImag),
	  "display only the imaginary part of the result")
         ("print,p", prog_opt::bool_switch(&Print), "Print the MPO to standard output")
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
         print_copyright(std::cerr);
         std::cerr << "usage: mp-iexpectation [options] <psi> <operator>\n";
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

      UnitCellMPO Op;
      InfiniteLattice Lattice;
      boost::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(OpStr);     

      CHECK(Op.GetSiteList() == Lattice.GetUnitCell().GetSiteList());

      if (Print)
      {
	 print_structure(Op.MPO(), std::cout);
	 //	 std::cout << Op << '\n';
	 //std::cout << "\nTransfer matrix:" << construct_transfer_matrix(herm(GenericMPO(Op.MPO())),
	 //							GenericMPO(Op.MPO())) << '\n';
      };

      // Check that Op is bosonic, otherwise it is not defined
      CHECK(Op.Commute() == LatticeCommute::Bosonic)("Cannot evaluate non-bosonic operator")(Op.Commute());

      LinearWavefunction PsiOrtho = get_orthogonal_wavefunction(*Psi);
      MatrixOperator Rho = scalar_prod(Psi->C_right, herm(Psi->C_right));
      MatrixOperator Identity = MatrixOperator::make_identity(PsiOrtho.Basis1());
      QuantumNumber QShift = Psi->shift();

      int PsiSize = PsiOrtho.size();

      // paranoid check the orthogonalization of the wavefunction
      //DEBUG_CHECK(norm_frob(delta_shift(inject_left(Identity, PsiOrtho, PsiOrtho), QShift) - Identity) < 1E-10)
      //	 (delta_shift(inject_left(Identity, PsiOrtho, PsiOrtho), QShift))(Identity)(norm_frob(delta_shift(inject_left(Identity, PsiOrtho, PsiOrtho), QShift) - Identity));
	 //DEBUG_CHECK(norm_frob(inject_right(Rho, PsiOrtho, PsiOrtho) - delta_shift(Rho, QShift)) < 1E-10);

      // extend Op1 to a multiple of the wavefunction size
      Op.ExtendToCoverUnitCell(PsiOrtho.size());

      // now calculate the actual expectation value
      MatrixOperator X = Identity;
      X = inject_left(X, PsiOrtho, QShift, Op.MPO(), PsiOrtho);

      std::complex<double> x = inner_prod(X, delta_shift(Rho, QShift));

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
