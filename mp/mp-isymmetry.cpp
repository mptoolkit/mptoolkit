// -*- C++ -*- $Id$

#include "mps/infinitewavefunction.h"
#include "mps/packunpack.h"
#include "lattice/latticesite.h"
#include "mps/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "lattice/unitcell.h"
#include "lattice/unitcell-parser.h"
#include "linearalgebra/arpack_wrapper.h"
#include <boost/algorithm/string/predicate.hpp>

namespace prog_opt = boost::program_options;

template <typename Func>
struct PackApplyFunc
{
   PackApplyFunc(PackMatrixOperator const& Pack_, Func f_) : Pack(Pack_), f(f_) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = f(x);
      Pack.pack(x, Out);
   } 
   PackMatrixOperator const& Pack;
   Func f;
};

template <typename Func>
PackApplyFunc<Func>
MakePackApplyFunc(PackMatrixOperator const& Pack_, Func f_)
{
   return PackApplyFunc<Func>(Pack_, f_);
}

std::pair<std::complex<double>, MatrixOperator>
get_left_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, 
                     QuantumNumber const& QShift, 
                     FiniteMPO const& StringOp,
                     double tol = 1E-14, int Verbose = 0)
{
   int ncv = 0;
   QuantumNumber Ident(Psi1.GetSymmetryList());
   PackMatrixOperator Pack(Psi1.Basis2(), Psi2.Basis2(), Ident);
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
							    LeftMultiplyString(Psi1, StringOp, Psi2, QShift)),
					  n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   MatrixOperator LeftVector = Pack.unpack(&(OutVec[0]));
      
   return std::make_pair(LeftEigen[0], LeftVector);
}

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      std::string LatticeFile;
      //      double Tol = 1e-10;
      int Verbose = 0;
      bool Quiet = false;
      std::vector<std::string> OperatorStr;
      std::vector<std::string> CommutatorStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("wavefunction,w", prog_opt::value(&PsiStr), "Wavefunction [required]")
	 ("lattice,l", prog_opt::value(&LatticeFile), "use this lattice file for the operators")
	 ("commutator,c", prog_opt::value(&CommutatorStr), 
	  "calculate the commutator phase angle, U X X^\\dagger = exp(i*theta) X")
	 ("quiet,q", prog_opt::bool_switch(&Quiet), "suppress informational preamble about each operator")
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("operator", prog_opt::value(&OperatorStr), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("operator", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1)
      {
         std::cerr << "usage: " << basename(argv[0]) << " -w <psi> [options] Operator1 [Operator2] ...\n";
         std::cerr << "If -l [--lattice] is specified, then the operators must all come from the specified lattice file\n";
         std::cerr << "Otherwise all operators must be of the form lattice:operator\n";
	 std::cerr << "Calculates the commutator phase of operator pairs <X Y X\u2020 Y\u2020>\n";
	 std::cerr << "For complex conjugation, prefix the operator expression with c&\n";
	 std::cerr << "For spatial reflection, prefix with r& (cr& or rc& for conjugate-reflection)\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // if the number of eigenvalues is specified but
      // the cutoff is not, then set the cutoff to zero
      //      if (vm.count("num-eigenvalues") == 1 && vm.count("eigen-cutoff") == 0)
      //      {
      //         EigenCutoff = 0;
      //      }


      // Load the wavefunction
      pvalue_ptr<InfiniteWavefunction> Psi 
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      // Load the lattice, if it was specified
      pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

      UnitCell Cell = Lattice->GetUnitCell();
      int UnitCellSize = Cell.size();
      int const NumUnitCells = Psi->size() / UnitCellSize;

      // orthogonalize the wavefunction
      LinearWavefunction Psi1 = get_orthogonal_wavefunction(*Psi);
      MatrixOperator Rho = scalar_prod(Psi->C_right, herm(Psi->C_right));
      MatrixOperator Identity = MatrixOperator::make_identity(Psi1.Basis1());
      double Dim = Psi1.Basis1().total_degree();

      // reflected and conjugated versions of the wavefunction - we leave them as null
      // until they are actually needed
      LinearWavefunction PsiR, PsiC, PsiRC;

      // The list of U matrices, for each operator
      std::vector<MatrixOperator> U;

      for (unsigned i = 0; i < OperatorStr.size(); ++i)
      {
	 std::string OpStr = OperatorStr[i];
	 bool Reflect = false;
	 bool Conjugate = false;
	 LinearWavefunction* Psi2 = &Psi1;
	 // Do we have time reversal or reflection?
	 if (boost::starts_with(OpStr, "r&"))
	 {
	    Reflect = true;
	    OpStr = std::string(OpStr.begin()+2, OpStr.end());
	    if (PsiR.empty())
	    {
	       InfiniteWavefunction PR = reflect(*Psi);
	       orthogonalize(PR);
	       PsiR = get_orthogonal_wavefunction(PR);
	       //TRACE(PR.C_right)(Psi->C_right);
	    }
	    Psi2 = &PsiR;
	 }
	 else if (boost::starts_with(OpStr,"c&"))
	 {
	    Conjugate = true;
	    OpStr = std::string(OpStr.begin()+2, OpStr.end());
	    if (PsiC.empty())
	    {
	       PsiC = conj(Psi1);
	    }
	    Psi2 = &PsiC;
	 }
	 else if (boost::starts_with(OpStr,"rc&") || boost::starts_with(OpStr,"cr&"))
	 {
	    Reflect = true;
	    Conjugate = true;
	    OpStr = std::string(OpStr.begin()+3, OpStr.end());
	    if (PsiR.empty())
	    {
	       InfiniteWavefunction PR = reflect(*Psi);
	       orthogonalize(PR);
	       PsiR = get_orthogonal_wavefunction(PR);
	    }
	    if (PsiRC.empty())
	    {
	       PsiRC = conj(PsiR);
	    }
	    Psi2 = &PsiRC;
	 }

	 FiniteMPO StringOperator = ParseUnitCellOperator(Cell, 0, OpStr).MPO();

	 StringOperator = repeat(StringOperator, Psi1.size() / StringOperator.size());

         std::complex<double> e;
         MatrixOperator v;
         boost::tie(e, v) = get_left_eigenvector(Psi1, *Psi2, Psi->shift(), StringOperator);

         v *= std::sqrt(Dim); // make it properly unitary
         U.push_back(v);

	 //	 TRACE(v); //(scalar_prod(herm(v),v));
	 TRACE(SingularValues(v));

	 TRACE(inner_prod(v, conj(v)));

#if 0
	 // Make v the closest approximation to a unitary
	 MatrixOperator U, D, Vh;
	 SingularValueDecomposition(v, U, D, Vh);
	 v = U*Vh;
#endif

	 // check if v is unitary; this will be close to zero
	 double u = norm_frob(scalar_prod(herm(v), v) - Identity);

	 if (!Quiet)
	    std::cout << "Operator: " << OperatorStr[i] 
		      << " eigenvalue=" << e
		      << " expectation=" << trace(v*Rho)
		      << " unitary=" << u << std::endl;
	 
	 // The eigenvalue should be nearly 1, or it isn't a unitary operator
	 if (LinearAlgebra::norm_frob(LinearAlgebra::norm_frob(e)-1.0) > 1E-5)
	 {
	    WARNING("Operator is not unitary - eigenvalue is not modulus 1!")(OperatorStr[i])(e);
	 }

	 // And the eigenvector should be a unitary matrix
	 if (u > 1E-5)
	 {
	    WARNING("Operator is not unitary - U U\\dagger is not identity!")(OperatorStr[i])(u);
	 }
      }
      if (!Quiet)
	 std::cout << '\n';

      // Now go through each operator pair
      std::cout << "#Op1                 #Op2                 #Commutator-Real  #Commutator-Imag\n";
      for (unsigned i = 0; i < U.size(); ++i)
      {
         for (unsigned j = i+1; j < U.size(); ++j)
         {
	    //            std::complex<double> x = inner_prod(U[i]*U[j], U[j]*U[i]) / Dim;
            std::complex<double> x = inner_prod(scalar_prod(herm(U[i]*U[j]), U[j]*U[i]), Rho);
	    //            std::complex<double> tr = trace(U[i]*U[j]*Rho);
	    //	    TRACE(scalar_prod(herm(U[i]*U[j]), U[i]*U[j]));
            std::cout << std::setw(20) << std::left << OperatorStr[i] << " "
                      << std::setw(20) << std::left << OperatorStr[j] << " "
                      << std::setw(17) << std::right << std::fixed << x.real() << " "
                      << std::setw(17) << std::right << std::fixed << x.imag() << " "
		      << std::endl;
         }
      }

#if 0
      // Now calculate the commutators
      for (unsigned j = 0; j < CommutatorStr.size(); ++j)
      {
	 FiniteMPO Op = ParseUnitCellOperator(Cell, NumUnitCells, CommutatorStr[j]).MPO();
	 Op = repeat(Op, Psi.size() / UnitCellSize);
	 MatrixOperator O = MatrixOperator::make_identity(Psi.Basis2());
	 O = inject_left_qshift(O, Op, Psi, Psi1->shift());

	 for (unsigned i = 0; i < U.size(); ++i)
	 {

	    MatrixOperator Obar = triple_prod(U[i], O, herm(U[i]));
	    std::complex<double> Phase = inner_prod(O, Obar) / inner_prod(O, O);
	    std::cout << "Phase " << CommutatorStr[j] << " with " << OperatorStr[i] << " = " << Phase << '\n';
	 }
      }
#endif

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
