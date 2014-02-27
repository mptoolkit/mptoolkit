
#include "mpo/triangular_mpo.h"
#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"
#include "common/polynomial.h"
#include "tensor/tensor_eigen.h"
#include "common/prog_opt_accum.h"
#include "mp/copyright.h"
#include "common/prog_options.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"
#include "models/bosehubbard-spinless.h"
#include "models/bosehubbard-spinless-u1.h"
#include "models/hubbard-u1su2.h"

#include "mps/momentum_operations.h"
#include "mp-algorithms/triangular_mpo_solver.h"

namespace prog_opt = boost::program_options;

// binomial coefficient n choose k
long Binomial(int n, int k)
{
   if (k > n/2)
      k = n-k;     // take advantage of symmetry
   double r = 1.0;
   for (int i = 1; i <= k; ++i)
   {
      r *= double(n-k+i) / double(i);
   }
   return long(r+0.5); // round to nearest
}

// calculate (1-sT)(x) where s is some complex scale factor (typically norm 1)
// The identity and rho matrices are used to remove any spurious component
// in the direction of the identity
struct OneMinusTransfer
{
   OneMinusTransfer(std::complex<double> ScaleFactor, LinearWavefunction const& Psi, QuantumNumber const& QShift, 
                    MatrixOperator const& Rho,
                    MatrixOperator const& Identity)
      : Scale_(ScaleFactor), Psi_(Psi), QShift_(QShift), Rho_(Rho), Identity_(Identity) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-Scale_*delta_shift(transfer_from_left(x, Psi_), QShift_);
      //      r = delta_shift(r, QShift_);
      if (is_scalar(r.TransformsAs()))
          r -= inner_prod(r, Rho_) * Identity_; // orthogonalize to the identity
      return r;
   }

   std::complex<double> Scale_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
   MatrixOperator const& Rho_;
   MatrixOperator const& Identity_;
};

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator Rhs)
{
   MatrixOperator Guess = Rhs;
   int m = 30;
   int max_iter = 10000;
   double tol = 1e-15;
   DEBUG_TRACE("running GMRES");
   GmRes(Guess, F, Rhs, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());
   return Guess;
}

template <typename T>
struct EigenPair
{
   T LeftVector;
   T RightVector;
   std::complex<double> Eigenvalue;
};

void remove_redundant(OperatorComponent& Op);

int main(int argc, char** argv)
{
   std::string FName;
   std::string Operator;

   double t = 1;
   double tc = 0.25;
   double U = 0;
   int Power = 1;
   bool Verbose = false;

   std::cout.precision(getenv_or_default("MP_PRECISION", 14));
 
   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("t", prog_opt::value(&t),
	  FormatDefault("nearest-neighbor hopping (for hubbard etc)", t).c_str())
	 ("tc", prog_opt::value(&tc),
	  FormatDefault("cluster hopping (for triangular cluster)", tc).c_str())
	 ("U", prog_opt::value(&U),
	  FormatDefault("coulomb repulsion", U).c_str())
	 ("power", prog_opt::value(&Power),
	  FormatDefault("Calculate expectation value of operator to this power", Power).c_str())
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "extra debug output")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&FName), "wavefunction")
         ("operator", prog_opt::value<std::string>(&Operator), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("operator", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("operator") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-expectation3 <psi1> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }
      
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      InfiniteWavefunction Psi = *PsiPtr;
      int UnitCellSize = Psi.Psi.size();

      TriangularMPO Op;

      
      if (Operator == "su2spin")
      {
	 double Spin = 0.5;
	 LatticeSite Site = CreateSU2SpinSite(Spin);	 
	 Op = TriangularOneSite(Site["S"]);
	 Op = -sqrt(3.0) * prod(Op, Op, QuantumNumber(Op.GetSymmetryList()));
      }
      else
      {
	 std::cerr << "Unknown operator!\n";
	 exit(1);
      }

      // construct the operator to the given power
      TriangularMPO Temp = Op;
      while (Power > 1)
      {
	 Op = Op * Temp;
	 --Power;
      }

      // Make a LinearWavefunction in the symmetric orthogonality constraint
      MatrixOperator Identity = MatrixOperator::make_identity(Psi.Psi.Basis1());
      MatrixOperator Rho = scalar_prod(Psi.C_right, herm(Psi.C_right));
      LinearWavefunction Phi = Psi.Psi; // no need to bugger around with C_old,C_right
 
      MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
      MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);

      Phi.set_front(prod(LambdaInvSqrt, Phi.get_front()));
      Phi.set_back(prod(Phi.get_back(), delta_shift(LambdaSqrt, adjoint(Psi.QShift))));
      Rho = Psi.C_old;
      Identity = Rho;

      // Rho and Identity are the same matrix in this representation

      //   TRACE(norm_frob(operator_prod(herm(Phi), Identity, Phi) - Identity));
      //   TRACE(norm_frob(operator_prod(Phi, Rho, herm(Phi)) - Rho));

      // make Op the same size as our unit cell
      Op = repeat(Op, UnitCellSize / Op.size());

      std::map<std::complex<double>, Polynomial<MatrixOperator>, CompareComplex> 
	 E = SolveMPO_Left(Phi, Psi.QShift, Op, Rho, Identity, Verbose);
      Polynomial<std::complex<double> > aNorm = ExtractOverlap(E[1.0], Rho);

      std::cout << "#degree #real #imag\n";
      for (int i = 0; i <= aNorm.degree(); ++i)
	 {
	    std::cout << i << ' ' << aNorm[i].real() << ' ' << aNorm[i].imag() << '\n';
	 }
      std::cout << std::endl;

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
