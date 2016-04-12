
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
#include <fstream>
#include "lattice/infinitelattice-parser.h"
#include "mps/momentum_operations.h"
#include "mp-algorithms/triangular_mpo_solver.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   std::string FName;
   std::string OpStr;

   double t = 0;
   double t2 = 0;
   double tc = 0;
   double U = 0;
   double Lambda = 0;
   int Power = 1;
   bool Verbose = false;
   int NMax = 3;
   int NLegs = 1;
   double Spin = 0.5;
   std::string CouplingFile;

   std::cout.precision(getenv_or_default("MP_PRECISION", 14));
 
   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("power", prog_opt::value(&Power),
	  FormatDefault("Calculate expectation value of operator to this power", Power).c_str())
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "extra debug output")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&FName), "wavefunction")
         ("operator", prog_opt::value<std::string>(&OpStr), "operator")
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
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-expectation3 <psi1> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }
      
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      InfiniteWavefunction Psi = *PsiPtr;
      int WavefuncUnitCellSize = Psi.Psi.size();

      TriangularMPO Op;

      InfiniteLattice Lattice;
      boost::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

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
 
      Rho = delta_shift(Rho, Psi.QShift);

#if 0
      MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
      MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);

      Phi.set_front(prod(LambdaInvSqrt, Phi.get_front()));
      Phi.set_back(prod(Phi.get_back(), delta_shift(LambdaSqrt, adjoint(Psi.QShift))));
      Rho = Psi.C_old;
      Identity = Rho;
#endif

      // Rho and Identity are the same matrix in this representation

      //   TRACE(norm_frob(operator_prod(herm(Phi), Identity, Phi) - Identity));
      //   TRACE(norm_frob(operator_prod(Phi, Rho, herm(Phi)) - Rho));

      // make Op the same size as our unit cell
      if (UnitCellSize % Op.size() != 0)
      {
         std::cout << "mp-iexpectation: fatal: the wavefunction unit cell "
            "must be a multiple of the operator unit cell.\n";
         return 1;
      }

      Op = repeat(Op, WavefuncUnitCellSize / Op.size());

      KMatrixPolyType E = SolveMPO_Left(Phi, Psi.QShift, Op, Identity, Rho, Verbose);
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
