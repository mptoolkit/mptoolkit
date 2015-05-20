// -*- C++ -*- $Id$

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
#include "lattice/triangular-parser.h"
#include "mps/momentum_operations.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "common/prog_opt_accum.h"

namespace prog_opt = boost::program_options;

void PrintFormat(std::complex<double> const& Value, bool ShowRealPart, bool ShowImagPart, 
		 bool ShowMagnitude, bool ShowArgument,
		 bool ShowRadians)
{
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "    ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "    ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "    ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
	 Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
}

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
   int Verbose = 0;
   int NMax = 3;
   int NLegs = 1;
   double Spin = 0.5;
   bool Print = false;
   bool ShowRealPart = false;
   bool ShowImagPart = false;
   bool ShowMagnitude = false;
   bool ShowArgument = false;
   bool ShowRadians = false;
   bool ShowPolar = false;
   bool ShowCartesian = false;
   bool Quiet = false;
   bool Columns = false;

   std::cout.precision(getenv_or_default("MP_PRECISION", 14));
   std::cerr.precision(getenv_or_default("MP_PRECISION", 14));
 
   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("power", prog_opt::value(&Power),
	  FormatDefault("Calculate expectation value of operator to this power", Power).c_str())
	 ("cart,c", prog_opt::bool_switch(&ShowCartesian),
	  "show the result in cartesian coordinates [equivalent to --real --imag]")
	 ("polar,p", prog_opt::bool_switch(&ShowPolar),
	  "show the result in polar coodinates [equivalent to --mag --arg]")
	 ("real,r", prog_opt::bool_switch(&ShowRealPart),
	  "display the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImagPart),
	  "display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
	 ("radians", prog_opt::bool_switch(&ShowRadians),
	  "display the argument in radians instead of degrees")
	 ("columns", prog_opt::bool_switch(&Columns),
	  "Show the prefactors of each degree in columns rather than rows")
	 ("quiet,q", prog_opt::bool_switch(&Quiet), "Don't show column headings")
         ("print,p", prog_opt::bool_switch(&Print), "Print the MPO to standard output")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "extra debug output (can be used more than once)")
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
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " <psi1> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }
      
      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
	  && !ShowCartesian && !ShowPolar && !ShowArgument
	  && !ShowRadians)
      {
	 ShowCartesian = true;
	 ShowPolar = true;
      }
      
      if (ShowCartesian)
      {
	 ShowRealPart = true;
	 ShowImagPart = true;
      }
      if (ShowPolar)
      {
	 ShowMagnitude = true;
	 ShowArgument = true;
      }
      if (ShowRadians)
	 ShowArgument = true;      

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

      if (Print)
      {
	 print_structure(Op, std::cout);
	 //	 std::cout << Op << '\n';
	 //std::cout << "\nTransfer matrix:" << construct_transfer_matrix(herm(GenericMPO(Op)),
	 //								GenericMPO(Op)) << '\n';
      };


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
      if (WavefuncUnitCellSize % Op.size() != 0)
      {
         std::cout << "mp-iexpectation: fatal: the wavefunction unit cell "
            "must be a multiple of the operator unit cell.\n";
         return 1;
      }

      Op = repeat(Op, WavefuncUnitCellSize / Op.size());

      KMatrixPolyType E = SolveMPO_Left(Phi, Psi.QShift, Op, Identity, Rho, Verbose);
      Polynomial<std::complex<double> > aNorm = ExtractOverlap(E[1.0], Rho);

      if (!Quiet)
      {
	 std::cout << "#" << argv[0];
	 for (int i = 1; i < argc; ++i)
	    std::cout << ' ' << argv[i];
	 std::cout << "\n#quantities are calculated per unit cell size of " << WavefuncUnitCellSize 
		   << (WavefuncUnitCellSize == 1 ? " site\n" : " sites\n");
	 if (Columns)
	 {
	    int d = aNorm.degree();
	    for (int i = 1; i <= d; ++i)
	    {
	       if (ShowRealPart)
		  std::cout << '#' << i << "-real                 ";
	       if (ShowImagPart)
		  std::cout << '#' << i << "-imag                 ";
	       if (ShowMagnitude)
		  std::cout << '#' << i << "-magnitude            ";
	       if (ShowArgument)
		  std::cout << '#' << i << "-argument" << (ShowRadians ? "(rad)" : "(deg)") << "        ";
	    }
	 }
	 else
	 {
	    std::cout << "#degree     ";
	    if (ShowRealPart)
	       std::cout << "#real                   ";
	    if (ShowImagPart)
	       std::cout << "#imag                   ";
	    if (ShowMagnitude)
	       std::cout << "#magnitude              ";
	    if (ShowArgument)
	       std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
	 }
         std::cout << '\n';
         std::cout << std::left;
      }

      for (int i = 1; i <= aNorm.degree(); ++i)
      {
	 if (!Columns)
	    std::cout << std::setw(11) << i << ' ';
	 PrintFormat(aNorm[i], ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
	 if (!Columns)
	    std::cout << std::endl;
      }
      if (Columns)
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
