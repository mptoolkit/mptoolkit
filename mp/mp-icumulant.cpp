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

void ShowMoments(std::vector<Polynomial<std::complex<double> > > const& Moments, 
		 bool Quiet, bool Columns, bool ShowRealPart, bool ShowImagPart, 
		 bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
{
   if (!Quiet)
   {
      std::cout << "#moment #degree ";
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

   for (unsigned n = 0; n < Moments.size(); ++n)
   {
      int m = Moments[n].degree();
      for (int i = 1; i <= m; ++i)
      {
	 std::cout << std::setw(7) << m << ' ' << std::setw(7) << i << ' ';
	 PrintFormat(Moments[n][i], ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
	 std::cout << '\n';
      }
   }
}

void ShowCumulants(std::vector<std::complex<double> > const& Cumulants, 
		   bool Quiet, bool Columns, bool ShowRealPart, bool ShowImagPart, 
		   bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
{
   if (!Quiet)
   {
      std::cout << "#cumulant ";
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
   for (unsigned n = 1; n < Cumulants.size(); ++n)
   {
      std::cout << std::setw(9) << n << ' ';
      PrintFormat(Cumulants[n], ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
      std::cout << '\n';
   }
}

// given an array of moment polynomials (ordered by degree, and in multiples of the
// degree of the first moment), calculate the corresponding cumulants
std::vector<std::complex<double> >
MomentsToCumulants(std::vector<Polynomial<std::complex<double> > > const& Moments)
{
   int const Degree = Moments.back().degree();
   int const FirstMoment = Moments.front().degree();
   std::vector<std::complex<double> > Cumulants(Degree+1, 0.0);
   
   if (FirstMoment == 1)
   {
      // This is the easy case, we can calculate the n'th cumulant from the n'th cumulant
      // The complication that we handle is that possibly kappa_1 is zero but
      // kappa_1^2 is non-zero.
      CHECK_EQUAL(int(Moments.size()), Degree);
      for (int n = 0; n < Degree; ++n)
      {
	 Cumulants[n+1] = Moments[n][1];
      }
      // special case for kappa_1
      if (Degree >= 2)
      {
	 std::complex<double> k1 = std::sqrt(Moments[1][2]);
	 if (norm_2_sq(k1) > 10.0*norm_2_sq(Cumulants[1]))
	 {
	    Cumulants[1] = k1;
	 }
      }
   }
   else if (FirstMoment == 2)
   {
      // we have only every second moment
      // mu_2 = kappa_2 L + kappa_1^2 L^2
      std::cout << "#WARNING: sign of kappa_1 is unspecified.\n";
      Cumulants[1] = std::sqrt(Moments[0][2]); 
      Cumulants[2] = Moments[0][1];

      if (Moments.size() > 1)
      {
	 // next two cumulants
	 // mu_4 = kappa_4 L + (4 \kappa_3 \kappa_1 + 3 \kappa_2^2) L^2 
	 //        + 6 \kappa_2 \kappa_1^2 L^3 + \kappa_1^4 L^4
	 // NOTE: we cannot get the sign of kappa_3 correct in this case
	 std::cout << "#WARNING: sign of kappa_3 is unspecified.\n";
	 Cumulants[3] = (Moments[1][2] - 3.0*Cumulants[2]*Cumulants[2]) / (4.0 * Cumulants[1]);
	 Cumulants[4] = Moments[1][1];

	 if (Moments.size() > 3)
	 {
	    // next two cumulants
	    // mu_6 = kappa_6 L 
	    //        + (6 kappa_5 kappa_1 + 15 kappa_2 kappa_4 + 10 kappa_3^2) L^2
	    //        + (15 kappa_4 kappa_1^2 + 60 kappa_3 kappa_2 kappa_1 + 15 kappa_2^3) L^3
	    //        + 45 kappa_2^2 kappa_1^2 L^4
	    //        + 15 kappa_2 kappa_1^4 L^5
	    //        + kappa_1^6 L^6  
	    std::cout << "#WARNING: sign of kappa_5 is unspecified.\n";
	    Cumulants[5] = (Moments[2][2] - 15.0 * Cumulants[2]*Cumulants[4] 
			    - 10.0*Cumulants[3]*Cumulants[3]) / (6.0 * Cumulants[1]);
	    Cumulants[6] = Moments[2][1];
	    
	    if (Moments.size() > 4)
	    {
	       std::cout << "#WARNING: Cumulants > 6 are not yet implemented!\n";
	    }
	 }
      }
   }
   else
   {
      PANIC("First moment is higher than degree 2, not yet implemented!");
   }
   return Cumulants;
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
   bool CalculateCumulants = false;
   double UnityEpsilon = DefaultEigenUnityEpsilon;

   std::cout.precision(getenv_or_default("MP_PRECISION", 14));
   std::cerr.precision(getenv_or_default("MP_PRECISION", 14));
 
   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("power", prog_opt::value(&Power),
	  FormatDefault("Calculate expectation value of operator to this power", Power).c_str())
         ("cumulants,u", prog_opt::bool_switch(&CalculateCumulants),
          "calculate the commulants kappa")
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
	 ("unityepsilon", prog_opt::value(&UnityEpsilon),
	  FormatDefault("Epsilon value for testing eigenvalues near unity", UnityEpsilon).c_str())
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

      if (!Quiet)
      {
	 std::cout << "#" << argv[0];
	 for (int i = 1; i < argc; ++i)
	    std::cout << ' ' << argv[i];
	 std::cout << std::endl;
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

      if (!Quiet)
      {
	 std::cout << "#quantities are calculated per unit cell size of " << WavefuncUnitCellSize 
		   << (WavefuncUnitCellSize == 1 ? " site\n" : " sites\n");
      }
      
      TriangularMPO Op;

      InfiniteLattice Lattice;
      boost::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

      if (Print)
      {
	 print_structure(Op, std::cout, UnityEpsilon);
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

      TriangularMPO OriginalOp = Op;  // keep a copy so we can do repeated powers

      std::vector<Polynomial<std::complex<double> > > Moments;

      // first power
      std::vector<KMatrixPolyType> 
	 E = SolveMPO_Left(Phi, Psi.QShift, Op, Identity, Rho, Power > 1, UnityEpsilon, Verbose);
      Moments.push_back(ExtractOverlap(E.back()[1.0], Rho));

      // loop over the powers of the operator
      for (int p = 1; p < Power; ++p)
      {
	 // construct the operator to the given power
	 // The way we have defined prod() for MPO's, is A*B is
	 // [A_00 B   A_01 B ... ]
	 // [A_10 B   A_11 B ... ]
	 // [ ..       ..        ]
	 // That is, in order to re-use the E matrices for a higher power, we need
	 // to multiply the new operator on the left, where the identity in the top right corner (A_00)
	 // will correspond to the already calculated terms.
	 Op = OriginalOp * Op;
	 E = SolveMPO_Left(Phi, Psi.QShift, E, Op, Identity, Rho, p < Power-1, UnityEpsilon, Verbose);
	 Moments.push_back(ExtractOverlap(E.back()[1.0], Rho));
      }

      if (CalculateCumulants)
      {
	 std::vector<std::complex<double> > Cumulants = MomentsToCumulants(Moments);

	 ShowCumulants(Cumulants, Quiet, Columns, ShowRealPart, 
		       ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
      }
      else
	 ShowMoments(Moments, Quiet, Columns, ShowRealPart, ShowImagPart, 
		     ShowMagnitude, ShowArgument, ShowRadians);


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
