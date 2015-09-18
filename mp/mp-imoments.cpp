// -*- C++ -*-

#include "mpo/triangular_mpo.h"
#include "mps/infinitewavefunction.h"
#include "mps/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"
#include "common/polynomial.h"
#include "tensor/tensor_eigen.h"
#include "mp/copyright.h"
#include "common/prog_options.h"
#include "lattice/infinite-parser.h"
#include "mps/momentum_operations.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "common/prog_opt_accum.h"
#include <boost/algorithm/string.hpp>

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

void ShowMomentsHeading(bool ShowRealPart, bool ShowImagPart, 
			bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
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
   std::cout << std::endl;
}

void ShowMoments(Polynomial<std::complex<double> > const& Moments, 
		 bool ShowRealPart, bool ShowImagPart, 
		 bool ShowMagnitude, bool ShowArgument, bool ShowRadians,
		 double ScaleFactor)
{
   std::cout << std::left;
   int m = Moments.degree();
   for (int i = 1; i <= m; ++i)
   {
      std::cout << std::setw(7) << m << ' ' << std::setw(7) << i << ' ';
      PrintFormat(Moments[i] * pow(ScaleFactor, double(i)), 
		  ShowRealPart, ShowImagPart, ShowMagnitude, 
		  ShowArgument, ShowRadians);
      std::cout << std::endl;
   }
}

void ShowCumulants(std::vector<std::complex<double> > const& Cumulants,
		   bool Quiet, bool ShowRealPart, bool ShowImagPart, 
		   bool ShowMagnitude, bool ShowArgument, bool ShowRadians,
		   double ScaleFactor)
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
      std::cout << '\n';
   }
   std::cout << std::left;
   for (unsigned n = 1; n < Cumulants.size(); ++n)
   {
      std::cout << std::setw(9) << n << ' ';
      PrintFormat(Cumulants[n]*ScaleFactor, ShowRealPart, 
		  ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
      std::cout << '\n';
   }
}

// given an array of moment polynomials (ordered by degree, and in multiples of the
// degree of the first moment), calculate the corresponding cumulants
std::vector<std::complex<double> >
MomentsToCumulants(std::vector<Polynomial<std::complex<double> > > const& Moments,
		   double Epsilon = 1E-15, bool Quiet = false)
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
      (Quiet ? std::cerr : std::cout)
	 << "#WARNING: sign of kappa_1 is unspecified, choosing +ve value.\n";
      Cumulants[1] = std::sqrt(Moments[0][2]); 
      Cumulants[2] = Moments[0][1];

      if (Moments.size() > 1)
      {
	 // next two cumulants
	 // mu_4 = kappa_4 L + (4 \kappa_3 \kappa_1 + 3 \kappa_2^2) L^2 
	 //        + 6 \kappa_2 \kappa_1^2 L^3 + \kappa_1^4 L^4
	 // NOTE: we cannot get the sign of kappa_3 correct in this case
	 (Quiet ? std::cerr : std::cout)
	    << "#WARNING: sign of kappa_3 is relative to the sign of kappa_1.\n";
	 // if kappa_1 is very small then we can have a catastrophic loss of precision here.
	 // The subtraction mu_1(2) - 3*kappa_2^2 may be very small.
	 std::complex<double> Diff = Moments[1][2] - 3.0*Cumulants[2]*Cumulants[2];
	 double Eps = std::abs(Diff) / (std::abs(Moments[1][2]) 
					+ std::abs(3.0*Cumulants[2]*Cumulants[2]));
	 if (Eps < Epsilon*1000)
	 {
	    (Quiet ? std::cerr : std::cout)
	       << "# *** WARNING *** catastrophic loss of precision in kappa_3 *** \n";
	 }
	 DEBUG_TRACE(Diff)(Eps)(Epsilon);
	 Cumulants[3] = (Diff) / (4.0 * Cumulants[1]);
	 Cumulants[4] = Moments[1][1];

	 if (Moments.size() > 2)
	 {
	    // next two cumulants
	    // mu_6 = kappa_6 L 
	    //        + (6 kappa_5 kappa_1 + 15 kappa_2 kappa_4 + 10 kappa_3^2) L^2
	    //        + (15 kappa_4 kappa_1^2 + 60 kappa_3 kappa_2 kappa_1 + 15 kappa_2^3) L^3
	    //        + 45 kappa_2^2 kappa_1^2 L^4
	    //        + 15 kappa_2 kappa_1^4 L^5
	    //        + kappa_1^6 L^6  
	    (Quiet ? std::cerr : std::cout)
	       << "#WARNING: sign of kappa_5 is relative to the sign of kappa_1.\n";
	    Cumulants[5] = (Moments[2][2] - 15.0 * Cumulants[2]*Cumulants[4] 
			    - 10.0*Cumulants[3]*Cumulants[3]) / (6.0 * Cumulants[1]);
	    Cumulants[6] = Moments[2][1];
	    
	    if (Moments.size() > 3)
	    {
	       (Quiet ? std::cerr : std::cout)
		  << "#WARNING: Cumulants > 6 are not yet implemented!\n";
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
   int UnitCellSize = 0;
   int Degree = 0;
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
   bool CalculateMoments = false;
   bool CalculateCumulants = false;
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   double Tol = 1E-15;

   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("power", prog_opt::value(&Power),
	  FormatDefault("Calculate expectation value of operator to this power", Power).c_str())
         ("moments", prog_opt::bool_switch(&CalculateMoments),
          "calculate the moments [default, unless --cumulants is specified]")
         ("cumulants,t", prog_opt::bool_switch(&CalculateCumulants),
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
	 ("unitcell,u", prog_opt::value(&UnitCellSize),
	  "scale the results to use this unit cell size [default wavefunction unit cell]")
	 ("degree,d", prog_opt::value(&Degree),
	  "force setting the degree of the MPO")
	 ("quiet,q", prog_opt::bool_switch(&Quiet), "Don't show column headings")
         ("print,p", prog_opt::bool_switch(&Print), "Print the MPO to standard output")
	 ("tol", prog_opt::value(&Tol),
	  FormatDefault("Linear solver convergence tolerance", Tol).c_str())
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
         std::cerr << "usage: " << basename(argv[0]) << " <psi1> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));
 
      if (!Quiet)
	 print_preamble(std::cout, argc, argv);
      
      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
	  && !ShowCartesian && !ShowPolar && !ShowArgument)
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

      // If none of --cumulants and --moments are specified, then the default is to calculate moments
      if (!CalculateMoments && !CalculateCumulants)
	 CalculateMoments = true;

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      InfiniteWavefunction Psi = *PsiPtr;
      int WavefuncUnitCellSize = Psi.Psi.size();

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
	 UnitCellSize = WavefuncUnitCellSize;
      double ScaleFactor = double(UnitCellSize) / double(WavefuncUnitCellSize);

      if (!Quiet)
      {
	 std::cout << "#quantities are calculated per unit cell size of " << UnitCellSize 
		   << (UnitCellSize == 1 ? " site\n" : " sites\n");
      }
      
      TriangularMPO Op;

      InfiniteLattice Lattice;
      boost::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

      if (Print)
      {
	 print_structure(Op, std::cout, UnityEpsilon);
	 if (Verbose > 2)
	    std::cout << Op << '\n';
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
         std::cout << "mp-icumulant: fatal: the wavefunction unit cell "
            "must be a multiple of the operator unit cell.\n";
         return 1;
      }

      Op = repeat(Op, WavefuncUnitCellSize / Op.size());

      // Check that the local basis for the wavefunction and hamiltonian are compatible
      local_basis_compatible_or_abort(Psi.Psi, Op);
      
      TriangularMPO OriginalOp = Op;  // keep a copy so we can do repeated powers

      std::vector<Polynomial<std::complex<double> > > Moments;

      // If we're not calculating the cumulants, then we can print the moments as we calculate them.
      // BUT, if we have Verbose > 0, then don't print anything until the end, so that it doesn't get
      // mixed up with the verbose output.
      if (CalculateMoments && !Quiet && Verbose <= 0)
	 ShowMomentsHeading(ShowRealPart, ShowImagPart, 
			    ShowMagnitude, ShowArgument, ShowRadians);

      // first power
      std::vector<KMatrixPolyType> E;
      SolveMPO_Left(E, Phi, Psi.QShift, Op, Identity, Rho, Power > 1, Tol, UnityEpsilon, Verbose);
      Moments.push_back(ExtractOverlap(E.back()[1.0], Rho));
      // Force the degree of the MPO
      if (Degree != 0)
	 Moments.back()[Degree] += 0.0;
      if (CalculateMoments && Verbose <= 0)
	 ShowMoments(Moments.back(), ShowRealPart, ShowImagPart, 
		     ShowMagnitude, ShowArgument, ShowRadians,
		     ScaleFactor);

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
	 SolveMPO_Left(E, Phi, Psi.QShift, Op, Identity, Rho, p < Power-1, 
		       Tol, UnityEpsilon, Verbose);
	 Moments.push_back(ExtractOverlap(E.back()[1.0], Rho));
	 // Force the degree of the MPO
	 if (Degree != 0)
	    Moments.back()[std::pow(Degree, p+1)] += 0.0;
	 if (CalculateMoments && Verbose <= 0)
	    ShowMoments(Moments.back(), ShowRealPart, ShowImagPart, 
			ShowMagnitude, ShowArgument, ShowRadians,
			ScaleFactor);
      }

      // if we had verbose output, then we delay printing the moments until now
      if (CalculateMoments && Verbose > 0)
      {
	 if (!Quiet)
	    ShowMomentsHeading(ShowRealPart, ShowImagPart, 
			       ShowMagnitude, ShowArgument, ShowRadians);
	 for (unsigned i = 0; i < Moments.size(); ++i)
	 {
	    ShowMoments(Moments[i], ShowRealPart, ShowImagPart, 
			ShowMagnitude, ShowArgument, ShowRadians,
			ScaleFactor);
	 }
      }

      if (CalculateCumulants)
      {
	 // If we calculated the moments, then put in a blank line to separate them
	 if (CalculateMoments)
	    std::cout << '\n';

	 std::vector<std::complex<double> > Cumulants = MomentsToCumulants(Moments, Tol, Quiet);
	 ShowCumulants(Cumulants, Quiet, ShowRealPart, 
		       ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians,
		       ScaleFactor);
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
