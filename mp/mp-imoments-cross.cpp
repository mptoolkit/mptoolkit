// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-imoments.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

//
// Consistency of mp-ioverlap and mp-imoments.
//
// Given a Hamiltonian H, we can calculate the Loschmidt echo
// r(z) = (1/N) ln <psi|exp[zH]|psi> = (1/N) ln <psi_0 | psi_t>
// This is simply the log of the largest magnitude transfer matrix eigenvalue
// between psi_0 and psi_t.
// Note that we didn't put a minus sign in the definition of r(z), although
// this is commonly done.
//
// The derivatives of r(z) with respect to z are
// d^n r(z) / d z^n = n^th cumulant c_n of <psi_0 | H^n | psi_t>
//
// Given these cumulants we can expand r(z) as
// r(z_0 + z) = r(z_0) + \sum_n c_n z^n / n!
//
// where we define z = -(beta + it) is the complex time

#include "mpo/basic_triangular_mpo.h"
#include "mp-algorithms/transfer.h"
#include "wavefunction/mpwavefunction.h"
#include "wavefunction/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "common/polynomial.h"
#include "tensor/tensor_eigen.h"
#include "mp/copyright.h"
#include "common/prog_options.h"
#include "lattice/infinite-parser.h"
#include "wavefunction/momentum_operations.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "common/prog_opt_accum.h"
#include <boost/algorithm/string.hpp>
#include "common/openmp.h"

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

inline
int ipow(int x, int p)
{
   DEBUG_CHECK(p >= 0);
   if (p == 0)
      return 1;
   else if (p == 1)
      return p;
   else if (p % 2 == 0)
      return ipow(x*x, p/2);
   else
      return x*ipow(x*x, (p-1)/2);
}

void ShowCumulants(std::vector<std::complex<double> > const& Cumulants,
                   bool OneLine, bool Quiet, bool ShowRealPart, bool ShowImagPart,
                   bool ShowMagnitude, bool ShowArgument, bool ShowRadians,
                   double ScaleFactor)
{
   std::cout << std::left;
   if (!Quiet)
   {
      if (OneLine)
      {
         for (int i = 1; i < Cumulants.size(); ++i)
         {
               std::string Suffix = "_" + std::to_string(i);
            if (ShowRealPart)
               std::cout << std::setw(24) << std::string("#real" + Suffix);
            if (ShowImagPart)
               std::cout << std::setw(24) << std::string("#imag" + Suffix);
            if (ShowMagnitude)
               std::cout << std::setw(24) << std::string("#magnitude" + Suffix);
            if (ShowArgument)
               std::cout << std::setw(24) << std::string("#argument" + Suffix + (ShowRadians ? "(rad)" : "(deg)"));
            }
      }
      else
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
   }
   for (unsigned n = 1; n < Cumulants.size(); ++n)
   {
      if (!OneLine)
         std::cout << std::setw(9) << n << ' ';
      PrintFormat(Cumulants[n]*ScaleFactor, ShowRealPart,
                  ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
      if (!OneLine)
         std::cout << '\n';
   }
   if (OneLine)
      std::cout << '\n';
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
      // This is the easy case, we can calculate the n'th cumulant from the n'th moment
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

void ShowAllComponents(std::vector<KMatrixPolyType> const& E, MatrixOperator Rho)
{
   std::cout << "Showing all components.  Total columns = " << E.size() << "\n";
   for (unsigned i = 0; i < E.size(); ++i)
   {
      for (auto const& x : E[i])
      {
	 std::cout << "Column " << (i+1) << ": momentum: " << x.first << ": ";
	 Polynomial<std::complex<double>> m = ExtractOverlap(x.second, Rho);
	 std::cout << m;
      }
   }
}

int main(int argc, char** argv)
{
   omp::initialize();
   std::string Psi1Str;
   std::string OpStr;
   std::string Psi2Str;

   int Power = 1;
   int Verbose = 0;
   int UnitCellSize = 0;
   int Degree = 0;
   int Rotate = 0;
   bool Reflect = false;
   bool Conj = false;
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
   bool ShouldShowAllComponents = false;
   std::string Sector;
   bool OneLine = false;

   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("quantumnumber,q", prog_opt::value(&Sector), "quantum number sector of the transfer matrix")
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
         ("quiet", prog_opt::bool_switch(&Quiet), "don't show column headings")
         ("oneline", prog_opt::bool_switch(&OneLine), "Show all output on one line (currently only works with --cumulants)")
         ("rotate", prog_opt::value(&Rotate),
          "rotate the unit cell of psi1 this many sites to the left before calculating the overlap [default 0]")
         ("reflect", prog_opt::bool_switch(&Reflect),
          "reflect psi1 (gives parity eigenvalue)")
         ("conj", prog_opt::bool_switch(&Conj),
          "complex conjugate psi1")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Linear solver convergence tolerance", Tol).c_str())
         ("unityepsilon", prog_opt::value(&UnityEpsilon),
          FormatDefault("Epsilon value for testing eigenvalues near unity", UnityEpsilon).c_str())
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "extra debug output (can be used more than once)")
	      ("showall", prog_opt::bool_switch(&ShouldShowAllComponents), "show all columns of the fixed-point "
          "solutions (mostly for debugging)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value<std::string>(&Psi1Str), "wavefunction1")
         ("operator", prog_opt::value<std::string>(&OpStr), "operator")
         ("psi2", prog_opt::value<std::string>(&Psi2Str), "wavefunction2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi1", 1);
      p.add("operator", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi1> <operator> <psi2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

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
      pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::OpenPersistent(Psi2Str, CacheSize, true);
      pvalue_ptr<MPWavefunction> Psi1Ptr = pheap::ImportHeap(Psi1Str);
      InfiniteWavefunctionLeft Psi2 = Psi2Ptr->get<InfiniteWavefunctionLeft>();
      InfiniteWavefunctionLeft Psi1 = Psi1Ptr->get<InfiniteWavefunctionLeft>();
      auto q = QuantumNumber(Psi1.GetSymmetryList(), Sector);


      // Rotate as necessary.  Do this BEFORE determining the quantum number sectors!
      if (Verbose)
      {
         std::cout << "Rotating Psi1 right by" << Rotate << " sites\n";
      }
      Psi1.rotate_right(Rotate);
      if (Reflect)
      {
         if (Verbose)
            std::cout << "Reflecting psi1..." << std::endl;
         inplace_reflect(Psi1);
      }
      if (Conj)
      {
         if (Verbose)
            std::cout << "Conjugating psi1..." << std::endl;
         inplace_conj(Psi1);
      }

      // We require that the wafefunctions have the same quantum number per unit cell.  If this isn't the case,
      // the transfer matrix eigenmatrices have a weird form that would require changing 'q' (sector) of the operator.
      // In principle this is possible, but is beyond what delta_shift currently does.
      if (Psi1.qshift() != Psi2.qshift())
      {
         std::cerr << "fatal: wavefunctions have a different qshift, which is not supported.\n";
         return 1;
      }

      // Since the operator is a positional argument we need to include it.  But allow it to be empty, and
      // that will take it to be the Psi2 Hamiltonian.  If that is empty, look at the psi1 Hamiltonian.
      if (OpStr.empty())
      {
         OpStr = Psi2Ptr->Attributes()["Hamiltonian"].get_or_default(std::string());
         if (OpStr.empty())
         {
            OpStr = Psi1Ptr->Attributes()["Hamiltonian"].get_or_default(std::string());
         }
         else
         {
            if (Verbose > 1)
               std::cerr << "Taking operator from psi2 Hamiltonian attribute.\n";
         }
         if (OpStr.empty())
         {
            std::cerr <<  basename(argv[0]) << ": fatal: no operator specified, and wavefunction "
               "attribute Hamiltonian does not exist or is empty.\n";
            return 1;
         }
         else
         {
            if (Verbose > 1)
               std::cerr << "Taking operator from psi1 Hamiltonian attribute.\n";
         }
      }

      BasicTriangularMPO Op;

      InfiniteLattice Lattice;
      std::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

      int Size = statistics::lcm(Psi1.size(), Psi2.size(), Op.size());

      // TODO: a better approach would be to get SolveMPO to understand how to do
      // a multiple unit-cell operator.  But actually that is difficult.
      Psi1 = repeat(Psi1, Size / Psi1.size());
      Psi2 = repeat(Psi2, Size / Psi2.size());
      Op = repeat(Op, Size / Op.size());

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
         UnitCellSize = Psi1.size();
      double ScaleFactor = double(UnitCellSize) / double(Size);

      if (!Quiet)
      {
         print_preamble(std::cout, argc, argv);
         std::cout << "#operator " << EscapeArgument(OpStr) << '\n';
         std::cout << "#quantities are calculated per unit cell size of " << UnitCellSize
                   << (UnitCellSize == 1 ? " site\n" : " sites\n");
         std::cout << "#quantum number sector is " << q << '\n';
      }

      std::complex<double> lambda;
      MatrixOperator TLeft, TRight;
      std::tie(lambda, TLeft, TRight) = get_transfer_eigenpair(Psi1, Psi2, q);

      if (!Quiet)
      {
         std::cout << "#transfer matrix eigenvalue = " << formatting::format_complex(lambda) << '\n';
      }

      TRight = delta_shift(TRight, Psi2.qshift());

      // Check that the local basis for the wavefunction and hamiltonian are compatible
      if (ExtractLocalBasis(Psi1) != ExtractLocalBasis1(Op))
      {
         std::cerr << "fatal: operator is defined on a different local basis to the wavefunction.\n";
         return 1;
      }

      if (ExtractLocalBasis1(Op) != ExtractLocalBasis2(Op))
      {
         std::cerr << "fatal: operator has different domain and co-domain.\n";
         return 1;
      }

      BasicTriangularMPO OriginalOp = Op;  // keep a copy so we can do repeated powers

      std::vector<Polynomial<std::complex<double> > > Moments;

      LinearWavefunction Phi1 = get_left_canonical(Psi1).first;
      LinearWavefunction Phi2 = get_left_canonical(Psi2).first;

      // This does the equivalent.  Probably would have been easier to implement it this way in the first place,
      // but maybe more general to have lambda as an explicit argument?
      Phi2 *= (1.0 / lambda);
      lambda = 1.0;

      // first power
      std::vector<KMatrixPolyType> E;
      SolveMPO_Left_Cross(E, Phi1, Phi2, Psi1.qshift(), Op, TLeft, TRight, lambda, Power > 1, Degree, Tol, UnityEpsilon, Verbose);
      Moments.push_back(ExtractOverlap(E.back()[1.0], TRight));
      if (ShouldShowAllComponents)
      {
         ShowAllComponents(E, TRight);
      }
      // If we're not calculating the cumulants, then we can print the moments as we calculate them.
      // BUT, if we have Verbose > 0, then don't print anything until the end, so that it doesn't get
      // mixed up with the verbose output.
      if (CalculateMoments && !Quiet && Verbose <= 0)
         ShowMomentsHeading(ShowRealPart, ShowImagPart,
                            ShowMagnitude, ShowArgument, ShowRadians);

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

         // 2016-08-15: Due to a misunderstanding/bug this reuse doesnt work.
         E = std::vector<KMatrixPolyType>();

         Op = OriginalOp * Op;
         SolveMPO_Left_Cross(E, Phi1, Phi2, Psi1.qshift(), Op, TLeft, TRight, lambda, p < Power-1, Degree*(p+1),
                       Tol, UnityEpsilon, Verbose);
         Moments.push_back(ExtractOverlap(E.back()[1.0], TRight));
         // Force the degree of the MPO
         if (Degree != 0)
            Moments.back()[Degree*(p+1)] += 0.0;
	 //            Moments.back()[ipow(Degree, p+1)] += 0.0;
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
         ShowCumulants(Cumulants, OneLine, Quiet, ShowRealPart,
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
