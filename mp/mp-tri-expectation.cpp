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

#include "mpo/basic_triangular_mpo.h"
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
#include "common/randutil.h"

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
   std::string FName;
   std::string OpStr;
   int N = 100;
   bool ShowAll = false;
   bool Randomize = false;

   int Verbose = 0;

   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
	 ("sites,n", prog_opt::value(&N), "maximum length to calculate (sites)")
	 ("showall", prog_opt::bool_switch(&ShowAll), "show all columns of the MPO")
	 ("randomize", prog_opt::bool_switch(&Randomize), "randomize boundary tensors")
         ("help", "show this help message")
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

      if (vm.count("help") > 0 || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " <psi1> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();
      int WavefuncUnitCellSize = Psi.size();

      if (!vm.count("operator"))
      {
         OpStr = PsiPtr->Attributes()["Hamiltonian"].get_or_default(std::string());
         if (OpStr.empty())
         {
            std::cerr <<  basename(argv[0]) << ": fatal: no operator specified, and wavefunction "
               "attribute Hamiltonian does not exist or is empty.\n";
            return 1;
         }
      }

      BasicTriangularMPO Op;
      InfiniteLattice Lattice;
      std::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

      RealDiagonalOperator D;
      LinearWavefunction Phi;
      std::tie(Phi, D) = get_left_canonical(Psi);

      MatrixOperator Rho = D;
      Rho = scalar_prod(Rho, herm(Rho));
      Rho = delta_shift(Rho, Psi.qshift());

      MatrixOperator Identity = MatrixOperator::make_identity(Phi.Basis1());

      // Check that the local basis for the wavefunction and hamiltonian are compatible
      if (ExtractLocalBasis(Psi) != ExtractLocalBasis1(Op))
      {
         std::cerr << "fatal: operator is defined on a different local basis to the wavefunction.\n";
         return 1;
      }

      if (ExtractLocalBasis1(Op) != ExtractLocalBasis2(Op))
      {
         std::cerr << "fatal: operator has different domain and co-domain.\n";
         return 1;
      }

      StateComponent E = Initial_E(Op, Psi.Basis1());
      if (Randomize)
      {
	 srand(randutil::crypto_rand());
	 for (unsigned j = 1; j < E.size(); ++j)
	 {
	    E[j] = MakeRandomMatrixOperator(E[j].Basis1(), E[j].Basis2(), E[j].TransformsAs());
	 }
      }
      LinearWavefunction::const_iterator PsiI = Phi.cbegin();
      for (unsigned i = 0; i < N; ++i)
      {
	 E = contract_from_left(Op[i%Op.size()], herm(*PsiI), E, *PsiI);
	 ++PsiI;
	 if (PsiI == Phi.end())
	    PsiI = Phi.begin();
	 if (ShowAll)
	 {
	    for (unsigned j = 0; j < E.size(); ++j)
	    {
	       if (E[j].TransformsAs() == Rho.TransformsAs())
	       {
		  std::cout << (i+1) << " column " << (j+1) << ' ' << format_complex(inner_prod(E[j], Rho)) << '\n';
	       }
	    }
	 }
	 else
	 {
	    if (E.back().TransformsAs() == Rho.TransformsAs())
	    {
	       std::cout << (i+1) << ' ' << format_complex(inner_prod(E.back(), Rho)) << '\n';
	    }
	 }
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
