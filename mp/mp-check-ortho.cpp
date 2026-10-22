// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-check-ortho.cpp
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "wavefunction/mpwavefunction.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace prog_opt = boost::program_options;

struct Residual
{
   double Abs = 0.0;
   double Rel = 0.0;
};

Residual
residual(MatrixOperator const& X, MatrixOperator const& Ref)
{
   double const RefNorm = norm_frob(Ref);
   double const Abs = norm_frob(X - Ref);
   if (!std::isfinite(RefNorm) || !std::isfinite(Abs))
      throw std::runtime_error("non-finite orthogonality residual");
   double const Rel = RefNorm != 0.0 ? Abs / RefNorm
      : (Abs == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
   return {Abs, Rel};
}

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      int Verbose = 0;
      bool Quiet = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("quiet,q", prog_opt::bool_switch(&Quiet), "only print machine-readable values")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "show per-site residuals; use twice for extra structure checks")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&PsiStr), "psi")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi>\n";
         std::cerr << "Report canonical-form residuals of a stored InfiniteWavefunctionLeft without modifying it.\n"
                   << "Columns are Frobenius norms: abs = ||actual-reference||, rel = abs/||reference||.\n"
                   << "max_left_orthogonality: maximum site residual against the identity.\n"
                   << "max_stored_density_bond: maximum mismatch with stored lambda^2 while\n"
                   << "  propagating the right boundary density left through the unit cell.\n"
                   << "right_density_fixed_point: full-cell density residual, including qshift.\n"
                   << "For a zero reference, rel is zero if abs is zero, otherwise infinity.\n"
                   << "Residuals exclude the separate log amplitude; no pass/fail tolerance is applied.\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);
      if (!PsiPtr->is<InfiniteWavefunctionLeft>())
      {
         std::cerr << "mp-check-ortho: fatal: only InfiniteWavefunctionLeft is currently supported.\n";
         pheap::Shutdown();
         return 1;
      }

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();
      if (Psi.empty())
         throw std::runtime_error("empty wavefunction");
      if (Verbose > 1)
         Psi.check_structure();

      // This exposes the stored tensors and lambda, without reorthogonalizing.
      LinearWavefunction PsiLinear;
      RealDiagonalOperator Lambda;
      std::tie(PsiLinear, Lambda) = get_left_canonical(Psi);

      double MaxLocalAbs = 0.0;
      double MaxLocalRel = 0.0;
      int Site = 0;
      for (auto I = PsiLinear.begin(); I != PsiLinear.end(); ++I, ++Site)
      {
         MatrixOperator Input = MatrixOperator::make_identity(I->Basis1());
         MatrixOperator Output = operator_prod(herm(*I), Input, *I);
         MatrixOperator Ident = MatrixOperator::make_identity(I->Basis2());
         Residual R = residual(Output, Ident);
         MaxLocalAbs = std::max(MaxLocalAbs, R.Abs);
         MaxLocalRel = std::max(MaxLocalRel, R.Rel);
         if (Verbose > 0)
         {
            std::cout << "site " << Site
                      << " left_orthogonality_abs " << R.Abs
                      << " left_orthogonality_rel " << R.Rel << '\n';
         }
      }

      double MaxBondDensityAbs = 0.0;
      double MaxBondDensityRel = 0.0;
      MatrixOperator Rho = Lambda * Lambda;
      MatrixOperator BondRho = Rho;
      int Bond = Psi.size();
      auto I = PsiLinear.end();
      while (I != PsiLinear.begin())
      {
         --I;
         --Bond;
         BondRho = operator_prod(*I, BondRho, herm(*I));
         MatrixOperator StoredRho = Psi.lambda(Bond) * Psi.lambda(Bond);
         Residual R = residual(BondRho, StoredRho);
         MaxBondDensityAbs = std::max(MaxBondDensityAbs, R.Abs);
         MaxBondDensityRel = std::max(MaxBondDensityRel, R.Rel);
         if (Verbose > 0)
         {
            std::cout << "bond " << Bond
                      << " right_density_abs " << R.Abs
                      << " right_density_rel " << R.Rel << '\n';
         }
      }

      MatrixOperator RhoOut = delta_shift(BondRho, adjoint(Psi.qshift()));
      Residual Right = residual(RhoOut, Rho);

      if (!Quiet)
      {
         std::cout << "#quantity #abs #rel\n";
      }
      std::cout << "max_left_orthogonality " << MaxLocalAbs << ' ' << MaxLocalRel << '\n';
      std::cout << "max_stored_density_bond " << MaxBondDensityAbs << ' ' << MaxBondDensityRel << '\n';
      std::cout << "right_density_fixed_point " << Right.Abs << ' ' << Right.Rel << '\n';

      pheap::Shutdown();
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
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
