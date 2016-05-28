// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-debug.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>

void ShowWavefunc(CenterWavefunction Psi, std::ostream& out, int Max = -1)
{
   out << "Symmetry list is " << Psi.GetSymmetryList() << '\n';
   out << "State transforms as " << Psi.TransformsAs() << '\n';
   out << "Number of sites = " << Psi.size() << '\n';
   while (Psi.RightSize() > 1)
   {
      out << "\nMatrices at partition (" << Psi.LeftSize() << ", " << Psi.RightSize() << ")\n";
      out << "=========================================\n";

      out << "Left basis: " << Psi.Left().Basis2() 
          << "Right basis: " << Psi.Right().Basis1() << '\n';
      out << "Center matrix:\n" << Psi.Center() << '\n';
      out << "|Center| = " << norm_frob(Psi.Center()) << '\n';

      out << "\nLeft matrices:\n" << Psi.Left();
      out << "\n\nRight matrices:\n" << Psi.Right();

      out << "\n\nLeft norm matrix:\n"
          << scalar_prod(herm(Psi.Left()), Psi.Left()) << '\n';

      out << "\nRight norm matrix:\n"
          << scalar_prod(Psi.Right(), herm(Psi.Right())) << '\n';

      DensityMatrix<MatrixOperator> DM(scalar_prod(Psi.Center(), herm(Psi.Center())));
      out << "\nReduced density matrix:\n";
      DM.DensityMatrixReport(out, Max);
      Psi.RotateRight();
   }

   DensityMatrix<MatrixOperator> DM(scalar_prod(Psi.Center(), herm(Psi.Center())));
   out << "\nReduced density matrix at partition (" 
       << Psi.LeftSize() << "," << Psi.RightSize() << ") :\n";
   DM.DensityMatrixReport(out, Max);
}


int main(int argc, char** argv)
{
   if (argc < 2 || argc > 3)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-info <psi> [max-eigenvalues]\n";
      return 1;
   }

   int Max = -1;
   if (argc == 3) Max = boost::lexical_cast<int>(argv[2]);

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], 655360, true);

   std::cout.precision(13);
   ShowWavefunc(CenterWavefunction(*Psi), std::cout, Max);

   pheap::Shutdown();
}
