// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-chiral-correlator-su2.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/proccontrol.h"

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-chiral-correlator <lattice> <psi> <first> <last>\n";
      std::cerr << "Calculates the chiral correlator S(i-1)xS(i) . S(j)xS(j+1)!\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);

   std::string Op1 = "S";
   int FirstSite = boost::lexical_cast<int>(argv[3]);
   std::string Op2 = "S";
   int LastSite = boost::lexical_cast<int>(argv[4]);

   Lattice Lat = System->GetLattice();

   CenterWavefunction Psi = *Psi1;
   Psi1 = pvalue_ptr<MPWavefunction>();

   for (int i = 1; i < FirstSite; ++i)
      Psi.RotateRight();

   std::cout.precision(12);

   // The s=1 quantum number
   QuantumNumber qVector(Lat.GetSymmetryList(), "1");

   // This holds the E matrices for the left system
   typedef std::map<int, MatrixOperator> OpMapType;
   MatrixOperator SLast;
   OpMapType OpMap;
   for (int i = FirstSite; i < LastSite-1; ++i)
   {
      // Update all the existing E matrices
      for (OpMapType::iterator mI = OpMap.begin(); mI != OpMap.end(); ++mI)
      {
         mI->second = operator_prod(herm(Psi.Left()),
                                    mI->second,
                                    Psi.Left());
      }

      // Add another E matrix for the new site, only if the left operator exists in the block
      SiteBlock::const_iterator I = Lat[i].find(Op1);
      if (I != Lat[i].end())
      {
         SimpleOperator MyOp = I->second;

         if (!SLast.is_null())
         {
            OpMap[i] = operator_prod(herm(MyOp), herm(Psi.Left()), SLast, Psi.Left(), qVector);
         }
         MatrixOperator LeftIdentity = MatrixOperator::make_identity(Psi.Left().Basis1());
         SLast = operator_prod(herm(MyOp), herm(Psi.Left()), LeftIdentity, Psi.Left());
      }

      // For the right operator, construct the F matrix and the expectation value
      I = Lat[i+1].find(Op2);
      if (I != Lat[i+1].end())
      {
         MPStateComponent PsiR2 = Psi.LookupRight(Psi.RightSize()-2);
         SiteBlock::const_iterator J = Lat[i+2].find(Op2);
         if (J == Lat[i+2].end())
            break;

         SimpleOperator MyS = J->second;
         MatrixOperator Ident = MatrixOperator::make_identity(PsiR2.Basis2());
         MatrixOperator F2 = operator_prod(MyS, PsiR2, Ident, herm(PsiR2));

         SimpleOperator MyOp = I->second;
         MatrixOperator F = operator_prod(MyOp, Psi.Right(), F2, herm(Psi.Right()), qVector);
         for (OpMapType::iterator mI = OpMap.begin(); mI != OpMap.end(); ++mI)
         {
            std::complex<double> Res = inner_prod(Psi.Center(),
                                                  triple_prod(mI->second,
                                                              Psi.Center(),
                                                              herm(F)));
            // The 2.0 here for the sqrt(2)'s that arise in the cross product of S=1
            std::cout << std::setw(5) << Lat.coordinate_at_site(mI->first) << "   "
                      << std::setw(5) << Lat.coordinate_at_site(i+1) << "   "
                      << std::setw(18) << sqrt(12.0) * Res.real() << "   "
                      << std::setw(18) << sqrt(12.0) * Res.imag() << '\n';
         }
      }

      if (i != LastSite-1)
         Psi.RotateRight();
   }

   pheap::Shutdown();
}
