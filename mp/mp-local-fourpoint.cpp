// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-local-fourpoint.cpp
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
#include <boost/optional.hpp>
#include <boost/none.hpp>

int main(int argc, char** argv)
{
   if (argc < 9 || argc > 10)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-local-fourpoint <lattice> <psi> <Op1> <Op2> <Op3> <Op4> "
         "<first> <last> [qn]\n";
      std::cerr << "Calculates the four point correlator Op1(i-1) Op2(i) Op3(j) Op3(j+1)\n";
      std::cerr << "[qn] is the quantum number to project the coupling <Op3> <Op4>.\n";
      std::cerr << "This is only needed for non-abelian symmetries.\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);

   std::string Op1 = argv[3];
   std::string Op2 = argv[4];
   std::string Op3 = argv[5];
   std::string Op4 = argv[6];
   int FirstSite = boost::lexical_cast<int>(argv[7]);
   int LastSite = boost::lexical_cast<int>(argv[8]);
   boost::optional<QuantumNumbers::QuantumNumber> MiddleProjection = boost::none;
   if (argc >= 10)
      MiddleProjection = QuantumNumbers::QuantumNumber(Psi1->GetSymmetryList(), argv[9]);

   Lattice Lat = System->GetLattice();

   CenterWavefunction Psi = *Psi1;
   Psi1 = pvalue_ptr<MPWavefunction>();

   for (int i = 1; i < FirstSite; ++i)
      Psi.RotateRight();

   std::cout.precision(12);

   // This holds the E matrices for the left system
   typedef std::map<int, MatrixOperator> OpMapType;
   MatrixOperator SLast;
   OpMapType OpMap;
   //   int iLast = 0;
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
      SiteBlock::const_iterator I2 = Lat[i].find(Op2);
      if (I2 != Lat[i].end() && !SLast.is_null())
      {
         SimpleOperator MyOp2 = I2->second;
         if (MiddleProjection)
            OpMap[i] = operator_prod(herm(MyOp2), herm(Psi.Left()), SLast, Psi.Left(), 
                                     adjoint(*MiddleProjection));
         else
            OpMap[i] = operator_prod(herm(MyOp2), herm(Psi.Left()), SLast, Psi.Left());
      }
      SiteBlock::const_iterator I1 = Lat[i].find(Op1);
      if (I1 != Lat[i].end())
      {
         SimpleOperator MyOp1 = I1->second;
         MatrixOperator LeftIdentity = MatrixOperator::make_identity(Psi.Left().Basis1());
         SLast = operator_prod(herm(MyOp1), herm(Psi.Left()), LeftIdentity, Psi.Left());
      }

      // For the right operator, construct the F matrix and the expectation value
      SiteBlock::const_iterator I3 = Lat[i+1].find(Op3);
      if (I3 != Lat[i+1].end())
      {
         MPStateComponent PsiR2 = Psi.LookupRight(Psi.RightSize()-2);
         SiteBlock::const_iterator J = Lat[i+2].find(Op4);
         if (J == Lat[i+2].end())
            break;

         SimpleOperator MyOp4 = J->second;
         MatrixOperator Ident = MatrixOperator::make_identity(PsiR2.Basis2());
         MatrixOperator F2 = operator_prod(MyOp4, PsiR2, Ident, herm(PsiR2));

         SimpleOperator MyOp3 = I3->second;
         MatrixOperator F;
         if (MiddleProjection)
            F = operator_prod(MyOp3, Psi.Right(), F2, herm(Psi.Right()), *MiddleProjection);
         else
            F = operator_prod(MyOp3, Psi.Right(), F2, herm(Psi.Right()));

         for (OpMapType::iterator mI = OpMap.begin(); mI != OpMap.end(); ++mI)
         {
            std::complex<double> Res = inner_prod(Psi.Center(), 
                                                  triple_prod(mI->second, 
                                                              Psi.Center(), 
                                                              herm(F)));
            std::cout << std::setw(5) << Lat.coordinate_at_site(mI->first) << "   " 
                      << std::setw(5) << Lat.coordinate_at_site(i+1) << "   "
                      << std::setw(18) << Res.real() << "   " 
                      << std::setw(18) << Res.imag() << '\n';
         }
      }

      if (i != LastSite-1)
         Psi.RotateRight();
   }

   pheap::Shutdown();
}
