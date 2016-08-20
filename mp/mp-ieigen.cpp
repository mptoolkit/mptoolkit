// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ieigen.cpp
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

#include "matrixproduct/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "matrixproduct/triangularoperator.h"
#include "models/spin.h"
#include "matrixproduct/infinitewavefunction.h"

// this needs to be generalized to non-trivial unit cells

std::complex<double> SolveLeftEquations(MatrixOperator& Guess, MatrixOperator const& y,
                                        LinearWavefunction const& Psi)
{
   // simplest version, simple iteration
   MatrixOperator Ident = MatrixOperator::make_identity(Psi.Basis1());
   double N = trace(Ident).real();
   int Iterations = 100;
   std::complex<double> Offset;
   for (int i = 0; i < Iterations; ++i)
   {
      Guess = pass_from_left(Guess, Psi) + y;
      // remove the trace component
      Offset = trace(Guess) / N;
      Guess -= Offset * Ident;
   }
   return Offset;
}

MPStateComponent SolveRight(LinearWavefunction const& Psi, SimpleIrredMPOperator const& Op)
{
   CHECK_EQUAL(Op.size(), Psi.size());
   std::vector<std::list<MatrixOperator> > R(Op.size());
   // initial operator is always I
   int i = 0;
   for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
      R[i++].push_back(MatrixOperator::make_identity(I->Basis2()));
   }

   while (R.back().size() < Psi.Basis2().size())
   {
      bool OpAddedThisPass = false; // set to true if we add a new element to R this pass.
      i = R.size();
      LinearWavefunction::const_iterator I = Psi.end();
      while (I != Psi.begin())
      {
         --I;
         --i;   // the input
         int r = i-1; // the output
         if (r < 0) r = R.size()-1;

         // set as many terms of R[r] that we can from the elements in R[i].
         // There are two cases here; either we can set the element exactly
         // as Op[i] has non-zero elements only in the first R[i].size()
         // columns of the R[r].size() row.  Or, the element
         // Op[i](R[r].size(), R[i].size()) is the identity operator.
         // In this case, we must end up with a complete set at all sites
         // in the unit cell that satisfy this property, and we have a
         // nonhomogeneous equation to solve.
         // alternatively, we could handle only the first case here
         // and defer the nonhomogeneous equation to a second pass.
         bool done = false;
         while (!done)
         {
            int n = R[r].size();
            for (int j = R[i].size(); j < Op[i].Basis2().size(); ++j)
            {
               if (iterate_at(Op[i], n, j))
                  done = true;
            }
            if (done) continue;

            OpAddedThisPass = true;
            MatrixOperator m;
            for (int j = 0; j < R[i].size(); ++j)
            {
               if (iterate_at(Op[i], n, j))
               {
                  m += operator_prod(Op[i](n,j), *I, R[i][j], herm(*I));
               }
            }
            R[r].push_back(m);
         }

      }

      if (!OpAddedThisPass)
      {
         // if we get here, we must have a nonhomogeneous equation to solve
      }
   }
}

MPStateComponent SolveLeft(LinearWavefunction const& Psi, SimpleMPOperator const& Op, std::complex<double>& LastEigen)
{
   CHECK_EQUAL(Op.size(), Psi.size());
   int const Size = Op.front().Basis().size();
   std::vector<std::deque<MatrixOperator> > R(Op.size());
   // initial operator is always I
   int i = 0;
   for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
      R[i++].push_front(MatrixOperator::make_identity(I->Basis1()));
   }
   //   std::complex<double> LastEigen = 0.0;
   for (int n = 2; n < Size; ++n)
   {
      // Add terms to R
      i = 0;
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
      {
         int CurrentRank = R[i].size();
         int Cols = Op[i].Basis2().size();
         // while column
         R[i++].push_front(MatrixOperator::make_identity(I->Basis1()));
      }


      LinearAlgebra::Vector<SimpleOperator> c = project_row(Op.data(), Size-n);

      // get the nonhomogeneous term
      MatrixOperator m;
      for (int j = 1; j < n; ++j)
      {
         m += operator_prod(herm(c[Size-j]), herm(Psi.get_front()), R[j-1], Psi.get_front());
      }

      if (c[Size-n].is_null())
      {
         // no element on the diagonal.  This is the simple case
         R.push_front(m);
      }
      else
      {
         // an element on the diagonal.  Assume it is an identity operator
         // and solve the resulting nonhomogeneous equations
         MatrixOperator r = m;// MakeRandomMatrixOperator(Psi.Basis1(), Psi.Basis1(), Op.Basis()[Size-n]);
         LastEigen = SolveLeftEquations(r, m, Psi);
         R.push_front(r);
      }
   }

   MPStateComponent Result(Op.LocalBasis(), Psi.Basis1(), Psi.Basis1());
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = R[i];
   }
   return Result;
}

std::complex<double>
MPO_Eigenvalues(MPStateComponent& Guess, LinearWavefunction const& Psi, MpOpTriangular const& Op)
{
   std::complex<double> Eigen = 0.0;
   Guess = SolveLeft(Psi, Op, Eigen);
   return Eigen;
}

int main(int argc, char** argv)
{

   if (argc < 2 || argc > 4)
   {
      std::cout << "usage: mp-ieigen <wavefunction> [<psi2>] [rotate]\n";
      return 1;
   }

   std::string FName = argv[1];
   std::string FName2 = argc >= 3 ? argv[2] : FName;
   int Rotate = argc >= 4 ? boost::lexical_cast<int>(argv[3]) : 0;

   pvalue_ptr<InfiniteWavefunction> Psi1 = pheap::OpenPersistent(FName, mp_pheap::CacheSize(), true);
   pvalue_ptr<InfiniteWavefunction> Psi2 = pheap::ImportHeap(FName2);

   std::cout.precision(14);

   double Lambda = 1;
   SiteBlock Site = CreateSpinSite(0.5);

#if 0
   MpOpTriangular Ham = 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
      + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
#endif
#if 0
   MpOpTriangular Ham = TriangularOneSite(Site["Sx"]) + (0.9993747067577124 * 0.5 * TriangularOneSite(Site["I"]));
#endif

   // SU(2) Heisenberg model
#if 1
   MpOpTriangular Ham = TriangularTwoSite(-sqrt(3.0)*Site["S"],
                                          Site["S"],
                                          Site["I"].TransformsAs());
#endif

   //   Ham = Ham*Ham;

   //TRACE(Ham);

   LinearWavefunction p1 = get_orthogonal_wavefunction(*Psi1);

   if (Rotate > 0)
      *Psi2.mutate() = rotate_left(*Psi2, Rotate);

   LinearWavefunction p2 = get_orthogonal_wavefunction(*Psi2);

   MPStateComponent Guess = Initial_E(Ham, p.Basis1());

   std::complex<double> EVal = MPO_Eigenvalues(Guess, p, Ham);

   TRACE(EVal);

   pheap::Shutdown();
}
