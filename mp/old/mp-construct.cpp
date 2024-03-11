// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-construct.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp/copyright.h"
#include <fstream>
#include <sstream>
#include "common/environment.h"

// CreateRank1Wavefunction
// Creates a rank-1 state given a list of <QuantumNumber, std::string> pairs.
// The quantum number refers to the height, the string is the label of the
// basis state used to achieve that height.
template <typename FwdIter>
CenterWavefunction CreateRank1Wavefunction(Lattice const& L, FwdIter First, FwdIter Last);

template <typename FwdIter>
CenterWavefunction CreateRank1Wavefunction(Lattice const& L, FwdIter First, FwdIter Last)
{
   CHECK_EQUAL(std::distance(First, Last), L.size());

   BasisList Vac = make_vacuum_basis(L.GetSymmetryList());
   VectorBasis B2(Vac);
   QuantumNumber Ident(L.GetSymmetryList());  // the scalar quantum number

   CenterWavefunction Result;

   MatrixOperator Center(B2, B2, Ident);
   Center(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   int i = 1;
   while (First != Last)
   {
      VectorBasis B1(L.GetSymmetryList());
      B1.push_back(First->first, 1);
      MatrixOperator Next(B1, B2, L[i].Basis1().qn(First->second));
      Next(0,0) = LinearAlgebra::Matrix<double>(1,1,1);
      MPStateComponent R(L[i].Basis1().Basis(), B1, B2);
      R[L[i].Basis1().Lookup(First->second)] = Next;
      Result.PushRight(prod(R, Center));
      Center = TruncateBasis1(Result.Right());
      B2 = B1;
      ++First;
      ++i;
   }

   Result.PushLeft(prod(Center, Result.Right()));
   Result.PopRight();
   Result.Center() = TruncateBasis2(Result.Left());

   Result.normalize();
   return Result;
}

bool SiteBasisHasLabel(SiteBasis const& Basis, std::string const& s)
{
   for (std::size_t i = 0; i < Basis.size(); ++i)
   {
      if (Basis.Label(i) == s) return true;
   }
   return false;
}

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-construct <lattice> <states> <outfile>\n";
      return 1;
   }

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(argv[3], 1, PageSize, CacheSize);
   pvalue_ptr<OperatorList> OpList = pheap::ImportHeap(argv[1]);

   SymmetryList SList = OpList->GetSymmetryList();
   Lattice L = OpList->GetLattice();

   std::ifstream In(argv[2]);
   std::list<std::pair<QuantumNumber, std::string> > Heights;
   QuantumNumber Q = QuantumNumber(SList);
   int LineNum = 0;
   for (int i = L.size(); i >= 1; --i)
   {
      ++LineNum;
      SiteBasis Basis = L[i].Basis1();
      std::string Str;
      std::getline(In, Str);
      std::istringstream Line(Str);
      Line.seekg(0);
      std::string State, q;
      Line >> State;
      if (!SiteBasisHasLabel(Basis, State))
      {
         // currently we allow this case only when the site basis is one-dimensional.
         // we could support other cases and disambiguate on quantum number alone,
         // but that is a bit messy.
         if (Basis.size() != 1)
         {
            std::cerr << "Error in " << argv[2] << " at line " << LineNum
                      << ": '" << State
                      << "' is not a valid site basis element, possible choices are ";
            for (std::size_t s = 0; s < Basis.size(); ++s) std::cerr << Basis.Label(s) << ' ';
            std::cerr << '\n';
            exit(1);
         }
         q = State;
         State = Basis.Label(0);
      }
      else
      {
         Line >> q;
      }

      QuantumNumbers::QuantumNumberList PossibleQ
         = transform_targets(Q, Basis[Basis.Lookup(State)].second);
      if (q == "")
      {
         if (PossibleQ.size() != 1)
         {
            std::cerr << "Error in " << argv[2] << " at line " << LineNum
                      << ": a quantum number must be specified\n"
                      << "   to disambiguate between ";
            std::copy(PossibleQ.begin(), PossibleQ.end(),
                      std::ostream_iterator<QuantumNumber>(std::cerr, " "));
            std::cerr << '\n';
            exit(1);
         }
         Q = PossibleQ[0];
      }
      else
      {
         bool Found = false;
         for (std::size_t n = 0; n < PossibleQ.size() && !Found; ++n)
         {
            if (PossibleQ[n].ToString() == q)
            {
               Q = PossibleQ[n];
               Found = true;
            }
         }
         if (!Found)
         {
            std::cerr << "Error in " << argv[2] << " at line " << LineNum
                      << ": quantum number " << q << " is not possible; possible are ";
            std::copy(PossibleQ.begin(), PossibleQ.end(),
                      std::ostream_iterator<QuantumNumber>(std::cerr, " "));
            std::cerr << '\n';
            exit(1);
         }
      }
      Heights.push_back(std::make_pair(Q, State));
   }
   char c;
   In >> std::skipws >> c;
   if (In)
   {
      std::cerr << "warning: extraneous input in " << argv[2]
                << " beyond line " << LineNum << '\n';
   }

   pvalue_ptr<MPWavefunction> Psi
      = new MPWavefunction(CreateRank1Wavefunction(OpList->GetLattice(),
                                                   Heights.begin(),
                                                   Heights.end()).AsLinearWavefunction());
   pheap::ShutdownPersistent(Psi);
}
