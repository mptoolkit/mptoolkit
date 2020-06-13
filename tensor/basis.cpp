// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/basis.cpp
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

#include "basis.h"
#include "common/statistics.h"

namespace Tensor
{

//
// BasisList
//

bool
BasisList::is_identity() const
{
   return this->size() == 1 && is_scalar(this->front());
}

bool
BasisList::is_regular() const
{
   std::set<QuantumNumber> Used;
   for (const_iterator I = this->begin(); I != this->end(); ++I)
   {
      if (Used.find(*I) != Used.end())
         return false;

      Used.insert(*I);
   }
   return true;
}

int BasisList::total_degree() const
{
   int Deg = 0;
   for (std::size_t i = 0; i < this->size(); ++i)
   {
      Deg += degree(Q_[i]);
   }
   return Deg;
}

int
BasisList::find_first(QuantumNumber const& q) const
{
   unsigned n = 0;
   while (n < Q_.size() && Q_[n] != q)
      ++n;

   return (n < Q_.size()) ? int(n) : -1;
}

int
BasisList::find_next(QuantumNumber const& q, int n) const
{
   ++n;
   while (n < int(Q_.size()) && Q_[n] != q)
      ++n;

   return (n < int(Q_.size())) ? n : -1;
}

std::ostream& operator<<(std::ostream& Out, BasisList const& B)
{
   Out << "Basis has symmetry " << B.GetSymmetryList()
       << ", size = " << B.size() << ", degree = " << B.total_degree() << '\n';
   Out << "   N  QuantumNumber\n";
   int N = 0;
   std::map<QuantumNumbers::QuantumNumber, int> Dim;
   for (auto I = B.begin(); I != B.end(); ++I)
   {
      Out << std::setw(4) << N << "  "
          << std::setw(13) << *I << '\n';
      ++Dim[*I];
      ++N;
   }
   Out << "\nDimension of each quantum number:\n";
   Out << "QuantumNumber    Dimension\n";
   for (auto x : Dim)
   {
      Out << std::left << std::setw(16) << x.first << ' ' << x.second << '\n';
   }
   return Out;
}

std::ostream&
operator<<(std::ostream& out, std::vector<BasisList> const& u)
{
   out << "vector of BasisList, length " << u.size() << '\n';
   for (auto const& x : u)
   {
      out << x << '\n';
   }
   return out;
}

void
BasisList::delta_shift(QuantumNumbers::QuantumNumber const& q)
{
   Q_.delta_shift(q);
}

BasisList delta_shift(BasisList const& Orig, QuantumNumbers::QuantumNumber const& q)
{
   BasisList Result(Orig.GetSymmetryList());
   for (BasisList::const_iterator I = Orig.begin(); I != Orig.end(); ++I)
      Result.push_back(delta_shift(*I, q));

   return Result;
}


#if 0
BasisList RenameSymmetry(BasisList const& BL, SymmetryList const& NewSL)
{
   BasisList Result(NewSL);
   for (unsigned i = 0; i < BL.size(); ++i)
      Result.push_back(RenameSymmetry(BL[i], NewSL));
   return Result;
}
#endif

std::map<int, int>
LinearizeQuantumNumberSubspace(BasisList const& b, QuantumNumbers::QuantumNumber const& q)
{
   std::map<int, int> Result;
   int n = 0;
   for (unsigned i = 0; i < b.size(); ++i)
   {
      if (b[i] == q)
         Result[i] = n++;
   }
   return Result;
}

//
// VectorBasis
//

int VectorBasis::total_degree() const
{
   int Deg = 0;
   for (std::size_t i = 0; i < this->size(); ++i)
   {
      Deg += degree(Basis_[i]) * Dimension_[i];
   }
   return Deg;
}

std::ostream& operator<<(std::ostream& Out, VectorBasis const& B)
{
   Out << "Basis has symmetry " << B.GetSymmetryList()
       << ", subspace size = " << B.size()
       << ", dimension = " << B.total_dimension()
       << ", degree = " << B.total_degree() << '\n';
   Out << "   N  QuantumNumber  Dimension\n";
   for (std::size_t N = 0; N < B.size(); ++N)
   {
      Out << std::setw(4) << N << "  "
          << std::setw(13) << B[N] << "  "
          << std::setw(9) << B.dim(N) << '\n';
   }
   return Out;
}

void
VectorBasis::delta_shift(QuantumNumbers::QuantumNumber const& q)
{
   Basis_.delta_shift(q);
}

VectorBasis delta_shift(VectorBasis const& Orig, QuantumNumbers::QuantumNumber const& q)
{
   return VectorBasis(delta_shift(Orig.Basis(), q),
                      Orig.Dimension_.begin(),
                      Orig.Dimension_.end());
}

VectorBasis
CoerceSymmetryList(VectorBasis const& b, SymmetryList const& sl)
{
   VectorBasis Result;
   Result.Basis_ = CoerceSymmetryList(b.Basis_, sl);
   Result.Dimension_ = b.Dimension_;
   return Result;
}

#if 0
VectorBasis RenameSymmetry(VectorBasis const& BL, SymmetryList const& NewSL)
{
   return VectorBasis(RenameSymmetry(BL.Basis_, NewSL), BL.Dimension_);
}
#endif

bool
is_compatible(std::vector<BasisList> const& a, std::vector<BasisList> const& b)
{
   int Size = statistics::lcm(a.size(), b.size());
   for (int i = 0; i < Size; ++i)
   {
      if (a[i%a.size()] != b[i%b.size()])
	 return false;
   }
   return true;
}

} // namespace Tensor
