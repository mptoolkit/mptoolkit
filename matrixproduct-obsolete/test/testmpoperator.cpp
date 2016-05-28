// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/test/testmpoperator.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "models/hubbard-so4.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/lattice.h"
#include "common/trace.h"
#if 0
template <typename T>
struct PrintRunLengthCompressed : public boost::static_visitor<>
{
   PrintRunLengthCompressed(std::ostream& out_, int NestLevel_) : out(out_), NestLevel(NestLevel_) {}

   void DoNest() const
   {
      out << std::string(NestLevel*4, ' ');
   }

   void operator()(T const& x) const
   {
      this->DoNest();
      out << "Value:\n";
      this->DoNest();
      out << x << '\n';
   };

   void operator()(run_length_array<T> const& x) const
   {
      this->DoNest();
      out << "Array: size " << x.size() << '\n';
      int i = 0;
      for (typename run_length_array<T>::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         this->DoNest();
         out << "Array element " << i++ << '\n';
         I->apply_visitor(PrintRunLengthCompressed<T>(out, NestLevel+1));
      }
   }

   void operator()(run_length_repeat<T> const& x) const
   {
      this->DoNest();
      out << "Repeat: size " << x.size() << '\n';
      x.nested().apply_visitor(PrintRunLengthCompressed<T>(out, NestLevel+1));
   }

   std::ostream& out;
   int NestLevel;
};

template <typename T>
std::ostream& operator<<(std::ostream& out, run_length_compressed<T> const& x)
{
   x.apply_visitor(PrintRunLengthCompressed<T>(out, 0));
   return out;
}

std::ostream& operator<<(std::ostream& out, MPOperator const& x)
{
   return out << x.data();
}
#endif
int main()
{
   SymmetryList Symmetry("Q:SU(2),S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2,QuantumNumbers::SU2> QN(Symmetry);

   SiteBlock A = CreateSO4HubbardSiteA();
   SiteBlock B = CreateSO4HubbardSiteB();

   Lattice L = join(Lattice(A), Lattice(B));
   L = repeat(L, 6);

   MPOperator C6 = CreateMPOperator(L, "C", 6);
   MPOperator C8 = CreateMPOperator(L, "C", 8);

   TRACE(std::complex<double>(3,2)*(C6*C8)*(C8*C6));
}
