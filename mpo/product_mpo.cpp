// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/product_mpo.cpp
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

#include "product_mpo.h"
#include "common/statistics.h"

ProductMPO::ProductMPO(GenericMPO const& Other)
   : Data(Other)
{
   CHECK_EQUAL(Data.Basis1(), Data.Basis2());
}

std::ostream&
operator<<(std::ostream& out, ProductMPO const& x)
{
   return out << x.data();
}

void print_structure(ProductMPO const& Op, std::ostream& out, double UnityEpsilon)
{
   out << "ProductMPO has " << Op.size() << " sites\n";
   for (int i = 0; i < Op.size(); ++i)
   {
      out << "Site " << i << " dimension " << Op[i].size1() << " x " << Op[i].size2() << '\n';
      print_structure(Op[i], out, UnityEpsilon);
   }
}

bool
ProductMPO::is_identity() const
{
   OperatorClassification c = classify(this->data());
   return c.is_identity();
}

bool
ProductMPO::is_string() const
{
   return this->Basis1().size() == 1 && this->Basis2().size() == 1;
}

ProductMPO
ProductMPO::make_identity(std::vector<BasisList> const& Basis)
{
   // delegate to the BasicFiniteMPO version
   return ProductMPO(BasicFiniteMPO::make_identity(Basis));
}

ProductMPO
ProductMPO::make_identity(std::vector<BasisList> const& Basis, QuantumNumber const& q)
{
   // delegate to the BasicFiniteMPO version
   return ProductMPO(BasicFiniteMPO::make_identity(Basis, q));
}

PStream::opstream&
operator<<(PStream::opstream& out, ProductMPO const& op)
{
   out << op.data();
   return out;
}

PStream::ipstream&
operator>>(PStream::ipstream& in, ProductMPO& op)
{
   in >> op.data();
   return in;
}

ProductMPO
join(ProductMPO const& Op1, ProductMPO const& Op2)
{
   if (Op1.is_null())
      return Op2;
   if (Op2.is_null())
      return Op1;
   CHECK_EQUAL(Op1.Basis2(), Op2.Basis1());
   ProductMPO Result(Op1.size() + Op2.size());
   for (int i = 0; i < Op1.size(); ++i)
   {
      Result[i] = Op1[i];
   }
   for (int i = 0; i < Op2.size(); ++i)
   {
      Result[i+Op1.size()] = Op2[i];
   }
   return Result;
}

ProductMPO
repeat(ProductMPO const& Op, int Count)
{
   CHECK_EQUAL(Op.Basis1(), Op.Basis2());
   ProductMPO Result(Op.size()*Count);
   for (int i = 0; i < Result.size(); ++i)
   {
      Result[i] = Op[i%Op.size()];
   }
   return Result;
}

ProductMPO prod(ProductMPO const& x, ProductMPO const& y)
{
   return x*y;
}

ProductMPO
operator*(ProductMPO const& x, ProductMPO const& y)
{
   if (x.empty())
      return x;
   if (y.empty())
      return y;

   // we handle also the case where x.size() != y.size()
   int sz =  statistics::lcm(x.size(), y.size());
   ProductMPO Result(sz);
   for (int i = 0; i < sz; ++i)
   {
      Result[i] = aux_tensor_prod(x[i%x.size()], y[i%y.size()]);
   }
   return Result;
}

ProductMPO& operator*=(ProductMPO& x, ProductMPO const& y)
{
   x = x*y;
   return x;
}

ProductMPO inner(ProductMPO const& x, ProductMPO const& y)
{
   return x*y;
}

ProductMPO outer(ProductMPO const& x, ProductMPO const& y)
{
   return x*y;
}

ProductMPO
pow(ProductMPO const& x, int n)
{
   if (n == 0)
   {
      return ProductMPO::make_identity(x.data().LocalBasis1List());
   }
   else if (n%2 == 0)
   {
      return pow(x*x, n/2);
   }
   else if (n == 1)
   {
      return x;
   }
   else
   {
      // n is odd
      return x*pow(x*x, (n-1)/2);
   }
}

ProductMPO
conj(ProductMPO const& x)
{
   ProductMPO Result(x);
   for (ProductMPO::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = conj(*I);
   }
   return Result;
}

ProductMPO
adjoint(ProductMPO const& x)
{
   ProductMPO Result(x);
   for (ProductMPO::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = adjoint(*I);
   }
   return Result;
}

ProductMPO
inv_adjoint(ProductMPO const& x)
{
   ProductMPO Result(x);
   for (ProductMPO::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = inv_adjoint(*I);
   }
   return Result;
}

// For the right to left product, we want
// .... * A(2)B(3) * A(1)B(2) * A(0)B(1) * ...
// or, as a quantum circuit,
//
//     A B | | |
//     | A B | |
//     | | A B |
//     | | | A B ...
// so that our unit cell is the tensor product B \otimes A
ProductMPO
prod_unit_right_to_left(BasicFiniteMPO const& Op, int UnitCellSize)
{
   CHECK(Op.size() % UnitCellSize == 0)("prod_unit: Operator must be a multiple of the unit cell size!")
      (Op.size())(UnitCellSize);
   CHECK_EQUAL(Op.Basis1(), Op.Basis2());

   ProductMPO Result(UnitCellSize);

   // initialize the first unit cell
   for (int i = 0; i < UnitCellSize; ++i)
   {
      Result[i] = Op[i];
   }

   // remaining unit cells
   for (int i = UnitCellSize; i < Op.size(); ++i)
   {
      Result[i%UnitCellSize] = aux_tensor_prod(Op[i], Result[i%UnitCellSize]);
   }

   return Result;
}

// For the left to right product, we want
//  ... * A(0)B(1) * A(1)B(2) * A(2)B(3) * ...
// or, as a quantum circuit,
//
//     | | | | A B ...
//     | | | A B |
//     | | A B | |
//     | A B | | |
// ... A B | | | |
// so that our unit cell is the tensor product A \otimes B
ProductMPO
prod_unit_left_to_right(BasicFiniteMPO const& Op, int UnitCellSize)
{
   CHECK(Op.size() % UnitCellSize == 0)("prod_unit: Operator must be a multiple of the unit cell size!")
      (Op.size())(UnitCellSize);
   CHECK_EQUAL(Op.Basis1(), Op.Basis2());

   ProductMPO Result(UnitCellSize);

   // initialize the first unit cell
   for (int i = 0; i < UnitCellSize; ++i)
   {
      Result[i] = Op[i];
   }

   // remaining unit cells
   for (int i = UnitCellSize; i < Op.size(); ++i)
   {
      Result[i%UnitCellSize] = aux_tensor_prod(Result[i%UnitCellSize], Op[i]);
   }

   return Result;
}

ProductMPO translate_right(std::vector<BasisList> const& LocalBasis)
{
   ProductMPO Result(LocalBasis.size());
   Result[0] = translate_right(LocalBasis.back(), LocalBasis.front());
   for (unsigned i = 1; i < LocalBasis.size(); ++i)
   {
      Result[i] = translate_right(LocalBasis[i-1], LocalBasis[i]);
   }
   return Result;
}
