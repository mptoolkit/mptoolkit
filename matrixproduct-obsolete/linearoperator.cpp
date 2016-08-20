// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/linearoperator.cpp
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

#include "linearoperator.h"
#include "pstream/variant.h"
#include <algorithm>

using QuantumNumbers::QuantumNumber;

bool LinearOperator::is_irreducible() const
{
   //   return (!this->is_null() && this->Basis1().size() == 1);
   return (this->is_null() || this->Basis1().size() == 1);
}

QuantumNumbers::QuantumNumber
LinearOperator::TransformsAs() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   CHECK(this->is_irreducible());
   return this->Basis1()[0];
}

PStream::opstream& operator<<(PStream::opstream& out, LinearOperator const& op)
{
   return out << op.data();
}

PStream::ipstream& operator>>(PStream::ipstream& in, LinearOperator& op)
{
   return in >> op.data();
}

LinearOperator& operator*=(LinearOperator& x, double a)
{
   SimpleOperator c = SimpleOperator::make_identity(x.Basis1());
   c *= a;
   x.data() = multiply_left(c, x.data());
   return x;
}

LinearOperator& operator*=(LinearOperator& x, std::complex<double> a)
{
   SimpleOperator c = SimpleOperator::make_identity(x.Basis1());
   c *= a;
   x.data() = multiply_left(c, x.data());
   return x;
}

LinearOperator& operator+=(LinearOperator& x, LinearOperator const& y)
{
   x = x+y;
   return x;
}

LinearOperator& operator-=(LinearOperator& x, LinearOperator const& y)
{
   x = x + (-1.0 * y);
   return x;
}

LinearOperator operator+(LinearOperator const& x, LinearOperator const& y)
{
   if (x.empty())
      return y;
   if (y.empty())
      return x;

   SumBasis<BasisList> sb(x.Basis1(), y.Basis1());
   SimpleOperator C = CollapseBasis(sb.Basis());
   MPOpCompressed Temp = do_sum(C, sb, x.data(), y.data());
   C = C * herm(CollapseBasis(C.Basis2()));
   Temp = inject_right(Temp, C);
   C = RemoveEmptyRows(C);
   Temp = multiply_left(C, Temp);
   return LinearOperator(Temp);
}

LinearOperator operator-(LinearOperator const& x, LinearOperator const& y)
{
   return x + (-1.0 * y);
}

LinearOperator operator-(LinearOperator const& x)
{
   return -1.0 * x;
}

LinearOperator operator*(double a, LinearOperator const& x)
{
   LinearOperator Result(x);
   Result.data() = multiply_left(a * SimpleOperator::make_identity(Result.Basis1()), Result.data());
   return Result;
}

LinearOperator operator*(LinearOperator const& x, double a)
{
   LinearOperator Result(x);
   Result.data() = multiply_left(a * SimpleOperator::make_identity(Result.Basis1()), Result.data());
   return Result;
}

LinearOperator operator*(std::complex<double> a, LinearOperator const& x)
{
   LinearOperator Result(x);
   Result.data() = multiply_left(a * SimpleOperator::make_identity(Result.Basis1()), Result.data());
   return Result;
}

LinearOperator operator*(LinearOperator const& x, std::complex<double> a)
{
   LinearOperator Result(x);
   Result.data() = multiply_left(a * SimpleOperator::make_identity(Result.Basis1()), Result.data());
   return Result;
}

LinearOperator prod(LinearOperator const& x, LinearOperator const& y, QuantumNumbers::QuantumNumber const& q)
{
   if (x.is_null())
      return x;
   else if (y.is_null())
      return y;
   // else
   ProductBasis<BasisList, BasisList> pb(x.Basis1(), y.Basis1(), q);
   SimpleOperator C = CollapseBasis(pb.Basis());
   MPOpCompressed Temp = do_prod(C, pb, x.data(), y.data());
   C = C * herm(CollapseBasis(C.Basis2()));
   Temp = inject_right(Temp, C);
   C = RemoveEmptyRows(C);
   Temp = multiply_left(C, Temp);
   return LinearOperator(Temp);
}

LinearOperator prod(LinearOperator const& x, LinearOperator const& y)
{
   if (x.is_null())
      return x;
   else if (y.is_null())
      return y;
   // else
   ProductBasis<BasisList, BasisList> pb(x.Basis1(), y.Basis1());
   SimpleOperator C = CollapseBasis(pb.Basis());
   MPOpCompressed Temp = do_prod(C, pb, x.data(), y.data());
   C = C * herm(CollapseBasis(C.Basis2()));
   Temp = inject_right(Temp, C);
   C = RemoveEmptyRows(C);
   Temp = multiply_left(C, Temp);
   return LinearOperator(Temp);
}

LinearOperator operator*(LinearOperator const& x, LinearOperator const& y)
{
   return prod(x,y);
}

LinearOperator dot(LinearOperator const& x, LinearOperator const& y)
{
   CHECK_EQUAL(x.TransformsAs(), adjoint(y.TransformsAs()));
   double f = std::sqrt(double(degree(x.TransformsAs())));
   return f*prod(x,y,QuantumNumber(x.GetSymmetryList()));
}

LinearOperator project(LinearOperator const& x, QuantumNumbers::QuantumNumber const& q)
{
   unsigned i = std::distance(x.Basis1().begin(), std::find(x.Basis1().begin(), x.Basis1().end(), q));
   if (i == x.Basis1().size())  // in this case, q does not exist in the basis anyway
      return LinearOperator();

   BasisList NewBasis(x.GetSymmetryList());
   NewBasis.push_back(x.Basis1()[i]);
   SimpleOperator C(NewBasis, x.Basis1());
   C(0,i) = 1.0;
   MPOpCompressed Temp = inject_left(C, x.data());
   C = C * herm(CollapseBasis(C.Basis2()));
   Temp = inject_right(Temp, C);
   C = RemoveEmptyRows(C);
   Temp = multiply_left(C, Temp);
   return LinearOperator(Temp);
}

LinearOperator pow(LinearOperator const& x, int n)
{
   // a binary algorithm could do this better
   CHECK(n >= 1);
   LinearOperator Result = x;
   while (n > 1)
   {
      Result = Result * x;
      --n;
   }
   return Result;
}

LinearOperator conj(LinearOperator const& x)
{
   return LinearOperator(conj(x.data()));
}

LinearOperator adjoint(LinearOperator const& x)
{
   return LinearOperator(adjoint(x.data()));
}
