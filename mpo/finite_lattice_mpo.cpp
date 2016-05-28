// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/finite_lattice_mpo.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "finite_lattice_mpo.h"

FiniteMPO
identity_mpo(UnitCell const& c)
{
   FiniteMPO Result(c.size());
   QuantumNumbers::QuantumNumber Ident(c.GetSymmetryList());
   BasisList b(c.GetSymmetryList());
   b.push_back(Ident);
   for (int i = 0; i < c.size(); ++i)
   {
      Result[i] = OperatorComponent(c.LocalBasis(i), b, b);
      Result[i](0,0) = c[i]["I"];
   }
   return Result;
}

FiniteMPO
identity_mpo(UnitCell const& c, int Size)
{
   if (Size == 0)
      return FiniteMPO();

   PRECONDITION_EQUAL(Size % c.size(), 0);
   return repeat(identity_mpo(c), Size / c.size());
}

FiniteLatticeMPO::FiniteLatticeMPO()
{
}

FiniteLatticeMPO::FiniteLatticeMPO(FiniteLatticeMPO const& Other)
   :
   UnitCell_(Other.UnitCell_),
   JWString_(Other.JWString_),
   Operator_(Other.Operator_)
{
}

FiniteLatticeMPO&
FiniteLatticeMPO::operator=(FiniteLatticeMPO const& Other)
{
   UnitCell_ = Other.UnitCell_;
   JWString_ = Other.JWString_;
   Operator_ = Other.Operator_;
   return *this;
}

int
FiniteLatticeMPO::size() const
{
   return Operator_.size();
}

FiniteMPO
FiniteLatticeMPO::JWString(int Size) const
{
   PRECONDITION_EQUAL(Size % UnitCell_.size(), 0);
   return repeat(JWString_, Size / UnitCell_.size());
}

FiniteMPO
FiniteLatticeMPO::ApplyJW(FiniteLatticeMPO const& f)
{
   return Operator_ * f.JWString(Operator_.size());
}

FiniteMPO
FiniteLatticeMPO::AsFiniteMPO(int Size) const
{
   PRECONDITION_EQUAL(Size % UnitCell_.size(), 0);
   return join(Operator_, identity_mpo(UnitCell_, Size - Operator_.size()));
}

PStream::opstream&
operator<<(PStream::opstream& out, FiniteLatticeMPO const& op)
{
   out << op.UnitCell_;
   out << op.JWString_;
   out << op.Operator_;
   return out;
}

PStream::ipstream&
operator>>(PStream::ipstream& in, FiniteLatticeMPO& op)
{
   in >> op.UnitCell_;
   in >> op.JWString_;
   in >> op.Operator_;
   return in;
}

FiniteLatticeMPO& operator*=(FiniteLatticeMPO& x, double a)
{
   x.Operator_ *= a;
   return x;
}

FiniteLatticeMPO& operator*=(FiniteLatticeMPO& x, std::complex<double> a)
{
   x.Operator_ *= a;
   return x;
}

FiniteLatticeMPO& operator+=(FiniteLatticeMPO& x, FiniteLatticeMPO const& y)
{
   DEBUG_CHECK_EQUAL(x.JWString(), y.JWString());
   int Size = std::max(x.size(), y.size());
   x.Operator_ = x.AsFiniteMPO(Size) + y.AsFiniteMPO(Size);
   return x;
}

FiniteLatticeMPO& operator-=(FiniteLatticeMPO& x, FiniteLatticeMPO const& y)
{
   DEBUG_CHECK_EQUAL(x.JWString(), y.JWString());
   int Size = std::max(x.size(), y.size());
   x.Operator_ = x.AsFiniteMPO(Size) - y.AsFiniteMPO(Size);
   return x;
}

FiniteLatticeMPO operator+(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y)
{
   DEBUG_CHECK_EQUAL(x.JWString(), y.JWString());
   int Size = std::max(x.size(), y.size());
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), x.AsFiniteMPO(Size) + y.AsFiniteMPO(Size));
}

FiniteLatticeMPO operator-(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y)
{
   DEBUG_CHECK_EQUAL(x.JWString(), y.JWString());
   int Size = std::max(x.size(), y.size());
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), x.AsFiniteMPO(Size) - y.AsFiniteMPO(Size));
}

FiniteLatticeMPO operator-(FiniteLatticeMPO const& x)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), -x.AsFiniteMPO());
}

FiniteLatticeMPO operator*(double a, FiniteLatticeMPO const& x)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), a * x.AsFiniteMPO());
}   

FiniteLatticeMPO operator*(FiniteLatticeMPO const& x, double a)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), x.AsFiniteMPO() * a);
}   


FiniteLatticeMPO operator*(std::complex<double> a, FiniteLatticeMPO const& x)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), a * x.AsFiniteMPO());
}   


FiniteLatticeMPO operator*(FiniteLatticeMPO const& x, std::complex<double> a)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), x.AsFiniteMPO() * a);
}   

FiniteLatticeMPO prod(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y, 
		      QuantumNumbers::QuantumNumber const& q)
{
   CHECK_EQUAL(x.GetUnitCell(), y.GetUnitCell());
   int Size = std::max(x.size(), y.size());
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString() * y.JWString(), 
			   prod(x.AsFiniteMPO(Size), y.AsFiniteMPO(Size), q));
}

FiniteLatticeMPO prod(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y)
{
   CHECK_EQUAL(x.GetUnitCell(), y.GetUnitCell());
   int Size = std::max(x.size(), y.size());
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString() * y.JWString(), 
			   prod(x.AsFiniteMPO(Size), y.AsFiniteMPO(Size)));
}

FiniteLatticeMPO operator*(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y)
{
   CHECK_EQUAL(x.GetUnitCell(), y.GetUnitCell());
   int Size = std::max(x.size(), y.size());
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString() * y.JWString(), 
			   prod(x.AsFiniteMPO(Size), y.AsFiniteMPO(Size)));
}

FiniteLatticeMPO dot(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y)
{
   CHECK_EQUAL(x.GetUnitCell(), y.GetUnitCell());
   int Size = std::max(x.size(), y.size());
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString() * y.JWString(), 
			   dot(x.AsFiniteMPO(Size), y.AsFiniteMPO(Size)));
}

FiniteLatticeMPO project(FiniteLatticeMPO const& x, QuantumNumbers::QuantumNumber const& q)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), project(x.AsFiniteMPO(), q));
}

FiniteLatticeMPO pow(FiniteLatticeMPO const& x, int n)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), pow(x.AsFiniteMPO(), n));
}
   
FiniteLatticeMPO conj(FiniteLatticeMPO const& x)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), conj(x.AsFiniteMPO()));
}

FiniteLatticeMPO adjoint(FiniteLatticeMPO const& x)
{
   return FiniteLatticeMPO(x.GetUnitCell(), x.JWString(), adjoint(x.AsFiniteMPO()));
}
