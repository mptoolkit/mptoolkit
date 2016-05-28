// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/dense.cc
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

inline
DenseOperator::DenseOperator()
{
}

inline
DenseOperator::DenseOperator(VectorBasis const& b1, VectorBasis const& b2, QuantumNumber const& q)
   : Basis1_(b1), Basis2_(b2), Data_(b,q)
{
}

inline
DenseOperator::DenseOperator(VectorBasis const& b1, VectorBasis const& b2, data_type const& d)
   : Basis1_(b1), Basis2_(b2), SparsePart_(b, sdTransformsAs()), Data_(d)
{
   DEBUG_PRECONDITION(d.Basis1() == b1);
   DEBUG_PRECONDITION(d.Basis2() == b2);
}

inline
DenseOperator& DenseOperator::operator+=(DenseOperator const& x)
{
   Data_ += x.Data_;
   return *this;
}

DenseOperator& DenseOperator::operator-=(DenseOperator const& x)
{
   Data_ -= x.Data_;
   return *this;
}

inline
DenseOperator& DenseOperator::operator*=(value_type a)
{
   Dta_ *= a;
   return *this;
}

inline
bool operator==(DenseOperator const& x, DenseOperator const& y)
{
   return x.data() == y.data();
}

inline
bool operator!=(DenseOperator const& x, DenseOperator const& y)
{
   return x.data() != y.data();
}

inline
bool equal(DenseOperator const& x, DenseOperator const& y, double tol)
{
   return equal(x.data(), y.data(), tol);
}

inline
DenseOperator operator-(DenseOperator const& x)
{
   return DenseOperator(x.Basis1(), x.Basis2(), -x.data());
}

inline
DenseOperator operator+(DenseOperator const& x, DenseOperator const& y)
{
   DEBUG_PRECONDITION(x.Basis() == y.Basis());
   DenseOperator Result(x.Basis1(), x.Basis2(), x.data() + y.data());
   return Result;
}

inline
DenseOperator operator-(DenseOperator const& x, DenseOperator const& y)
{
   DEBUG_PRECONDITION(x.Basis() == y.Basis());
   DenseOperator Result(x.Basis1(), x.Basis2(), x.data() - y.data());
   return Result;
}

inline
DenseOperator::value_type trace(DenseOperator const& x)
{
   return trace(x.data());
}

inline
double norm_2_sq(DenseOperator const& x)
{
   return norm_2_sq(x.data());
}

inline
double norm_2(DenseOperator const& x)
{
   return std::sqrt(norm_2_sq(x));
}

inline
DenseOperator operator*(double a, DenseOperator const& x)
{
   return DenseOperator(x.Basis(), a * x.data());
}

inline
DenseOperator operator*(DenseOperator const& x, double a)
{
   return DenseOperator(x.Basis(), x.data() * a);
}

inline
DenseOperator::value_type
scalar_prod(DenseOperator const& x, DenseOperator const& y)
{
   return LinearAlgebra::scalar_prod(x,y);
}

inline
DenseOperator prod(DenseOperator const& x, DenseOperator const& y, 
		   QuantumNumber const& q = QuantumNUmber())
{
   DenseOperator Result(x.Basis1(), y.Basis2(),
			prod(x.data(), y.data(), q));
}

inline
HermitianProxy<DenseOperator const> herm(DenseOperator const& x)
{
   return HermitianProxy<DenseOperator const>(x);
}

inline
DenseOperator scalar_product(DenseOperator const& x, HermitianProxy<DenseOperator const> y)
{
   return DenseOperator(x.Basis1(), y.base().Basis1(), 
			scalar_product(x.data(), herm(y.base().data())));
}

inline
DenseOperator scalar_product(HermitianProxy<DenseOperator const> x, DenseOperator const& y)
{
   return DenseOperator(x.base().Basis2(), y.Basis2(),
			scalar_product(herm(x.base().data()), y.data()));
}

inline
DenseOperator adjoint(DenseOperator const& x)
{
   return DenseOperator(x.Basis2(), x.Basis1(), adjoint(x.data()));
}

inline
DenseOperator inv_adjoint(DenseOperator const& x)
{
   return DenseOperator(x.Basis2(), x.Basis1(), inv_adjoint(x.data()));
}
