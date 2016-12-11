// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/unitcell_mpo.cpp
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

#include "unitcell_mpo.h"
#include "unitcell.h"
#include "common/statistics.h"

UnitCellMPO::UnitCellMPO(SiteListPtrType const& SiteList_, FiniteMPO Op_, LatticeCommute Com_, int Offset_,
                         std::string Description_, int CoarseGrain_)
   : SiteList(SiteList_), Op(std::move(Op_)), Com(Com_), Offset(Offset_), Description(std::move(Description_)), CoarseGrain(CoarseGrain_)
{
}

UnitCellMPO&
UnitCellMPO::operator=(UnitCellMPO const& c)
{
   SiteList = c.SiteList;
   Op = c.Op;
   Com = c.Com;
   Offset = c.Offset;
   if (Description.empty())
      Description = c.Description;
   CoarseGrain = c.CoarseGrain;
}

UnitCellMPO&
UnitCellMPO::operator=(UnitCellMPO&& c)
{
   SiteList = std::move(c.SiteList);
   Op = std::move(c.Op);
   Com = c.Com;
   Offset = c.Offset;
   if (Description.empty())
      Description = std::move(c.Description);
   CoarseGrain = c.CoarseGrain;
   return *this;
}

extern PStream::VersionTag LatticeVersion;

PStream::opstream& operator<<(PStream::opstream& out, UnitCellMPO const& L)
{
   out << L.SiteList << L.Op << L.Com << L.Offset << L.Description << L.CoarseGrain;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, UnitCellMPO& L)
{
   in >> L.SiteList >> L.Op >> L.Com >> L.Offset;
   if (in.version_of(LatticeVersion) >= 3)
   {
      in >> L.Description;
   }
   else
   {
      L.Description = "";
   }
   if (in.version_of(LatticeVersion) >= 4)
   {
      in >> L.CoarseGrain;
   }
   return in;
}

std::ostream& operator<<(std::ostream& out, UnitCellMPO const& Op)
{
   out << "Unit cell operator starts at offset " << Op.offset() << ", size " << Op.size() << '\n';
   out << "Commutes: " << Op.Commute() << '\n';
   out << "Coarse grain: " << Op.coarse_grain_factor() << '\n';
   out << Op.MPO();
   return out;
}

UnitCellMPO& operator*=(UnitCellMPO& x, double a)
{
   x.MPO() *= a;
   return x;
}

UnitCellMPO& operator*=(UnitCellMPO& x, std::complex<double> a)
{
   x.MPO() *= a;
   return x;
}

void
UnitCellMPO::ExtendToCover(int OtherSize, int OtherOffset)
{
   // do we need to extend the operator on the left?
   if (Offset > OtherOffset)
   {
      // need to extend this operator at the front with JW strings
      Op = join(repeat(string_mpo(*SiteList, Com.SignOperator(), Op.qn1()),
                       (Offset-OtherOffset) / SiteList->size()), Op);
      Offset = OtherOffset;
   }

   // do we need to extend the operator on the right?
   if (Offset+Op.size() < OtherOffset+OtherSize)
   {
      Op = join(Op, repeat(identity_mpo(*SiteList, Op.qn2()),
                           (OtherOffset+OtherSize-Offset-Op.size())/SiteList->size()));
   }
}

// round n up to the nearest multiple of m, assuming m positive
int round_up(int n, int m)
{
   int nn = n - n%m;
   if (nn < n)
      nn += m;
   return nn;
}

// round n down to the nearest multiple of m, assuming m positive
int round_down(int n, int m)
{
   int nn = n - n%m;
   if (nn > n)
      nn -= m;
   return nn;
}

void
UnitCellMPO::ExtendToCoverUnitCell(int OtherSize)
{
   int NewOffset = Offset;
   // on the left hand side, we want to decrease the offset to a multiple of OtherSize
   if (Offset % OtherSize != 0)
   {
      NewOffset = round_down(Offset, OtherSize);
   }
   int NewSize = round_up(Op.size()+(Offset-NewOffset), OtherSize);
   this->ExtendToCover(NewSize, NewOffset);
}

FiniteMPO
UnitCellMPO::GetJWStringUnit()
{
   return string_mpo(*SiteList, Com.SignOperator(), Op.qn1());
}

// Many of these implementations are not the most efficient possible.
// We simpy copy the argument and extend it so we can use the FiniteMPO operations.
// In many cases we could avoid the copy.

UnitCellMPO& operator+=(UnitCellMPO& x, UnitCellMPO const& y)
{
   if (x.is_null())
   {
      x = y;
      return x;
   }
   if (y.is_null())
   {
      return x;
   }
   CHECK_EQUAL(x.Commute(), y.Commute());

   x.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   x.MPO() += yCopy.MPO();
   return x;
}

UnitCellMPO& operator-=(UnitCellMPO& x, UnitCellMPO const& y)
{
   if (x.is_null())
   {
      x = -y;
      return x;
   }
   if (y.is_null())
   {
      return x;
   }
   CHECK_EQUAL(x.Commute(), y.Commute());

   x.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   x.MPO() -= yCopy.MPO();
   return x;
}


UnitCellMPO operator+(UnitCellMPO const& x, UnitCellMPO const& y)
{
   if (x.is_null())
      return y;
   if (y.is_null())
      return x;
   CHECK_EQUAL(x.Commute(), y.Commute());

   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), xCopy.MPO()+yCopy.MPO(), xCopy.Commute(), xCopy.offset());
}

UnitCellMPO operator-(UnitCellMPO const& x, UnitCellMPO const& y)
{
   if (x.is_null())
      return -y;
   if (y.is_null())
      return x;
   CHECK_EQUAL(x.Commute(), y.Commute());

   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), xCopy.MPO()-yCopy.MPO(), xCopy.Commute(), xCopy.offset());
}

UnitCellMPO operator-(UnitCellMPO const& x)
{
   if (x.is_null())
      return x;
   return UnitCellMPO(x.GetSiteList(), -x.MPO(), x.Commute(), x.offset());
}

UnitCellMPO operator*(double a, UnitCellMPO const& x)
{
   return UnitCellMPO(x.GetSiteList(), a*x.MPO(), x.Commute(), x.offset());
}

UnitCellMPO operator*(UnitCellMPO const& x, double a)
{
   return UnitCellMPO(x.GetSiteList(), x.MPO()*a, x.Commute(), x.offset());
}

UnitCellMPO operator*(std::complex<double> a, UnitCellMPO const& x)
{
   return UnitCellMPO(x.GetSiteList(), a*x.MPO(), x.Commute(), x.offset());
}

UnitCellMPO operator*(UnitCellMPO const& x, std::complex<double> a)
{
   return UnitCellMPO(x.GetSiteList(), x.MPO()*a, x.Commute(), x.offset());
}

UnitCellMPO prod(UnitCellMPO const& x, UnitCellMPO const& y, QuantumNumbers::QuantumNumber const& q)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), prod(xCopy.MPO(), yCopy.MPO(), q),
                      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO prod(UnitCellMPO const& x, UnitCellMPO const& y)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), prod(xCopy.MPO(), yCopy.MPO()),
                      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO operator*(UnitCellMPO const& x, UnitCellMPO const& y)
{
   return prod(x,y);
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), prod(xCopy.MPO(), yCopy.MPO()),
                      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO dot(UnitCellMPO const& x, UnitCellMPO const& y)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), dot(xCopy.MPO(), yCopy.MPO()),
                      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO inner(UnitCellMPO const& x, UnitCellMPO const& y)
{
   return dot(adjoint(x),y);
}

UnitCellMPO cross(UnitCellMPO const& x, UnitCellMPO const& y)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), cross(xCopy.MPO(), yCopy.MPO()),
                      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO outer(UnitCellMPO const& x, UnitCellMPO const& y)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(xCopy.GetSiteList(), outer(xCopy.MPO(), yCopy.MPO()),
                      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

// project a (reducible) operator onto an irreducible component
UnitCellMPO project(UnitCellMPO const& x, QuantumNumbers::QuantumNumber const& q)
{
   return UnitCellMPO(x.GetSiteList(), project(x.MPO(), q), x.Commute(), x.offset());
}

UnitCellMPO pow(UnitCellMPO const& x, int n)
{
   return UnitCellMPO(x.GetSiteList(), pow(x.MPO(), n), x.Commute()*n, x.offset());
}

// Exponential operator.
UnitCellMPO exp(UnitCellMPO const& x)
{
   CHECK_EQUAL(x.Commute(), LatticeCommute::Bosonic)("Operator must be bosonic to calculate the exponential");
   return UnitCellMPO(x.GetSiteList(), exp(x.MPO()), LatticeCommute::Bosonic, x.offset());
}

// Conjugate
UnitCellMPO conj(UnitCellMPO const& x)
{
   return UnitCellMPO(x.GetSiteList(), conj(x.MPO()), x.Commute(), x.offset());
}

// Adjoint
UnitCellMPO adjoint(UnitCellMPO const& x)
{
   return UnitCellMPO(x.GetSiteList(), adjoint(x.MPO()), x.Commute(), x.offset());
}

// Inverse Adjoint
UnitCellMPO inv_adjoint(UnitCellMPO const& x)
{
   return UnitCellMPO(x.GetSiteList(), inv_adjoint(x.MPO()), x.Commute(), x.offset());
}

UnitCellMPO MakeIdentityFrom(UnitCellMPO const& x)
{
   return UnitCellMPO(x.GetSiteList(), identity_mpo(*x.GetSiteList()), LatticeCommute::Bosonic, x.offset());
}

void optimize(UnitCellMPO& Op)
{
   optimize(Op.MPO());
}

void qr_optimize(UnitCellMPO& Op)
{
   qr_optimize(Op.MPO());
}

UnitCellMPO translate(UnitCellMPO x, int Sites)
{
   CHECK(Sites % x.unit_cell_size() == 0);
   UnitCellMPO Result(std::move(x));
   Result.translate(Sites);
   return Result;
}
