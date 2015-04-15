// -*- C++ -*- $Id$

#include "unitcell_mpo.h"
#include "unitcell.h"
#include "common/statistics.h"

UnitCellMPO::UnitCellMPO()
   : Cell(NULL)
{
}

UnitCellMPO::UnitCellMPO(UnitCell const* Cell_, FiniteMPO const& Op_, LatticeCommute Com_, int Offset_)
   : Cell(Cell_), Op(Op_), Com(Com_), Offset(Offset_)
{
}

#if 0
UnitCellMPO::UnitCellMPO(UnitCell const* Cell_, int Size, LatticeCommute Com_, int Offset_)
   : Cell(Cell_), Op(Size, Com_), Com(Com_), Offset(Offset_)
{
}
#endif

PStream::opstream& operator<<(PStream::opstream& out, UnitCellMPO const& L)
{
   out << L.Op << L.Com << L.Offset;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, UnitCellMPO& L)
{
   in >> L.Op >> L.Com >> L.Offset;
   L.Cell = NULL;
   return in;
}

std::ostream& operator<<(std::ostream& out, UnitCellMPO const& Op)
{
   out << "Unit cell operator starts at offset " << Op.offset() << ", size " << Op.size() << '\n';
   out << "Commutes: " << Op.Commute() << '\n';
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
      Op = join(repeat(Cell->string_mpo(Com, Op.qn1()), (Offset-OtherOffset) / Cell->size()), Op);
      Offset = OtherOffset;
   }

   // do we need to extend the operator on the right?
   if (Offset+Op.size() < OtherOffset+OtherSize)
   {
      Op = join(Op, repeat(Cell->identity_mpo(Op.qn2()), (OtherOffset+OtherSize-Offset-Op.size())/Cell->size()));
   }
}

void
UnitCellMPO::ExtendToCoverUnitCell(int OtherSize)
{
   int NewOffset = Offset;
   // on the left hand side, we want to decrease the offset to a multiple of OtherSize
   if (Offset % OtherSize != 0)
   {
      if (Offset < 0)
      {
	 // Offset is negative - increase it to a multiple of OtherSize
	 NewOffset = -statistics::lcm(-Offset, OtherSize);
      }
      else
      {
	 // Offset is positive - reduce it to be a multiple of OtherSize
	 NewOffset = Offset - (Offset%OtherSize);
      }
   }
   int NewSize = statistics::lcm(Op.size(), OtherSize);
   this->ExtendToCover(NewSize, NewOffset);
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
   return UnitCellMPO(&xCopy.GetUnitCell(), xCopy.MPO()+yCopy.MPO(), xCopy.Commute(), xCopy.offset());
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
   return UnitCellMPO(&xCopy.GetUnitCell(), xCopy.MPO()-yCopy.MPO(), xCopy.Commute(), xCopy.offset());
}

UnitCellMPO operator-(UnitCellMPO const& x)
{
   if (x.is_null())
      return x;
   return UnitCellMPO(&x.GetUnitCell(), -x.MPO(), x.Commute(), x.offset());
}

UnitCellMPO operator*(double a, UnitCellMPO const& x)
{
   return UnitCellMPO(&x.GetUnitCell(), a*x.MPO(), x.Commute(), x.offset());
}

UnitCellMPO operator*(UnitCellMPO const& x, double a)
{
   return UnitCellMPO(&x.GetUnitCell(), x.MPO()*a, x.Commute(), x.offset());
}

UnitCellMPO operator*(std::complex<double> a, UnitCellMPO const& x)
{
   return UnitCellMPO(&x.GetUnitCell(), a*x.MPO(), x.Commute(), x.offset());
}

UnitCellMPO operator*(UnitCellMPO const& x, std::complex<double> a)
{
   return UnitCellMPO(&x.GetUnitCell(), x.MPO()*a, x.Commute(), x.offset());
}

UnitCellMPO prod(UnitCellMPO const& x, UnitCellMPO const& y, QuantumNumbers::QuantumNumber const& q)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(&xCopy.GetUnitCell(), prod(xCopy.MPO(), yCopy.MPO(), q), 
		      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO prod(UnitCellMPO const& x, UnitCellMPO const& y)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(&xCopy.GetUnitCell(), prod(xCopy.MPO(), yCopy.MPO()), 
		      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO operator*(UnitCellMPO const& x, UnitCellMPO const& y)
{
   return prod(x,y);
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(&xCopy.GetUnitCell(), prod(xCopy.MPO(), yCopy.MPO()), 
		      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO dot(UnitCellMPO const& x, UnitCellMPO const& y)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(&xCopy.GetUnitCell(), dot(xCopy.MPO(), yCopy.MPO()), 
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
   return UnitCellMPO(&xCopy.GetUnitCell(), cross(xCopy.MPO(), yCopy.MPO()), 
		      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

UnitCellMPO outer(UnitCellMPO const& x, UnitCellMPO const& y)
{
   UnitCellMPO xCopy(x);
   xCopy.ExtendToCover(y.size(), y.offset());
   UnitCellMPO yCopy(y);
   yCopy.ExtendToCover(x.size(), x.offset());
   return UnitCellMPO(&xCopy.GetUnitCell(), outer(xCopy.MPO(), yCopy.MPO()), 
		      xCopy.Commute()*yCopy.Commute(), xCopy.offset());
}

// project a (reducible) operator onto an irreducible component
UnitCellMPO project(UnitCellMPO const& x, QuantumNumbers::QuantumNumber const& q)
{
   return UnitCellMPO(&x.GetUnitCell(), project(x.MPO(), q), x.Commute(), x.offset());
}

UnitCellMPO pow(UnitCellMPO const& x, int n)
{
   return UnitCellMPO(&x.GetUnitCell(), pow(x.MPO(), n), x.Commute()*n, x.offset());
}

// Exponential operator.
UnitCellMPO exp(UnitCellMPO const& x)
{
   CHECK_EQUAL(x.Commute(), LatticeCommute::Bosonic)("Operator must be bosonic to calculate the exponential");
   return UnitCellMPO(&x.GetUnitCell(), exp(x.MPO()), LatticeCommute::Bosonic, x.offset());
}

// Conjugate
UnitCellMPO conj(UnitCellMPO const& x)
{
   return UnitCellMPO(&x.GetUnitCell(), conj(x.MPO()), x.Commute(), x.offset());
}

// Adjoint
UnitCellMPO adjoint(UnitCellMPO const& x)
{
   return UnitCellMPO(&x.GetUnitCell(), adjoint(x.MPO()), x.Commute(), x.offset());
}

UnitCellMPO MakeIdentityFrom(UnitCellMPO const& x)
{
   return UnitCellMPO(&x.GetUnitCell(), x.GetUnitCell().identity_mpo(), LatticeCommute::Bosonic, x.offset());
}
