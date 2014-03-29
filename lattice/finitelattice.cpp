// -*- C++ -*- $Id$

#include "unitcell.h"

LatticeSite
flip_conj(LatticeSite const& A)
{
   LatticeSite Result;

   if (A.empty())
      return Result;

   SiteBasis ReflectedBasis = adjoint(A.Basis1());
   for (LatticeSite::const_iterator ai = A.begin(); ai != A.end(); ++ai)
   {
      Result[ai->first] = flip_conj(ai->second, ReflectedBasis);
   }
   return Result;
}

//
// UnitCell members
//

UnitCell::UnitCell()
{
}

UnitCell::UnitCell(LatticeSite const& s)
   : Data_(new UnitCellType(s))
{
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t)
   : Data_(new UnitCellType(s))
{
   Data_.mutate()->push_back(t);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(new UnitCellType(s))
{
   pvalue_lock<UnitCellType> Lock(Data_);
   Lock->push_back(t);
   Lock->push_back(u);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u, LatticeSite const& v)
   : Data_(new UnitCellType(s))
{
   pvalue_lock<UnitCellType> Lock(Data_);
   Lock->push_back(t);
   Lock->push_back(u);
   Lock->push_back(v);
}

LatticeSite CoerceSL(SymmetryList const& sl, LatticeSite const& s)
{
   LatticeSite r(s);
   r.CoerceSymmetryList(sl);
   return r;
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s)
   : Data_(new UnitCellType(CoerceSL(sl,s)))
{
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t)
   : Data_(new UnitCellType(CoerceSL(sl, s)))
{
   Data_.mutate()->push_back(CoerceSL(sl, t));
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(new UnitCellType(CoerceSL(sl, s)))
{
   pvalue_lock<UnitCellType> Lock(Data_);
   Lock->push_back(CoerceSL(sl, t));
   Lock->push_back(CoerceSL(sl, u));
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, 
                 LatticeSite const& u, LatticeSite const& v)
   : Data_(new UnitCellType(CoerceSL(sl, s)))
{
   pvalue_lock<UnitCellType> Lock(Data_);
   Lock->push_back(CoerceSL(sl, t));
   Lock->push_back(CoerceSL(sl, u));
   Lock->push_back(CoerceSL(sl, v));
}

UnitCell::UnitCell(int RepeatCount, UnitCell const& l)
   : Data_(new UnitCellType(RepeatCount, l.data()))
{
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2)
   : Data_(x1.Data_)
{
   Data_.mutate()->push_back(*x2.Data_);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3)
   : Data_(x1.Data_)
{
   pvalue_lock<UnitCellType> Lock(Data_);
   Lock->push_back(*x2.Data_);
   Lock->push_back(*x3.Data_);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3, UnitCell const& x4)
   : Data_(x1.Data_)
{
   pvalue_lock<UnitCellType> Lock(Data_);
   Lock->push_back(*x2.Data_);
   Lock->push_back(*x3.Data_);
   Lock->push_back(*x4.Data_);
}

UnitCell::UnitCell(int Size, LatticeSite const& s)
   : Data_(new UnitCellType(Size, run_length_compressed<LatticeSite>(s)))
{
}

UnitCell::UnitCell(LatticeSite const& s, std::string const& Coord)
   : Data_(new UnitCellType(s))
{
}

LatticeSite const& 
UnitCell::operator[](int n) const
{
   CHECK(n >= 1 && n <= this->size())(n)(this->size());
   const_iterator I = this->begin();
   std::advance(I, n-1);
   return *I;
}

PStream::opstream& operator<<(PStream::opstream& out, UnitCell const& L)
{
   return out << L.Data_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, UnitCell& L)
{
   in >>  L.Data_;
   return in;
}

UnitCell repeat(UnitCell const& x, int RepeatCount)
{
   return UnitCell(RepeatCount, x);
}

UnitCell join(UnitCell const& x, UnitCell const& y)
{
   return UnitCell(x, y);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z)
{
   return UnitCell(x,y,z);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w)
{
   return UnitCell(x,y,z,w);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w,
	     UnitCell const& v)
{
   return join(UnitCell(x,y,z,w), v);
}

bool
operator==(UnitCell const& u1, UnitCell const& u2)
{
   if (u1.Data_ == u2.Data_)
      return true;
   // else
   if (u1.Data_.is_null() || u2.Data_.is_null())
      return false;
   // else
   if (u1.size() != u2.size())
      return false;
   // else
   // Here we need to do a deep comparison of the unit cells.  Since we don't have comparison
   // on LatticeSite objects (which could be added but is relatively expensive) we'll cheat
   // and just compare the local hilbert space.
   int Sz = u1.size();
   for (int i = 0; i < Sz; ++i)
   {
      if (u1[i].Basis1() != u2[i].Basis1() || u1[i].Basis2() != u2[i].Basis2())
	 return false;
   }
   return true;
}

bool
operator!=(UnitCell const& u1, UnitCell const& u2)
{
   return !(u1 == u2);
}

std::ostream&
operator<<(std::ostream& out, UnitCell const& u)
{
   out << "UnitCell object size " << u.size() << "\n";
   return out;
}
