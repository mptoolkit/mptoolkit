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
   : Data_(s)
{
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t)
   : Data_(s)
{
   Data_.push_back(t);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(s)
{
   Data_.push_back(t);
   Data_.push_back(u);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u, LatticeSite const& v)
   : Data_(s)
{
   Data_.push_back(t);
   Data_.push_back(u);
   Data_.push_back(v);
}

LatticeSite CoerceSL(SymmetryList const& sl, LatticeSite const& s)
{
   LatticeSite r(s);
   r.CoerceSymmetryList(sl);
   return r;
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s)
   : Data_(CoerceSL(sl,s))
{
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t)
   : Data_(CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, 
                 LatticeSite const& u, LatticeSite const& v)
   : Data_(CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
   Data_.push_back(CoerceSL(sl, v));
}

UnitCell::UnitCell(int RepeatCount, UnitCell const& l)
   : Data_(RepeatCount, l.data())
{
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2)
   : Data_(x1.Data_)
{
   Data_.push_back(x2.Data_);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3)
   : Data_(x1.Data_)
{
   Data_.push_back(x2.Data_);
   Data_.push_back(x3.Data_);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3, UnitCell const& x4)
   : Data_(x1.Data_)
{
   Data_.push_back(x2.Data_);
   Data_.push_back(x3.Data_);
   Data_.push_back(x4.Data_);
}

UnitCell::UnitCell(int Size, LatticeSite const& s)
   : Data_(Size, run_length_compressed<LatticeSite>(s))
{
}

UnitCell::UnitCell(LatticeSite const& s, std::string const& Coord)
   : Data_(s)
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

