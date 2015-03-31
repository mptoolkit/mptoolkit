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
   : Data_(1, s)
{
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t)
   : Data_(1, s)
{
   Data_.push_back(t);
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(1, s)
{
   Data_.push_back(t);
   Data_.push_back(u);
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u, LatticeSite const& v)
   : Data_(1, s)
{
   Data_.push_back(t);
   Data_.push_back(u);
   Data_.push_back(v);
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s)
   : Data_(1, CoerceSL(sl,s))
{
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t)
   : Data_(1, CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(1, CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, 
                 LatticeSite const& u, LatticeSite const& v)
   : Data_(1, CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
   Data_.push_back(CoerceSL(sl, v));
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(int RepeatCount, UnitCell const& l)
{
   Data_.reserve(l.size()*RepeatCount);
   for (int i = 0; i < RepeatCount; ++i)
   {
      Data_.insert(Data_.end(), l.begin(), l.end());
   }
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2)
   : Data_(x1.Data_)
{
   Data_.insert(Data_.end(), x2.begin(), x2.end());
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3)
   : Data_(x1.Data_)
{
   Data_.insert(Data_.end(), x2.begin(), x2.end());
   Data_.insert(Data_.end(), x3.begin(), x3.end());
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3, UnitCell const& x4)
   : Data_(x1.Data_)
{
   Data_.insert(Data_.end(), x2.begin(), x2.end());
   Data_.insert(Data_.end(), x3.begin(), x3.end());
   Data_.insert(Data_.end(), x4.begin(), x4.end());
   OperatorMap_["I"] = identity_mpo(*this);
}

UnitCell::UnitCell(int Size, LatticeSite const& s)
   : Data_(Size, s)
{
   OperatorMap_["I"] = identity_mpo(*this);
}

LatticeSite const& 
UnitCell::operator[](int n) const
{
   return Data_[n];
}

bool
UnitCell::operator_exists(std::string const& s) const
{
   if (OperatorMap_.find(s) != OperatorMap_.end())
      return true;
   if (Data_.size() == 1 && Data_.front().operator_exists(s))
      return true;
   return false;
}

FiniteMPO
UnitCell::Operator(std::string const& Op) const
{
   operator_map_type::const_iterator I = OperatorMap_.find(Op);
   if (I != OperatorMap_.end())
      return I->second;

   if (Data_.size() != 1)
   {
      PANIC("Operator not found in unit cell")(Op);
   }

   return this->Operator(Op, 0);
}

FiniteMPO&
UnitCell::Operator(std::string const& Op)
{
   return OperatorMap_[Op];
}

bool
UnitCell::operator_exists(std::string const& Op, int n) const
{
   if (n < 0 || n >= int(Data_.size()))
      return false;
   return Data_[n].operator_exists(Op);
}

LatticeCommute
UnitCell::Commute(std::string const& Op, int n) const
{
   CHECK(0 <= n && n < int(Data_.size()))("Site index is out of range")(n)(Data_.size());
   return Data_[n][Op].Commute();
}

FiniteMPO
UnitCell::Operator(std::string const& Op, int n) const
{
   CHECK(0 <= n && n < int(Data_.size()))("Site index is out of range")(n)(Data_.size());

   SiteOperator Operator = Data_[n][Op];
   std::string SignOperator = Operator.Commute().SignOperator();

   FiniteMPO Result(Data_.size(), Operator.Commute());

   BasisList Vacuum(Operator.GetSymmetryList());
   Vacuum.push_back(QuantumNumber(Operator.GetSymmetryList()));

   BasisList Basis(Operator.GetSymmetryList());
   Basis.push_back(Operator.TransformsAs());

   // Assemble the JW-string.
   for (int i = 0; i < n; ++i)
   {
      if (!Data_[i].operator_exists(SignOperator))
      {
	 WARNING("JW-string operator doesn't exist at a lattice site, using the identity")(i)(SignOperator);
      }
      SimpleOperator Op = Data_[i].operator_exists(SignOperator) ? Data_[i][SignOperator] : Data_[i]["I"];
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Basis, Basis);
      Result[i](0,0) = Op;
   }
   Result[n] = OperatorComponent(Operator.Basis1(), Operator.Basis2(), Basis, Vacuum);
   Result[n](0,0) = Operator;
   for (int i = n+1; i < this->size(); ++i)
   {
      SimpleOperator I = Data_[i]["I"];
      Result[i] = OperatorComponent(I.Basis1(), I.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = I;
   }

   return Result;
}

FiniteMPO
UnitCell::JWString(std::string const& Op, int n) const
{
   CHECK(0 <= n && n < int(Data_.size()))("Site index is out of range")(n)(Data_.size());

   SiteOperator Operator = Data_[n][Op];
   std::string SignOperator = Operator.Commute().SignOperator();

   // All JW strings are in the identity quantum number sector, so bosonic
   FiniteMPO Result(Data_.size(), LatticeCommute::Bosonic);

   BasisList Vacuum(Operator.GetSymmetryList());
   Vacuum.push_back(QuantumNumber(Operator.GetSymmetryList()));

   // Assemble the JW-string
   for (int i = 0; i < this->size(); ++i)
   {
      if (!Data_[i].operator_exists(SignOperator))
      {
	 WARNING("JW-string operator doesn't exist at a lattice site, using the identity")(i)(SignOperator);
      }
      SimpleOperator Op = Data_[i].operator_exists(SignOperator) ? Data_[i][SignOperator] : Data_[i]["I"];
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = Op;
   }

   return Result;
}


std::string
UnitCell::JWStringName(std::string const& Op, int n) const
{
   CHECK(0 <= n && n < int(Data_.size()))("Site index is out of range")(n)(Data_.size());

   SiteOperator Operator = Data_[n][Op];
   std::string SignOperator = Operator.Commute().SignOperator();
   return SignOperator;
}

PStream::opstream& operator<<(PStream::opstream& out, UnitCell const& L)
{
   out << L.Data_;
   out << L.OperatorMap_;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, UnitCell& L)
{
   in >>  L.Data_;
   in >>  L.OperatorMap_;
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

#if 0
FiniteMPO
UnitCell::Parse(std::string const& s)
{
}
#endif

FiniteMPO
identity_mpo(UnitCell const& c, QuantumNumbers::QuantumNumber const& q)
{
   FiniteMPO Result(c.size(), LatticeCommute::Bosonic);
   BasisList b(q.GetSymmetryList());
   b.push_back(q);
   for (int i = 0; i < c.size(); ++i)
   {
      Result[i] = OperatorComponent(c.LocalBasis(i), b, b);
      Result[i](0,0) = c[i]["I"];
   }
   return Result;
}

FiniteMPO
identity_mpo(UnitCell const& c)
{
   return identity_mpo(c, QuantumNumbers::QuantumNumber(c.GetSymmetryList()));
}

