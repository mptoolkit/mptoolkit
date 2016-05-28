// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/lattice.cpp
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

#include "lattice.h"

SiteBlock flip_conj(SiteBlock const& A)
{
   SiteBlock Result;

   if (A.empty())
      return Result;

   SiteBasis ReflectedBasis = adjoint(A.Basis1());
   for (SiteBlock::const_iterator ai = A.begin(); ai != A.end(); ++ai)
   {
      Result[ai->first] = flip_conj(ai->second, ReflectedBasis);
   }
   return Result;
}

//
// Lattice members
//

Lattice::Lattice()
   : CoordinatesFixed_(false)
{
}

Lattice::Lattice(SiteBlock const& s)
   : Data_(s), Coordinates_(1), CoordinatesFixed_(false)
{
}

Lattice::Lattice(SiteBlock const& s, SiteBlock const& t)
   : Data_(s), Coordinates_(2), CoordinatesFixed_(false)
{
   Data_.push_back(t);
}

Lattice::Lattice(SiteBlock const& s, SiteBlock const& t, SiteBlock const& u)
   : Data_(s), Coordinates_(3), CoordinatesFixed_(false)
{
   Data_.push_back(t);
   Data_.push_back(u);
}

Lattice::Lattice(SiteBlock const& s, SiteBlock const& t, SiteBlock const& u, SiteBlock const& v)
   : Data_(s), Coordinates_(4), CoordinatesFixed_(false)
{
   Data_.push_back(t);
   Data_.push_back(u);
   Data_.push_back(v);
}

SiteBlock CoerceSL(SymmetryList const& sl, SiteBlock const& s)
{
   SiteBlock r(s);
   r.CoerceSymmetryList(sl);
   return r;
}

Lattice::Lattice(SymmetryList const& sl, SiteBlock const& s)
   : Data_(CoerceSL(sl,s)), Coordinates_(1), CoordinatesFixed_(false)
{
}

Lattice::Lattice(SymmetryList const& sl, SiteBlock const& s, SiteBlock const& t)
   : Data_(CoerceSL(sl, s)), Coordinates_(2), CoordinatesFixed_(false)
{
   Data_.push_back(CoerceSL(sl, t));
}

Lattice::Lattice(SymmetryList const& sl, SiteBlock const& s, SiteBlock const& t, SiteBlock const& u)
   : Data_(CoerceSL(sl, s)), Coordinates_(3), CoordinatesFixed_(false)
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
}

Lattice::Lattice(SymmetryList const& sl, SiteBlock const& s, SiteBlock const& t, 
                 SiteBlock const& u, SiteBlock const& v)
   : Data_(CoerceSL(sl, s)), Coordinates_(4), CoordinatesFixed_(false)
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
   Data_.push_back(CoerceSL(sl, v));
}

Lattice::Lattice(int RepeatCount, Lattice const& l)
   : Data_(RepeatCount, l.data()), 
     CoordinatesFixed_(false)
{
   Coordinates_.reserve(l.Coordinates_.size() * RepeatCount);
   for (int i = 0; i < RepeatCount; ++i)
      Coordinates_.insert(Coordinates_.end(), l.Coordinates_.begin(), l.Coordinates_.end());
}

Lattice::Lattice(Lattice const& x1, Lattice const& x2)
   : Data_(x1.Data_),
     Coordinates_(x1.Coordinates_), 
     CoordinatesFixed_(false)
{
   Data_.push_back(x2.Data_);
   Coordinates_.insert(Coordinates_.end(), x2.Coordinates_.begin(), x2.Coordinates_.end());
}

Lattice::Lattice(Lattice const& x1, Lattice const& x2, Lattice const& x3)
   : Data_(x1.Data_),
     Coordinates_(x1.Coordinates_),
     CoordinatesFixed_(false)
{
   Data_.push_back(x2.Data_);
   Data_.push_back(x3.Data_);
   Coordinates_.insert(Coordinates_.end(), x2.Coordinates_.begin(), x2.Coordinates_.end());
   Coordinates_.insert(Coordinates_.end(), x3.Coordinates_.begin(), x3.Coordinates_.end());
}

Lattice::Lattice(Lattice const& x1, Lattice const& x2, Lattice const& x3, Lattice const& x4)
   : Data_(x1.Data_),
     Coordinates_(x1.Coordinates_),
     CoordinatesFixed_(false)
{
   Data_.push_back(x2.Data_);
   Data_.push_back(x3.Data_);
   Data_.push_back(x4.Data_);
   Coordinates_.insert(Coordinates_.end(), x2.Coordinates_.begin(), x2.Coordinates_.end());
   Coordinates_.insert(Coordinates_.end(), x3.Coordinates_.begin(), x3.Coordinates_.end());
   Coordinates_.insert(Coordinates_.end(), x4.Coordinates_.begin(), x4.Coordinates_.end());
}

Lattice::Lattice(int Size, SiteBlock const& s)
   : Data_(Size, run_length_compressed<SiteBlock>(s)),
     Coordinates_(Size),
     CoordinatesFixed_(false)
{
}

Lattice::Lattice(SiteBlock const& s, std::string const& Coord)
   : Data_(s),
     Coordinates_(1, Coord),
     CoordinatesFixed_(false)
{
}

std::string const& 
Lattice::coordinate_at_site(int Site) const
{
   CHECK(CoordinatesFixed_);
   return Coordinates_[Site-1];
}

SiteBlock const& 
Lattice::operator[](int n) const
{
   const_iterator I = this->begin();
   std::advance(I, n-1);
   return *I;
}

int 
Lattice::site_at_coordinate(std::string const& s) const
{
   CHECK_EQUAL(CoordinatesFixed_, true)("The lattice coordinates are not fixed!");
   std::map<std::string, int>::const_iterator I = SiteAtCoord_.find(s);
   if (I == SiteAtCoord_.end())
      return -1;
   return I->second;
}

void 
Lattice::fix_coordinates()
{
   //   CHECK(!CoordinatesFixed_);
   std::map<std::string, int> UseCount;
   std::vector<std::string>::iterator IEnd = Coordinates_.end();
   for (std::vector<std::string>::iterator I = Coordinates_.begin(); I != IEnd; ++I)
   {
      int c = ++UseCount[*I];
      if (!I->empty())
         *I += ',';
      *I += boost::lexical_cast<std::string>(c);
   }
   this->Fixate();
}

void 
Lattice::fix_coordinates_starting_from(int n)
{
   //   CHECK(!CoordinatesFixed_);
   std::map<std::string, int> UseCount;
   std::vector<std::string>::iterator IEnd = Coordinates_.end();
   for (std::vector<std::string>::iterator I = Coordinates_.begin(); I != IEnd; ++I)
   {
      int c = (UseCount[*I]++) + n;
      if (!I->empty())
         *I += ',';
      *I += boost::lexical_cast<std::string>(c);
   }
   this->Fixate();
}

void 
Lattice::fix_coordinates_reverse()
{
   //   CHECK(!CoordinatesFixed_);
   std::map<std::string, int> UseCount;
   std::vector<std::string>::iterator IBegin = Coordinates_.begin();
   std::vector<std::string>::iterator I = Coordinates_.end();
   while (I != IBegin)
   {
      --I;
      int c = ++UseCount[*I];
      if (!I->empty())
         *I += ',';
      *I += boost::lexical_cast<std::string>(c);
   }
   this->Fixate();
}

void 
Lattice::fix_coordinates_prepend()
{
   //   CHECK(!CoordinatesFixed_);
   std::map<std::string, int> UseCount;
   std::vector<std::string>::iterator IEnd = Coordinates_.end();
   for (std::vector<std::string>::iterator I = Coordinates_.begin(); I != IEnd; ++I)
   {
      int c = ++UseCount[*I];
      if (I->empty())
         *I = boost::lexical_cast<std::string>(c);
      else
         *I = boost::lexical_cast<std::string>(c) + ',' + (*I);
   }
   this->Fixate();
}

void 
Lattice::fix_coordinates_prepend_starting_from(int n)
{
   //   CHECK(!CoordinatesFixed_);
   std::map<std::string, int> UseCount;
   std::vector<std::string>::iterator IEnd = Coordinates_.end();
   for (std::vector<std::string>::iterator I = Coordinates_.begin(); I != IEnd; ++I)
   {
      int c = (UseCount[*I]++) + n;
      if (I->empty())
         *I = boost::lexical_cast<std::string>(c);
      else
         *I = boost::lexical_cast<std::string>(c) + ',' + (*I);
   }
   this->Fixate();
}

void Lattice::fix_coordinates_unique()
{
   //   CHECK(!CoordinatesFixed_);
   std::set<std::string> Used(Coordinates_.begin(), Coordinates_.end());
   CHECK_EQUAL(Used.size(), Coordinates_.size())("Coordinates are not distinct!");
   this->Fixate();
}

void Lattice::Fixate()
{
   int Site = 1;
   SiteAtCoord_.clear();
   std::vector<std::string>::const_iterator IEnd = Coordinates_.end();
   for (std::vector<std::string>::const_iterator I = Coordinates_.begin(); I != IEnd; ++I)
      SiteAtCoord_[*I] = Site++;
   CoordinatesFixed_ = true;
}

void Lattice::fix_coordinates(char const* c1)
{
   CHECK_EQUAL(Data_.size(), 1);
   Coordinates_.clear();
   Coordinates_.push_back(std::string(c1));
   this->Fixate();
}

void Lattice::fix_coordinates(char const* c1, char const* c2)
{
   CHECK_EQUAL(Data_.size(), 2);
   Coordinates_.clear();
   Coordinates_.push_back(std::string(c1));
   Coordinates_.push_back(std::string(c2));
   this->Fixate();
}

void Lattice::fix_coordinates(char const* c1, char const* c2, char const* c3)
{
   CHECK_EQUAL(Data_.size(), 3);
   Coordinates_.clear();
   Coordinates_.push_back(std::string(c1));
   Coordinates_.push_back(std::string(c2));
   Coordinates_.push_back(std::string(c3));
   this->Fixate();
}

void Lattice::fix_coordinates(char const* c1, char const* c2, char const* c3, char const* c4)
{
   CHECK_EQUAL(Data_.size(), 4);
   Coordinates_.clear();
   Coordinates_.push_back(std::string(c1));
   Coordinates_.push_back(std::string(c2));
   Coordinates_.push_back(std::string(c3));
   Coordinates_.push_back(std::string(c4));
   this->Fixate();
}

PStream::opstream& operator<<(PStream::opstream& out, Lattice const& L)
{
   return out << L.Data_ << L.Coordinates_ << L.CoordinatesFixed_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, Lattice& L)
{
   in >>  L.Data_ >> L.Coordinates_ >> L.CoordinatesFixed_;
   if (L.CoordinatesFixed_)
      L.Fixate();
   return in;
}

Lattice repeat(Lattice const& x, int RepeatCount)
{
   return Lattice(RepeatCount, x);
}

Lattice join(Lattice const& x, Lattice const& y)
{
   return Lattice(x, y);
}

Lattice join(Lattice const& x, Lattice const& y, Lattice const& z)
{
   return Lattice(x,y,z);
}

Lattice join(Lattice const& x, Lattice const& y, Lattice const& z, Lattice const& w)
{
   return Lattice(x,y,z,w);
}

Lattice join(Lattice const& x, Lattice const& y, Lattice const& z, Lattice const& w,
	     Lattice const& v)
{
   return join(Lattice(x,y,z,w), v);
}

//
// CreateTensorProductOperator
//

struct DoCreateTensorProductOperator : boost::static_visitor<MPOpCompressed>
{
   DoCreateTensorProductOperator(std::string const& Op_, BasisList const& b_)
      : Op(Op_), b(b_) { DEBUG_CHECK_EQUAL(b.size(), 1)(b); }

   MPOpCompressed operator()(SiteBlock const& s) const
   {
      SiteBlock::const_iterator I = s.find(Op);
      if (I == s.end())     // fallback case; if the operator doesn't exist at this site use I
         I = s.find("I");
      CHECK(I != s.end())("site operator is missing both the required operator and the identity!")(Op);
      MPOpComponent Result(I->second.Basis().Basis(), b, b);
      SimpleOperator Element(b, b, QuantumNumber(b.GetSymmetryList()));
      Element(0,0) = 1.0;
      Result[QuantumNumber(b.GetSymmetryList())] = scalar(Element) * I->second.base();
      return MPOpCompressed(Result);
   }

   MPOpCompressed operator()(run_length_array<SiteBlock> const& s) const
   {
      MPOpArray Result;
      std::transform(s.begin(), s.end(), std::back_inserter(Result),
                     boost::apply_visitor(*this));
      return Result;
   }

   MPOpCompressed operator()(run_length_repeat<SiteBlock> const& s) const
   {
      return MPOpRepeat(s.size(), s.nested().apply_visitor(*this));
   }

   std::string Op;
   BasisList b;
};

MPOpCompressed 
CreateTensorProductOpCompressed(Lattice::data_type const& L, 
                                std::string const& Operator, 
                                BasisList const& b)
{
   return L.apply_visitor(DoCreateTensorProductOperator(Operator, b));
}

//
// CreateMPOpCompressed
//

MPOpCompressed CreateMPOpCompressed(Lattice const& L, std::string const& Operator, int Site)
{
   Lattice::data_type Left, Right;
   SiteBlock s;
   std::tie(Left, s, Right) = split_lcr(L.data(), Site);
   SiteOperator SiteOp = s[Operator];
   CHECK(SiteOp.Commute() != LatticeCommute::None);

   BasisList LeftBasis(SiteOp.GetSymmetryList());
   LeftBasis.push_back(SiteOp.TransformsAs());

   // vacuum basis
   BasisList RightBasis(SiteOp.GetSymmetryList());
   RightBasis.push_back(QuantumNumber(SiteOp.GetSymmetryList()));

   std::string CommuteOperator = SiteOp.Commute().SignOperator();

   MPOpCompressed LeftOp = CreateTensorProductOpCompressed(Left, CommuteOperator, LeftBasis);
   MPOpCompressed RightOp = CreateTensorProductOpCompressed(Right, "I", RightBasis);

   // The middle operator
   MPOpComponent MiddleOp(SiteOp.Basis().Basis(), LeftBasis, RightBasis);
   SimpleOperator Element(LeftBasis, RightBasis, SiteOp.TransformsAs());
   Element(0,0) = 1;
   MiddleOp[SiteOp.TransformsAs()] = scalar(Element) * SiteOp.base();

   return join(LeftOp, join(MPOpCompressed(MiddleOp), RightOp));
}

//
// CreateMPOperator
//

MPOperator CreateMPOperator(Lattice const& L, std::string const& Operator, int Site)
{
   return MPOperator(CreateMPOpCompressed(L, Operator, Site));
}


QuantumNumbers::QuantumNumber
LatticeSiteOperatorTransformsAs(Lattice const& L, std::string const& Operator)
{
   // dumb implementation
   for (Lattice::const_iterator I = L.begin(); I != L.end(); ++I)
   {
      SiteBlock::const_iterator Si = I->find(Operator);
      if (Si != I->end())
         return Si->second.TransformsAs();
   }
   // fallback: the operator does not exist at all, return identity quantum number
   return QuantumNumber(L.front().GetSymmetryList());
}

MPOpComponent ConstructFromSymbolic(LinearAlgebra::SparseMatrix<std::string> const& M,
                                    SiteBlock const& Block,
                                    BasisList const& B)
{
   typedef LinearAlgebra::SparseMatrix<std::string> MatType;
   MPOpComponent Res(Block.Basis1().Basis(), B, B);
   for (const_iterator<MatType>::type I = iterate(M); I; ++I)
   {
      for (const_inner_iterator<MatType>::type J = iterate(I); J; ++J)
      {
         SiteBlock::const_iterator Si = Block.find(*J);
         if (Si != Block.end())
         {
            SimpleOperator E(B, B, Si->second.TransformsAs());
            E(J.index1(), J.index2()) = 1.0;
            Res[Si->second.TransformsAs()] += scalar(E) * Si->second.base();
         }
      }
   }
   return Res;
}

// this is such a common pattern, it should be some kind of transform() functor
struct DoCreateRepeatedMPOpCompressed : boost::static_visitor<MPOpCompressed>
{
   typedef LinearAlgebra::SparseMatrix<std::string> MatType;
   DoCreateRepeatedMPOpCompressed(BasisList const& B_, MatType const& M_)
      : B(B_), M(M_) {}

   MPOpCompressed operator()(SiteBlock const& s) const
   {
      return MPOpCompressed(ConstructFromSymbolic(M, s, B));
   }

   MPOpCompressed operator()(run_length_array<SiteBlock> const& s) const
   {
      MPOpArray Result;
      std::transform(s.begin(), s.end(), std::back_inserter(Result),
                     boost::apply_visitor(*this));
      return Result;
   }

   MPOpCompressed operator()(run_length_repeat<SiteBlock> const& s) const
   {
      return MPOpRepeat(s.size(), s.nested().apply_visitor(*this));
   }

   BasisList B;
   MatType M;
};

// utility

MPOperator
CreateRepeatedOperator(Lattice const& L, std::string const& Op1)
{
   // Firstly, construct the basis.
   BasisList B(L.front().GetSymmetryList());
   B.push_back(QuantumNumber(B.GetSymmetryList()));
   B.push_back(QuantumNumber(B.GetSymmetryList()));
   typedef LinearAlgebra::SparseMatrix<std::string> MatType;
   MatType M(2,2);
   M(0,0) = "I";
   M(1,1) = "I";
   M(1,0) = Op1;
   return MPOperator(
     project_component(L.apply_visitor(DoCreateRepeatedMPOpCompressed(B, M)), 1, 0));
}

MPOperator
CreateRepeatedOperator(Lattice const& L, std::string const& Op1, std::string const& Op2)
{
   // Firstly, construct the basis.
   BasisList B(L.front().GetSymmetryList());
   B.push_back(QuantumNumber(B.GetSymmetryList()));
   QuantumNumbers::QuantumNumber q2 = LatticeSiteOperatorTransformsAs(L, Op2);
   CHECK_EQUAL(adjoint(LatticeSiteOperatorTransformsAs(L, Op1)), q2)(Op1)(Op2);
   B.push_back(q2);
   B.push_back(QuantumNumber(B.GetSymmetryList()));
   typedef LinearAlgebra::SparseMatrix<std::string> MatType;
   MatType M(3,3);
   M(0,0) = "I";
   M(2,2) = "I";
   M(1,0) = Op2;
   M(2,1) = Op1;
   MPOpCompressed Temp = L.apply_visitor(DoCreateRepeatedMPOpCompressed(B, M));
   // sqrt(degree(q2)) here for the dot product scale factor
   return std::sqrt(double(degree(q2)))*MPOperator(project_component(Temp, 2, 0));
}

MPOperator
CreateRepeatedOperator(Lattice const& L, std::string const& Op1, std::string const& Op2,
                       std::string const& Op3)
{
   // Firstly, construct the basis.
   BasisList B(L.front().GetSymmetryList());
   B.push_back(QuantumNumber(B.GetSymmetryList()));
   QuantumNumbers::QuantumNumber q2 = LatticeSiteOperatorTransformsAs(L, Op2);
   QuantumNumbers::QuantumNumber q3 = LatticeSiteOperatorTransformsAs(L, Op3);
   CHECK_EQUAL(q2, QuantumNumber(B.GetSymmetryList()))("Middle operator must be a scalar");
   CHECK_EQUAL(adjoint(LatticeSiteOperatorTransformsAs(L, Op1)), q3);
   B.push_back(q3);
   B.push_back(q3);
   B.push_back(QuantumNumber(B.GetSymmetryList()));
   typedef LinearAlgebra::SparseMatrix<std::string> MatType;
   MatType M(4,4);
   M(0,0) = "I";
   M(3,3) = "I";
   M(1,0) = Op3;
   M(2,1) = Op2;
   M(3,2) = Op1;
   MPOpCompressed Temp = L.apply_visitor(DoCreateRepeatedMPOpCompressed(B, M));
   // sqrt(degree(q3)) here for the dot product scale factor
   return std::sqrt(double(degree(q3)))*MPOperator(project_component(Temp, 3, 0));
}
