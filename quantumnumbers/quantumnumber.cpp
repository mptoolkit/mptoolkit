// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/quantumnumber.cpp
//
// Copyright (C) 2004-2020 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "common/trace.h"
#include "common/poolallocator.h"
#include "quantumnumber.h"
#include "symmetrylist.h"
#include <algorithm>

namespace QuantumNumbers
{

//
// QuantumNumber
//

QuantumNumber::QuantumNumber(SymmetryList const& qL, std::string const& Str)
  : RepLabelBase<QuantumNumber>(qL.GetSymmetryListImpl(), qL.QuantumNumberSize())
{
   QN_TRACE("Constructing from string ")(Str);
   std::string::const_iterator beg = Str.begin(), end = Str.end();

   for (int i = 0; i < qL.NumSymmetries(); ++i)
   {
      while (*beg == ' ')
         ++beg;
      size_t Offset = qL.QuantumNumberOffset(i);

      std::string::const_iterator next = std::find(beg, end, ',');
      qL.GetSymmetryBase(i)->FromString(std::string(beg, next), this->begin()+Offset);
      beg = next;
      if (beg != end) ++beg;
   }
   CHECK(beg == end);
}

QuantumNumber::QuantumNumber(SymmetryList const& qL, char const* s)
  : RepLabelBase<QuantumNumber>(qL.GetSymmetryListImpl(), qL.QuantumNumberSize())
{
   QN_TRACE("Constructing from char* ")(s);
   char const* beg = s;
   char const* end = s+strlen(s);

   for (int i = 0; i < qL.NumSymmetries(); ++i)
   {
      while (*beg == ' ')
         ++beg;
      size_t Offset = qL.QuantumNumberOffset(i);

      char const* next = std::find(beg, end, ',');
      qL.GetSymmetryBase(i)->FromString(std::string(beg, next), this->begin()+Offset);
      beg = next;
      if (beg != end) ++beg;
   }
   CHECK(beg == end);
}

QuantumNumber::QuantumNumber(SymmetryList const& qL, char* s)
  : RepLabelBase<QuantumNumber>(qL.GetSymmetryListImpl(), qL.QuantumNumberSize())
{
   QN_TRACE("Constructing from char* ")(s);
   char const* beg = s;
   char const* end = s+strlen(s);

   for (int i = 0; i < qL.NumSymmetries(); ++i)
   {
      while (*beg == ' ')
         ++beg;
      size_t Offset = qL.QuantumNumberOffset(i);

      char const* next = std::find(beg, end, ',');
      qL.GetSymmetryBase(i)->FromString(std::string(beg, next), this->begin()+Offset);
      beg = next;
      if (beg != end) ++beg;
   }
   CHECK(beg == end);
}

std::ostream& operator<<(std::ostream& out, QuantumNumber const& Q)
{
   return out << Q.ToString();
}

std::string QuantumNumber::ToString() const
{
   if (this->is_null()) return "<null>";
   SymmetryListImpl const* MyImpl = this->GetSymmetryListImpl();
   std::string Result;
   for (int i = 0; i < MyImpl->NumSymmetries(); ++i)
   {
      size_t Offset = MyImpl->QuantumNumberOffset(i);

      if (i > 0) Result += ',';
      Result += MyImpl->GetSymmetryBase(i)->ToString(this->begin()+Offset);
   }
   return Result;
}

int
QuantumNumber::degree() const
{
   PRECONDITION(!this->is_null());
   SymmetryListImpl const* MyImpl = this->GetSymmetryListImpl();
   int Result = 1;
   for (int i = 0; i < MyImpl->NumSymmetries(); ++i)
   {
      size_t Offset = MyImpl->QuantumNumberOffset(i);

      Result *= MyImpl->GetSymmetryBase(i)->degree(this->begin()+Offset);
   }
   return Result;
}

QuantumNumber CoerceSymmetryList(QuantumNumber const& q, SymmetryList const& SList)
{
   if (SList == q.GetSymmetryList()) return q;

   QuantumNumber Other(SList);
   int* Iter = Other.begin();
   int ThisCountCheck = 0; // count the number of symmetries from *this that we've inserted

   for (int i = 0; i < SList.NumSymmetries(); ++i)
   {
      SymmetryBase const* S = SList.GetSymmetryBase(i);
      size_t Size = S->QuantumNumberSize();
      // see if S is located in this symmetry list
      int ThisWhere = q.GetSymmetryList().WhichSymmetry(SList.SymmetryName(i));

      QN_TRACE(ThisWhere)(S)(Size);

      if (ThisWhere == -1)
      {
         // this symmetry is new, initialize it with the identity quantum number
         S->scalar_transforms_as(Iter);
      }
      else
      {
         // copy the quantum number from *this into Other
         CHECK_EQUAL(q.GetSymmetryList().GetSymmetryBase(ThisWhere), S);
         size_t ThisOffset = q.GetSymmetryList().QuantumNumberOffset(ThisWhere);

         QN_TRACE(ThisOffset);

         memcpy(Iter, q.begin()+ThisOffset, Size * sizeof(*Iter));
         ++ThisCountCheck;
      }
      Iter += Size;
   }

   // Remove the limitation that the symmetry must be included - this lets us remove
   // symmetries from the symmetry list
   // NOTE: we should do a check here that any symmetry we've removed is abelian.
   // Removing a non-abelian symmetry will do very weird things!  Or more precisely,
   // the quantum number should have degree 1 in the component that we removed.
#if 0
   // make sure that we've included all symmetries
   CHECK(ThisCountCheck == q.GetSymmetryList().NumSymmetries())
      ("Quantum number cannot be coerced into the symmetry list")
      (q)(q.GetSymmetryList())(SList);
#endif

   return Other;
}

PStream::opstream& operator<<(PStream::opstream& out, QuantumNumber const& L)
{
   if (L.is_null())
   {
      out << std::string("");
   }
   else
   {
      out << L.GetSymmetryList();
      std::copy(L.begin(), L.end(), PStream::opstream_iterator<QuantumNumber::value_type>(out));
   }
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, QuantumNumber& L)
{
   // We cannot combine the stream extractions and do
   // L = QuantumNumber(in.read<SymmetryList>(), PStream::ipstream_iterator<QuantumNumber::value_type>(in))
   // because the order of evaluation of function arguments is not specified.
   SymmetryList S;
   in >> S;
   if (S.is_null())
      L = QuantumNumber();
   else
      L = QuantumNumber(S, PStream::ipstream_iterator<QuantumNumber::value_type>(in));
   return in;
}

//
// Projection
//

Projection::Projection() noexcept
{
}

Projection::Projection(SymmetryList const& qL)
  : RepLabelBase<Projection>(qL.GetSymmetryListImpl(), qL.ProjectionSize())
{
}

Projection::Projection(SymmetryList const& qL, std::string const& Str)
  : RepLabelBase<Projection>(qL.GetSymmetryListImpl(), qL.ProjectionSize())
{
   std::string::const_iterator beg = Str.begin(), end = Str.end();

   for (int i = 0; i < qL.NumSymmetries(); ++i)
   {
      size_t Offset = qL.ProjectionOffset(i);

      std::string::const_iterator next = std::find(beg, end, ',');
      qL.GetSymmetryBase(i)->ProjectionFromString(std::string(beg, next), begin()+Offset);
      beg = next;
      if (beg != end) ++beg;
   }
   CHECK(beg == end);
}

Projection::Projection(SymmetryList const& qL, char const* s)
  : RepLabelBase<Projection>(qL.GetSymmetryListImpl(), qL.ProjectionSize())
{
   char const* beg = s;
   char const* end = s+strlen(s);

   for (int i = 0; i < qL.NumSymmetries(); ++i)
   {
      size_t Offset = qL.ProjectionOffset(i);

      char const* next = std::find(beg, end, ',');
      qL.GetSymmetryBase(i)->ProjectionFromString(std::string(beg, next), begin()+Offset);
      beg = next;
      if (beg != end) ++beg;
   }
   CHECK(beg == end);
}

Projection::Projection(SymmetryList const& qL, char* s)
  : RepLabelBase<Projection>(qL.GetSymmetryListImpl(), qL.ProjectionSize())
{
   char const* beg = s;
   char const* end = s+strlen(s);

   for (int i = 0; i < qL.NumSymmetries(); ++i)
   {
      size_t Offset = qL.ProjectionOffset(i);

      char const* next = std::find(beg, end, ',');
      qL.GetSymmetryBase(i)->ProjectionFromString(std::string(beg, next), begin()+Offset);
      beg = next;
      if (beg != end) ++beg;
   }
   CHECK(beg == end);
}

std::string Projection::ToString() const
{
   PRECONDITION(!this->is_null());
   SymmetryListImpl const* MyImpl = this->GetSymmetryListImpl();
   std::string Result;
   for (int i = 0; i < MyImpl->NumSymmetries(); ++i)
   {
      size_t Offset = MyImpl->ProjectionOffset(i);

      if (i > 0) Result += ',';
      Result += MyImpl->GetSymmetryBase(i)->ProjectionToString(this->begin()+Offset);
   }
   return Result;
}

Projection CoerceSymmetryList(Projection const& q, SymmetryList const& SList)
{
   if (SList == q.GetSymmetryList()) return q;

   Projection Other(SList);
   int* Iter = Other.begin();
   int ThisCountCheck = 0; // count the number of symmetries from *this that we've inserted

   for (int i = 0; i < SList.NumSymmetries(); ++i)
   {
      SymmetryBase const* S = SList.GetSymmetryBase(i);
      size_t Size = S->ProjectionSize();
      // see if S is located in this symmetry list
      int ThisWhere = q.GetSymmetryList().WhichSymmetry(SList.SymmetryName(i));

      QN_TRACE(ThisWhere)(S)(Size);

      if (ThisWhere == -1)
      {
         // this symmetry is new, initialize it with the identity quantum number
         S->scalar_transforms_as(Iter);
      }
      else
      {
         // copy the quantum number from *this into Other
         CHECK(q.GetSymmetryList().GetSymmetryBase(ThisWhere) == S);
         size_t ThisOffset = q.GetSymmetryList().ProjectionOffset(ThisWhere);

         QN_TRACE(ThisOffset);

         memcpy(Iter, q.begin()+ThisOffset, Size * sizeof(*Iter));
         ++ThisCountCheck;
      }
      Iter += Size;
   }
   // make sure that we've included all symmetries
   CHECK(ThisCountCheck == q.GetSymmetryList().NumSymmetries())
        ("Quantum number cannot be coerced into the symmetry list")
        (q)(q.GetSymmetryList())(SList);

   return Other;
}

void Projection::CoerceSymmetryList(SymmetryList const& SList)
{
   *this = QuantumNumbers::CoerceSymmetryList(*this, SList);
}

std::ostream& operator<<(std::ostream& out, Projection const& P)
{
   return out << P.ToString();
}

PStream::opstream& operator<<(PStream::opstream& out, Projection const& L)
{
   out << L.GetSymmetryList();
   std::copy(L.begin(), L.end(), PStream::opstream_iterator<Projection::value_type>(out));
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, Projection& L)
{
   SymmetryList S;
   in >> S;
   L = Projection(S, PStream::ipstream_iterator<Projection::value_type>(in));
   return in;
}

//
// QuantumNumberList
//

PStream::opstream& operator<<(PStream::opstream& out, QuantumNumberList const& q)
{
   int Size = q.size();
   out << Size;
   if (Size > 0)
   {
      out << q[0].GetSymmetryList();
      for (int i = 0; i < Size; ++i)
      {
         q[i].WriteRaw(out);
      }
   }
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, QuantumNumberList& q)
{
   q.clear();
   int Size;
   in >> Size;
   if (Size > 0)
   {
      SymmetryList SList;
      in >> SList;
      QuantumNumber QN(SList, QuantumNumber::NoInitialization());
      for (int i = 0; i < Size; ++i)
      {
         QN.ReadRaw(in);
         q.push_back(QN);
      }
   }
   return in;
}

// QuantumNumberList

void
QuantumNumberList::delta_shift(QuantumNumber const& q)
{
   for (iterator I = this->begin(); I != this->end(); ++I)
   {
      *I = ::QuantumNumbers::delta_shift(*I, q);
   }
}

} // namespace QuantumNumber
