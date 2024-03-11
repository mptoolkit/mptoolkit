// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/symmetrylist.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
#include "symmetrylist.h"
#include "quantumnumber.h"
#include <algorithm>
#include "common/stringutil.h" // for Split function

namespace QuantumNumbers
{

namespace
{

// helper function to return a SymmetryListImpl* from a given list of symmetry types
SymmetryListImpl* CreateSymmetryList(std::string const& List)
{
   // search the list of created quantum number lists for one that
   // matches the name
   SymmetryListImpl* pImpl = SymmetryListImpl::SearchForCreated(List);
   if (pImpl) return pImpl;

   // else create a new SymmetryListImpl with our quantum numbers
   QN_TRACE("Creating new SymmetryListImpl()")(List);

   pImpl = new SymmetryListImpl(List);
   std::vector<std::string> NameType;
   Split(List, ',', std::back_inserter(NameType));
   for (size_t i = 0U; i < NameType.size(); ++i)
   {
      std::string::iterator ColonLoc = std::find(NameType[i].begin(), NameType[i].end(), ':');
      CHECK(ColonLoc != NameType[i].end())("Colon missing in symmetry list")(List);
      pImpl->Append(std::string(NameType[i].begin(), ColonLoc),
                    SymmetryBase::Create(std::string(ColonLoc+1, NameType[i].end())));
   }

   // Check that the reconstructed full name agrees with pImpl->FullName()
   std::string CheckFullName;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      if (i != 0) CheckFullName += ',';
      CheckFullName += pImpl->SymmetryName(i) + ':' + pImpl->SymmetryType(i);
   }
   if (CheckFullName != pImpl->FullName())
   {
     PANIC("") << "Unexpected mismatch between supplied symmetry list name " << pImpl->FullName()
            << " and the reconstructed name " << CheckFullName;
   }

   return pImpl;
}

} // namespace

//
// SymmetryListImpl
//

// declared as a pointer so that it is statically initialized.  The
// real initialization is done in InitializeInstances via a nifty counter
SymmetryListImpl::InstanceListType* SymmetryListImpl::Instances = NULL;

void
SymmetryListImpl::InitializeInstances()
{
   Instances = new InstanceListType();
}

SymmetryListImpl::SymmetryListImpl(std::string const& FName_)
  : FName(FName_), count(0), projectionCount(0)
{
   QN_TRACE("SymmetryListImpl::SymmetryListImpl()")(this);

   // add this instance to the static list
   Instances->push_front(this);
   MyInstance = Instances->begin();
}

SymmetryListImpl::~SymmetryListImpl()
{
   QN_TRACE("SymmetryListImpl::~SymmetryListImpl()")(this);

   for (size_t i = 0U; i < Data.size(); ++i)
   {
      delete Data[i];
   }
   // remove this from the static list
   Instances->erase(MyInstance);
}

void
SymmetryListImpl::Append(std::string const& Name, SymmetryBase const* N)
{
   PRECONDITION(N != NULL);
   Data.push_back(N);
   Names.push_back(Name);
   Offsets.push_back(count);
   count += N->QuantumNumberSize();
   ProjectionOffsets.push_back(projectionCount);
   projectionCount += N->ProjectionSize();
}

SymmetryListImpl* SymmetryListImpl::SearchForCreated(std::string const& Name)
{
   for (InstanceListType::const_iterator I = Instances->begin(); I != Instances->end(); ++I)
   {
      if ((*I)->FullName() == Name) return *I;
   }
   return NULL;
}

int
SymmetryListImpl::NumCasimirOperators() const
{
   int Result = 0;
   for (unsigned i = 0; i < Data.size(); ++i)
      Result += Data[i]->num_casimir();
   return Result;
}

std::string
SymmetryListImpl::CasimirName(int n) const
{
   // find out which quantum number owns the n'th casimir
   int i = 0;
   while (n >= Data[i]->num_casimir())
      n -= Data[i++]->num_casimir();
   return Data[i]->casimir_name(Names[i], n);
}

//
// SymmetryList
//

std::ostream& operator<<(std::ostream& out, SymmetryList const& s)
{
   if (s.is_null())
      out << "(null SymmetryList)";
   else
      out << s.FullName();
#if 0
   out << ", QuantumNumberSize=" << s.QuantumNumberSize() << ", NumSymmetries="
       << s.NumSymmetries() << ", ";
   for (int i = 0; i < s.NumSymmetries(); ++i)
   {
      if (i != 0) out << ", ";
      out << '(' << s.SymmetryName(i) << ':' << s.SymmetryType(i) << '='
          << s.QuantumNumberOffset(i) << '+' << s.QuantumNumberSize(i) << ')';
   }
#endif
   return out;
}

PStream::opstream& operator<<(PStream::opstream& out, SymmetryList const& L)
{
   return out << L.FullName();
}

PStream::ipstream& operator>>(PStream::ipstream& in, SymmetryList& L)
{
   std::string S = in.read<std::string>();
   L = SymmetryList(S);
   return in;
}

SymmetryList::SymmetryList()
   //  : pImpl(CreateSymmetryList(""))
   : pImpl(NULL)
{
}

SymmetryList::SymmetryList(std::string const& List)
  : pImpl(CreateSymmetryList(List))
{
}

SymmetryList::SymmetryList(char const* List)
  : pImpl(CreateSymmetryList(List))
{
}

int SymmetryList::degree(sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   int d = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      d *= pImpl->GetSymmetryBase(i)->degree(q+pImpl->QuantumNumberOffset(i));
   }
   return d;
}

double SymmetryList::trace(sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   double d = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      d *= pImpl->GetSymmetryBase(i)->trace(q+pImpl->QuantumNumberOffset(i));
   }
   return d;
}

bool
SymmetryList::cross_product_exists(sc_iter q1, sc_iter q2) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      if (!(pImpl->GetSymmetryBase(i)->cross_product_exists(q1+pImpl->QuantumNumberOffset(i),
                                                            q2+pImpl->QuantumNumberOffset(i))))
         return false;
   }
   return true;
}

void
SymmetryList::cross_product_transforms_as(sc_iter q1, sc_iter q2, s_iter qOut) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);
      pImpl->GetSymmetryBase(i)->cross_product_transforms_as(q1+Offset, q2+Offset, qOut+Offset);
   }
}

std::complex<double>
SymmetryList::cross_product_factor(sc_iter q1, sc_iter q2) const
{
   DEBUG_CHECK(pImpl);
   // the result is the product of the cross product factor of each symmetry
   std::complex<double> d = 1.0;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      d *= pImpl->GetSymmetryBase(i)->cross_product_factor(q1+pImpl->QuantumNumberOffset(i),
                                                           q2+pImpl->QuantumNumberOffset(i));
   }
   return d;
}

double SymmetryList::recoupling(sc_iter q1, sc_iter q3, sc_iter q13,
                                sc_iter q2, sc_iter q,  sc_iter q23) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->recoupling(q1+Offset, q3+Offset, q13+Offset,
                                                     q2+Offset, q+Offset,  q23+Offset);
   }
   return Value;
}


double SymmetryList::recoupling_12_3__13_2(sc_iter q1, sc_iter q3, sc_iter q13,
                                           sc_iter q2, sc_iter q,  sc_iter q23) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->recoupling_12_3__13_2(q1+Offset, q3+Offset, q13+Offset,
                                                                q2+Offset, q+Offset,  q23+Offset);
   }
   return Value;
}

double SymmetryList::product_coefficient(sc_iter k1, sc_iter k2, sc_iter k,
                                        sc_iter qp, sc_iter q, sc_iter qpp) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->product_coefficient(k1+Offset, k2+Offset, k+Offset,
                                                             qp+Offset, q+Offset,  qpp+Offset);
   }
   return Value;
}

double SymmetryList::inverse_product_coefficient(sc_iter k1, sc_iter k2, sc_iter k,
                                                 sc_iter qp, sc_iter q, sc_iter qpp) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->inverse_product_coefficient(k1+Offset, k2+Offset, k+Offset,
                                                                      qp+Offset, q+Offset,  qpp+Offset);
   }
   return Value;
}

double SymmetryList::tensor_coefficient(sc_iter k1,  sc_iter k2,  sc_iter k,
                                       sc_iter q1p, sc_iter q2p, sc_iter qp,
                                       sc_iter q1,  sc_iter q2,  sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->tensor_coefficient(k1+Offset,  k2+Offset,  k+Offset,
                                                            q1p+Offset, q2p+Offset, qp+Offset,
                                                            q1+Offset,  q2+Offset,  q+Offset);
   }
   return Value;
}

double SymmetryList::inverse_tensor_coefficient(sc_iter k1,  sc_iter k2,  sc_iter k,
                                                sc_iter q1p, sc_iter q2p, sc_iter qp,
                                                sc_iter q1,  sc_iter q2,  sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->inverse_tensor_coefficient(k1+Offset,  k2+Offset,  k+Offset,
                                                                     q1p+Offset, q2p+Offset, qp+Offset,
                                                                     q1+Offset,  q2+Offset,  q+Offset);
   }
   return Value;
}

int
SymmetryList::num_transform_targets(sc_iter q1, sc_iter q2) const
{
   DEBUG_CHECK(pImpl);
   int n = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);
      n *= pImpl->GetSymmetryBase(i)->num_transform_targets(q1+Offset, q2+Offset);
   }
   return n;
}

bool
SymmetryList::is_transform_target(sc_iter q1, sc_iter q2, sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      if (!(pImpl->GetSymmetryBase(i)->
            is_transform_target(q1+Offset, q2+Offset, q+Offset))) return false;
   }
   return true;
}

int
SymmetryList::num_inverse_transform_targets(sc_iter q1, sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   int n = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);
      n *= pImpl->GetSymmetryBase(i)->num_inverse_transform_targets(q1+Offset, q+Offset);
   }
   return n;
}

void
SymmetryList::adjoint(sc_iter q, s_iter qOut) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      pImpl->GetSymmetryBase(i)->adjoint(q+Offset, qOut+Offset);
   }
}

double
SymmetryList::adjoint_coefficient(sc_iter qp, sc_iter k, sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->adjoint_coefficient(qp+Offset, k+Offset, q+Offset);
   }
   return Value;
}

double
SymmetryList::conj_phase(sc_iter qp, sc_iter k, sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->conj_phase(qp+Offset, k+Offset, q+Offset);
   }
   return Value;
}

double
SymmetryList::clebsch_gordan(sc_iter qp,  sc_iter k,  sc_iter q,
                             sc_iter qpm, sc_iter km, sc_iter qm) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);
      int NthProjectionOffset = pImpl->ProjectionOffset(i);

      Value *= pImpl->GetSymmetryBase(i)->
         clebsch_gordan(qp+Offset, k+Offset, q+Offset,
                        qpm+NthProjectionOffset, km+NthProjectionOffset, qm+NthProjectionOffset);
   }
   return Value;
}

std::string
SymmetryList::QuantumNumberToString(sc_iter q) const
{
   DEBUG_CHECK(pImpl);
   std::string Result;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      if (i > 0) Result += ',';
      Result += pImpl->GetSymmetryBase(i)->ToString(q+Offset);
   }
   return Result;
}

void
SymmetryList::StringToQuantumNumber(std::string const& s, s_iter q) const
{
   DEBUG_CHECK(pImpl);
   std::string::const_iterator beg = s.begin(), end = s.end();

   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      std::string::const_iterator next = std::find(beg, end, ',');
      pImpl->GetSymmetryBase(i)->FromString(std::string(beg, next), q+Offset);
      beg = next;
      ++beg;
   }
}

std::string
SymmetryList::ProjectionToString(sc_iter p) const
{
   DEBUG_CHECK(pImpl);
   std::string Result;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int NthProjectionOffset = pImpl->ProjectionOffset(i);

      if (i > 0) Result += ',';
      Result += pImpl->GetSymmetryBase(i)->ToString(p+NthProjectionOffset);
   }
   return Result;
}

void
SymmetryList::StringToProjection(std::string const& s, s_iter p) const
{
   DEBUG_CHECK(pImpl);
   std::string::const_iterator beg = s.begin(), end = s.end();

   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int NthProjectionOffset = pImpl->ProjectionOffset(i);

      std::string::const_iterator next = std::find(beg, end, ',');
      pImpl->GetSymmetryBase(i)->FromString(std::string(beg, next), p+NthProjectionOffset);
      beg = next;
      ++beg;
   }
}

void
SymmetryList::scalar_transforms_as(s_iter q) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);

      pImpl->GetSymmetryBase(i)->scalar_transforms_as(q+Offset);
   }
}

// helper function to add all (Name, Type) pairs from a given symmetry list into
// a container.  If any names already exist in the map, then verify that
// the corresponding types are the same.

typedef std::map<std::string, std::string> NameTypeMapType;

void AppendSymmetryNameType(NameTypeMapType& NameTypeMap, SymmetryList const& L)
{
   for (int i = 0; i < L.NumSymmetries(); ++i)
   {
      if (NameTypeMap.count(L.SymmetryName(i)) != 0)
      {
         CHECK(NameTypeMap[L.SymmetryName(i)] == L.SymmetryType(i));
      }
      else
      {
         NameTypeMap[L.SymmetryName(i)] = L.SymmetryType(i);
      }
   }
}

SymmetryList FindSymmetrySuperset(SymmetryList const& L1,
                                  SymmetryList const& L2)
{
   // early return in the case where the symmetry lists are the same
   if (L1.GetSymmetryListImpl() == L2.GetSymmetryListImpl()) return L1;

   // assemble a map holding all (name, type) pairs
   NameTypeMapType AllSymmetries;
   AppendSymmetryNameType(AllSymmetries, L1);
   AppendSymmetryNameType(AllSymmetries, L2);

   // if either L1 or L2 is already a superset, then just return that one (otherwise we
   // might re-order the quantum numbers unnecessarily).  We can do this by simply
   // checking the size.
   if (L1.NumSymmetries() == int(AllSymmetries.size())) return L1;
   if (L2.NumSymmetries() == int(AllSymmetries.size())) return L2;

   // now assemble the full name of the combined symmetry list
   std::string FullName;
   for (NameTypeMapType::const_iterator I = AllSymmetries.begin(); I != AllSymmetries.end(); ++I)
   {
      if (I != AllSymmetries.begin()) FullName += ',';
      FullName += I->first + ':' + I->second;
   }

   return SymmetryList(FullName);
}

SymmetryList FindSymmetrySuperset(SymmetryList const& L1,
                                  SymmetryList const& L2,
                                  SymmetryList const& L3)
{
   // early return in the case where the symmetry lists are the same
   if (L1.GetSymmetryListImpl() == L2.GetSymmetryListImpl()
       && L1.GetSymmetryListImpl() == L3.GetSymmetryListImpl()) return L1;

   // assemble a map holding all (name, type) pairs
   NameTypeMapType AllSymmetries;
   AppendSymmetryNameType(AllSymmetries, L1);
   AppendSymmetryNameType(AllSymmetries, L2);
   AppendSymmetryNameType(AllSymmetries, L3);

   // if any symmetry list is already a superset, then just return that one (otherwise we
   // might re-order the quantum numbers unnecessarily).  We can do this by simply
   // checking the size.
   if (L1.NumSymmetries() == int(AllSymmetries.size())) return L1;
   if (L2.NumSymmetries() == int(AllSymmetries.size())) return L2;
   if (L3.NumSymmetries() == int(AllSymmetries.size())) return L3;

   // now assemble the full name of the combined symmetry list
   std::string FullName;
   for (NameTypeMapType::const_iterator I = AllSymmetries.begin(); I != AllSymmetries.end(); ++I)
   {
      if (I != AllSymmetries.begin()) FullName += ',';
      FullName += I->first + ':' + I->second;
   }

   return SymmetryList(FullName);
}

SymmetryList FindSymmetrySuperset(SymmetryList const& L1,
                                  SymmetryList const& L2,
                                  SymmetryList const& L3,
                                  SymmetryList const& L4)
{
   // early return in the case where the symmetry lists are the same
   if (L1.GetSymmetryListImpl() == L2.GetSymmetryListImpl()
       && L1.GetSymmetryListImpl() == L3.GetSymmetryListImpl()
       && L1.GetSymmetryListImpl() == L4.GetSymmetryListImpl()) return L1;

   // assemble a map holding all (name, type) pairs
   NameTypeMapType AllSymmetries;
   AppendSymmetryNameType(AllSymmetries, L1);
   AppendSymmetryNameType(AllSymmetries, L2);
   AppendSymmetryNameType(AllSymmetries, L3);
   AppendSymmetryNameType(AllSymmetries, L4);

   // if any symmetry list is already a superset, then just return that one (otherwise we
   // might re-order the quantum numbers unnecessarily).  We can do this by simply
   // checking the size.
   if (L1.NumSymmetries() == int(AllSymmetries.size())) return L1;
   if (L2.NumSymmetries() == int(AllSymmetries.size())) return L2;
   if (L3.NumSymmetries() == int(AllSymmetries.size())) return L3;
   if (L4.NumSymmetries() == int(AllSymmetries.size())) return L4;

   // now assemble the full name of the combined symmetry list
   std::string FullName;
   for (NameTypeMapType::const_iterator I = AllSymmetries.begin(); I != AllSymmetries.end(); ++I)
   {
      if (I != AllSymmetries.begin()) FullName += ',';
      FullName += I->first + ':' + I->second;
   }

   return SymmetryList(FullName);
}

SymmetryList ChangeSymmetryName(SymmetryList const& L,
                                std::string const& OldName,
                                std::string const& NewName)
{
   std::string NewList;
   for (int i = 0; i < L.NumSymmetries(); ++i)
   {
      std::string Name = L.SymmetryName(i);
      if (Name == OldName)
         Name = NewName;
      if (i != 0)
         NewList += ',';
      NewList += Name + ':' + L.SymmetryType(i);
   }
   return SymmetryList(NewList);
}

void Replicate(int const* Source, int SourceSize, int SourceCount, int* Out, int Stride, int RepeatCount, int ReplicaCount)
{
   int const* SourceEnd = Source + SourceSize * SourceCount;
   for (int Replica = 0; Replica < ReplicaCount; ++Replica)
   {
      for (int const* LocalPos = Source; LocalPos != SourceEnd; LocalPos += SourceSize)
      {
         for (int Repeat = 0; Repeat < RepeatCount; ++Repeat)
         {
            memcpy(Out, LocalPos, SourceSize * sizeof(int));
            Out += Stride;
         }
      }
   }
}

void
SymmetryList::transform_targets(sc_iter q1, sc_iter q2, std::vector<int>& RawList) const
{
   DEBUG_CHECK(pImpl);
   // Construct a linear buffer to accept the packed list of quantum numbers
   int NumTargets = this->num_transform_targets(q1, q2);
   int Size = this->QuantumNumberSize();

   // Assemble the quantum numbers into a raw vector
   RawList.resize(NumTargets * Size);
   if (NumTargets == 0) return;

   // make a vector of (conservatively) maximum size possible to keep the packed form for each symmetry
   std::vector<int> Targets(RawList.size());
   int* TargetsBegin = &Targets[0];

   // Now assemble the quantum numbers.
   int RepeatCount = 1;
   int ReplicaCount = NumTargets;
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ThisSize = this->pImpl->QuantumNumberSize(i);
      int ThisNum = this->pImpl->GetSymmetryBase(i)->num_transform_targets(q1+Offset, q2+Offset);
      ReplicaCount /= ThisNum;

      this->pImpl->GetSymmetryBase(i)->transform_targets(q1+Offset, q2+Offset, TargetsBegin);
      Replicate(TargetsBegin, ThisSize, ThisNum, &RawList[Offset], Size, RepeatCount, ReplicaCount);

      RepeatCount *= ThisNum;
   }
}

void
SymmetryList::inverse_transform_targets(sc_iter q1, sc_iter q, std::vector<int>& RawList) const
{
   DEBUG_CHECK(pImpl);
   // Construct a linear buffer to accept the packed list of quantum numbers
   int NumTargets = this->num_inverse_transform_targets(q1, q);
   int Size = this->QuantumNumberSize();

   RawList.resize(NumTargets * this->QuantumNumberSize());
   if (NumTargets == 0) return;

   // make a vector of maximum size possible to keep the packed form for each symmetry
   std::vector<int> Targets(RawList.size());
   int* TargetsBegin = &Targets[0];

   // Now assemble the quantum numbers.
   int RepeatCount = 1;
   int ReplicaCount = NumTargets;
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ThisSize = this->pImpl->QuantumNumberSize(i);
      int ThisNum = this->pImpl->GetSymmetryBase(i)->
         num_inverse_transform_targets(q1+Offset, q+Offset);
      ReplicaCount /= ThisNum;

      this->pImpl->GetSymmetryBase(i)->
         inverse_transform_targets(q1+Offset, q+Offset, TargetsBegin);
      Replicate(TargetsBegin, ThisSize, ThisNum, &RawList[Offset], Size, RepeatCount, ReplicaCount);

      RepeatCount *= ThisNum;
   }
}

void
SymmetryList::enumerate_projections(sc_iter q, std::vector<int>& RawList) const
{
   DEBUG_CHECK(pImpl);
   // Construct a linear buffer to accept the packed list of projections
   int degree = this->degree(q);
   int Size = this->ProjectionSize();

   RawList.resize(degree * this->ProjectionSize());

   // make a vector of maximum size possible to keep the packed form for each symmetry
   std::vector<int> Targets(RawList.size());
   int* TargetsBegin = &Targets[0];

   // Now assemble the quantum numbers.
   int RepeatCount = 1;
   int ReplicaCount = degree;
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int ThisSize = this->pImpl->ProjectionSize(i);
      if (ThisSize == 0) continue; // Thisdegree == 1 in this case

      int ProjectionOffset = this->pImpl->ProjectionOffset(i);
      int QuantumNumberOffset = this->pImpl->QuantumNumberOffset(i);
      int Thisdegree = this->pImpl->GetSymmetryBase(i)->degree(q+QuantumNumberOffset);
      ReplicaCount /= Thisdegree;

      this->pImpl->GetSymmetryBase(i)->enumerate_projections(q+QuantumNumberOffset, TargetsBegin);
      Replicate(TargetsBegin, ThisSize, Thisdegree, &RawList[ProjectionOffset], Size, RepeatCount, ReplicaCount);

      RepeatCount *= Thisdegree;
   }
}

bool
SymmetryList::is_delta(sc_iter q1, sc_iter Q, sc_iter P, sc_iter q2) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ProjOffset = this->pImpl->ProjectionOffset(i);

      if (!this->pImpl->GetSymmetryBase(i)->is_delta(q1+Offset, Q+Offset, P+ProjOffset, q2+Offset))
        return false;
   }
   return true;
}

void SymmetryList::difference(sc_iter q1, sc_iter q2, s_iter P) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      this->pImpl->GetSymmetryBase(i)->difference(q1+Offset, q2+Offset, P+ProjOffset);
   }
}

void SymmetryList::negate(s_iter p) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      this->pImpl->GetSymmetryBase(i)->negate(p+ProjOffset);
   }
}

void SymmetryList::sum(sc_iter p1, sc_iter p2, s_iter r) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      this->pImpl->GetSymmetryBase(i)->sum(p1+ProjOffset, p2+ProjOffset, r+ProjOffset);
   }
}

bool SymmetryList::is_projection(sc_iter q, sc_iter p) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      if (!this->pImpl->GetSymmetryBase(i)->is_projection(q+Offset, p+ProjOffset)) return false;
   }
   return true;
}

bool SymmetryList::is_possible(sc_iter q, sc_iter p) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      if (!this->pImpl->GetSymmetryBase(i)->is_possible(q+Offset, p+ProjOffset)) return false;
   }
   return true;
}

void SymmetryList::change(sc_iter q, sc_iter p, s_iter Q) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      this->pImpl->GetSymmetryBase(i)->change(q+Offset, p+ProjOffset, Q+Offset);
   }
}

void SymmetryList::heighest_weight(sc_iter p, s_iter q) const
{
   DEBUG_CHECK(pImpl);
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int Offset = this->pImpl->QuantumNumberOffset(i);
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      this->pImpl->GetSymmetryBase(i)->heighest_weight(p+ProjOffset, q+Offset);
   }
}

double SymmetryList::weight(sc_iter p) const
{
   DEBUG_CHECK(pImpl);
   double Ret = 0;
   for (int i = 0; i < this->pImpl->NumSymmetries(); ++i)
   {
      int ProjOffset = this->pImpl->ProjectionOffset(i);
      Ret += this->pImpl->GetSymmetryBase(i)->weight(p+ProjOffset);
   }
   return Ret;
}

double SymmetryList::delta_shift_coefficient(sc_iter qp, sc_iter k, sc_iter q,
                                             sc_iter Delta) const
{
   DEBUG_CHECK(pImpl);
   double Value = 1;
   for (int i = 0; i < pImpl->NumSymmetries(); ++i)
   {
      int Offset = pImpl->QuantumNumberOffset(i);
      Value *= pImpl->GetSymmetryBase(i)->delta_shift_coefficient(qp+Offset, k+Offset, q+Offset,
                                                                  Delta+Offset);
   }
   return Value;
}

double SymmetryList::casimir(sc_iter q, int n) const
{
   int i = 0;
   while (n >= pImpl->GetSymmetryBase(i)->num_casimir())
      n -= pImpl->GetSymmetryBase(i++)->num_casimir();

   return  pImpl->GetSymmetryBase(i)->casimir(q+pImpl->QuantumNumberOffset(i), n);
}

QuantumNumber map_projection_to_quantum(Projection const& p,  SymmetryList const& SL)
{
   CHECK_EQUAL(p.GetSymmetryList().ProjectionSize(), SL.QuantumNumberSize());
   QuantumNumber q(SL);
   std::copy(p.begin(), p.end(), q.begin());
   return q;
}

} // namespace QuantumNumbers
