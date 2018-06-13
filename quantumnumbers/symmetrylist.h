// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/symmetrylist.h
//
// Copyright (C) 2001-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  Defines the SymmetryList type, which uses the pImpl,
  but is not reference counted; the number of
  SymmetryListImpl instances is likely
  to remain very static, therefore the SymmetryListImpl
  maintains a static list of all instances which are deleted
  on termination of the program.

  Created 2001-11-19 Ian McCulloch
*/

#if !defined(QUANTUMNUMBERLIST_H_YR5765E78695E76UW273T67RQUFBVF5T78JFUIHKW)
#define QUANTUMNUMBERLIST_H_YR5765E78695E76UW273T67RQUFBVF5T78JFUIHKW

#include "common/trace.h"
#include "symmetrybase.h"
#include "pstream/pstream.h"
#include "common/niftycounter.h"
#include <vector>
#include <list>
#include <string>
#include <algorithm>

namespace QuantumNumbers
{

//
// SymmetryList
//
// This class maintains a list of SymmetryBase objects with the associated name.
// The SymmetryBase objects are owned by the SymmetryListImpl object
// and are destroyed by the SymmetryListImpl destructor.
// There is no attempt made to clean up the SymmetryListImpl objects
// before program termination.
//
// The name of a quantum number can by anything,
// EXCEPT it cannot include comma ',' or colon ':' characters.
//

class SymmetryListImpl
{
   public:
      // the ctor adds itself to the static list of instances
      SymmetryListImpl(std::string const& FName_);

      // appends a quantum number of the given name to the list.  N must have been allocated
      // with new, and is thereafter owned by this object.
      void Append(std::string const& Name, SymmetryBase const* N);

      // returns the size in ints required to store a quantum number
      int QuantumNumberSize() const { return count; }

      // returns the size in ints required to store a projection
      int ProjectionSize() const { return projectionCount; }

      // returns the size in ints required to store a multiplicity label.
      // This is a place-holder that always returns zero.
      int MultiplicitySize() const { return 0; }

      // returns the number of symmery groups in the list (ie, for tensor-products of symmetries)
      int NumSymmetries() const { return Data.size(); }

      // reuturns the size in ints of the i'th symmetry
      int QuantumNumberSize(int i) const
        { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Data.size())); return Data[i]->QuantumNumberSize(); }

      // returns the offset of the quantum numbers associated with the i'th symmetry
      int QuantumNumberOffset(int i) const
        { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Data.size())); return Offsets[i]; }

      // returns size in ints of the projection for the i'th symmetry
      int ProjectionSize(int i) const
        { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Data.size())); return Data[i]->ProjectionSize(); }

      // returns the offset of the projection numbers associated with the i'th symmetry
      int ProjectionOffset(int i) const
        { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Data.size())); return ProjectionOffsets[i]; }

      // returns the complete name of the quantum number list
      std::string FullName() const { return FName; }

      // returns the type of the i'th quantum number
      std::string SymmetryType(int i) const
        { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Data.size())); return Data[i]->Type(); }

      // returns the name of the i'th quantum number
      std::string SymmetryName(int i) const
        { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Names.size())); return Names[i]; }

      // returns the number of Casimir invariants of the symmetry group
      int NumCasimirOperators() const;

      // returns the name of the n'th casimir invariant.
      std::string CasimirName(int n) const;

      // returns the i'th SymmetryBase object
      SymmetryBase const* GetSymmetryBase(int i) const
        { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Data.size())); return Data[i]; }

      // returns the index of the symmetry with the given name, or -1 if there is no symmetry with that name.
      int WhichSymmetry(std::string const& Name) const
        { std::vector<std::string>::const_iterator Where = std::find(Names.begin(), Names.end(), Name);
        return Where == Names.end() ? -1 : Where-Names.begin(); }

      // deletes the owned SymmetryBase objects and detaches itself from the instance list
      ~SymmetryListImpl() noexcept;

      // searches the static instance list for an already created list of the same name.
      static SymmetryListImpl* SearchForCreated(std::string const& Name);

      // Initialize the Instances list -- called at startup via a NiftyCounter
      static void InitializeInstances();

   private:
      SymmetryListImpl(SymmetryListImpl const&); // not implemented
      SymmetryListImpl& operator=(SymmetryListImpl const&); // not implemented

      typedef std::list<SymmetryListImpl*> InstanceListType;

      std::string FName;     // full name of the symmetry list
      int count, projectionCount;  // size in ints of the memory required for a quantum number/projection
      std::vector<SymmetryBase const*> Data;
      std::vector<std::string> Names;  // the names of the symmetries
      // The Offset array is a set of offsets into the raw data representing each component of the quantum number
      std::vector<int> Offsets;
      std::vector<int> ProjectionOffsets;
      InstanceListType::iterator MyInstance;  // iterator into the list of all instances

      // maintain a static list of all instances
      static InstanceListType* Instances;
};

class SymmetryList
{
   public:
      // iterators into the quantum number/projection storage.  These MUST coincide
      // with RepLabelBase<T>::iterator and RepLabelBase<T>::const_iterator
      typedef int*       s_iter;
      typedef int const* sc_iter;

      SymmetryList();
      SymmetryList(SymmetryList const& QList) = default;
      SymmetryList(SymmetryList&& Other) noexcept;

      explicit SymmetryList(SymmetryListImpl const* Impl);

      ~SymmetryList() noexcept = default;

      // constructs a list from the given Name, which is of the form N:U(1),S:SU(2) etc.
      SymmetryList(std::string const& List);
      SymmetryList(char const* List);

      SymmetryList& operator=(SymmetryList const& QList) = default;
      SymmetryList& operator=(SymmetryList&& QList) noexcept;

      bool is_null() const { return pImpl == NULL; }

      // returns the size in ints required to store a quantum number
      int QuantumNumberSize() const
         { DEBUG_CHECK(pImpl); return pImpl->QuantumNumberSize(); }

      // returns the size in ints required to store a projection
      int ProjectionSize() const
         { DEBUG_CHECK(pImpl); return pImpl->ProjectionSize(); }

      // returns the size in ints required to store a multiplicity label.
      int MultiplicitySize() const
         { DEBUG_CHECK(pImpl); return pImpl->MultiplicitySize(); }

      // returns the number of symmery groups in the list (ie, for tensor-products of symmetries)
      int NumSymmetries() const
         { DEBUG_CHECK(pImpl); return pImpl->NumSymmetries(); }

      // returns the i'th SymmetryBase object
      SymmetryBase const* GetSymmetryBase(int i) const
         { DEBUG_CHECK(pImpl); return pImpl->GetSymmetryBase(i); }

      // returns the index of the symmetry with the given name,
      // or -1 if there is no symmetry with that name.
      int WhichSymmetry(std::string const& Name) const
         { DEBUG_CHECK(pImpl); return pImpl->WhichSymmetry(Name); }

      // reuturns the size in ints of the i'th quantum number
      int QuantumNumberSize(int i) const
         { DEBUG_CHECK(pImpl); return pImpl->QuantumNumberSize(i); }

      // returns the offset of the given quantum number
      int QuantumNumberOffset(int i) const
         { DEBUG_CHECK(pImpl); return pImpl->QuantumNumberOffset(i); }

      // reuturns the size in ints of the i'th projection
      int ProjectionSize(int i) const
         { DEBUG_CHECK(pImpl); return pImpl->ProjectionSize(i); }

      // returns the offset of the i'th projection
      int ProjectionOffset(int i) const
         { DEBUG_CHECK(pImpl); return pImpl->ProjectionOffset(i); }

      int NumCasimirOperators() const
         { DEBUG_CHECK(pImpl); return pImpl->NumCasimirOperators(); }

      std::string CasimirName(int n) const
         { DEBUG_CHECK(pImpl); return pImpl->CasimirName(n); }

      // returns the complete name of the quantum number list
      std::string FullName() const { if (!pImpl) return std::string(); return pImpl->FullName(); }

      // returns the type of the i'th quantum number
      std::string SymmetryType(int i) const { DEBUG_CHECK(pImpl); return pImpl->SymmetryType(i); }

      // returns the name of the i'th quantum number
      std::string SymmetryName(int i) const { DEBUG_CHECK(pImpl); return pImpl->SymmetryName(i); }

      SymmetryListImpl const* GetSymmetryListImpl() const { return pImpl; }

      int degree(sc_iter q) const;

      double trace(sc_iter q) const;

      bool is_transform_target(sc_iter q1, sc_iter q2, sc_iter q) const;

      bool cross_product_exists(sc_iter q1, sc_iter q2) const;

      void cross_product_transforms_as(sc_iter q1, sc_iter q2, s_iter qOut) const;

      std::complex<double> cross_product_factor(sc_iter q1, sc_iter q2) const;

      double recoupling(sc_iter q1, sc_iter q3, sc_iter q13,
                        sc_iter q2, sc_iter q, sc_iter q23) const;

      double recoupling_12_3__13_2(sc_iter q1, sc_iter q3, sc_iter q13,
                                   sc_iter q2, sc_iter q, sc_iter q23) const;

      double product_coefficient(sc_iter k1, sc_iter k2, sc_iter k,
                                sc_iter qp, sc_iter q, sc_iter qpp) const;

      double inverse_product_coefficient(sc_iter k1, sc_iter k2, sc_iter k,
                                         sc_iter qp, sc_iter q, sc_iter qpp) const;

      double tensor_coefficient(sc_iter k1,  sc_iter k2,  sc_iter k3,
                                sc_iter q1p, sc_iter q2p, sc_iter qp,
                                sc_iter q1,  sc_iter q2,  sc_iter q) const;

      double inverse_tensor_coefficient(sc_iter k1,  sc_iter k2,  sc_iter k3,
                                        sc_iter q1p, sc_iter q2p, sc_iter qp,
                                        sc_iter q1,  sc_iter q2,  sc_iter q) const;

      double reflection_coefficient(sc_iter q1, sc_iter q2, sc_iter q) const;

      int num_transform_targets(sc_iter q1, sc_iter q2) const;

      void transform_targets(sc_iter q1, sc_iter q2, std::vector<int>& RetList) const;

      int num_inverse_transform_targets(sc_iter q1, sc_iter q) const;

      void inverse_transform_targets(sc_iter q1, sc_iter q, std::vector<int>& RetList) const;

      void adjoint(sc_iter q, s_iter qOut) const;

      double adjoint_coefficient(sc_iter qp, sc_iter k, sc_iter q) const;

      double conj_phase(sc_iter qp, sc_iter k, sc_iter q) const;

      void enumerate_projections(sc_iter q, std::vector<int>& RetList) const;

      double clebsch_gordan(sc_iter qp,  sc_iter k,  sc_iter q,
                            sc_iter qpm, sc_iter km, sc_iter qm) const;

      std::string QuantumNumberToString(sc_iter q) const;

      void StringToQuantumNumber(std::string const& s, s_iter q) const;

      std::string ProjectionToString(sc_iter p) const;

      void StringToProjection(std::string const& s, s_iter p) const;

      void scalar_transforms_as(s_iter q) const;

      bool is_delta(sc_iter q1, sc_iter Q, sc_iter P, sc_iter q2) const;

      void difference(sc_iter q1, sc_iter q2, s_iter P) const;

      void negate(s_iter p) const;

      void sum(sc_iter p1, sc_iter p2, s_iter p) const;

      bool is_projection(sc_iter q, sc_iter p) const;

      bool is_possible(sc_iter q, sc_iter p) const;

      void change(sc_iter q, sc_iter p, s_iter Q) const;

      void heighest_weight(sc_iter p, s_iter q) const;

      double weight(sc_iter p) const;

      double delta_shift_coefficient(sc_iter qp, sc_iter k, sc_iter q, sc_iter Delta) const;

      double casimir(sc_iter q, int n) const;

   private:

      SymmetryListImpl const* pImpl;
};

// this is for debugging; to print a sensible name, use SymmetryList::FullName() instead
std::ostream& operator<<(std::ostream& out, SymmetryList const& s);

// Given two symmetry lists, finds a superset symmetry list that contains
// both L1 and L2.  This fails if L1 and L2 have a symmetry of the same name
// but a different type.
SymmetryList FindSymmetrySuperset(SymmetryList const& L1,
                                  SymmetryList const& L2);

SymmetryList FindSymmetrySuperset(SymmetryList const& L1,
                                  SymmetryList const& L2,
                                  SymmetryList const& L3);

SymmetryList FindSymmetrySuperset(SymmetryList const& L1,
                                  SymmetryList const& L2,
                                  SymmetryList const& L3,
                                  SymmetryList const& L4);

// return a SymmetryList that is the same as the input list, but with
// the name of a symmetry changed.
SymmetryList ChangeSymmetryName(SymmetryList const& L,
                                std::string const& OldName,
                                std::string const& NewName);

inline
bool operator==(SymmetryList const& L1, SymmetryList const& L2)
{
   return L1.GetSymmetryListImpl() == L2.GetSymmetryListImpl();
}

inline
bool operator!=(SymmetryList const& L1, SymmetryList const& L2)
{
   return L1.GetSymmetryListImpl() != L2.GetSymmetryListImpl();
}

PStream::opstream& operator<<(PStream::opstream& out, SymmetryList const& L);
PStream::ipstream& operator>>(PStream::ipstream& in, SymmetryList& L);

//
// inlines
//

inline
SymmetryList::SymmetryList(SymmetryList&& QList) noexcept
  : pImpl(QList.pImpl)
{
   QList.pImpl = nullptr;
}

inline
SymmetryList::SymmetryList(SymmetryListImpl const* Impl)
  : pImpl(Impl)
{
}

inline
SymmetryList&
SymmetryList::operator=(SymmetryList&& QList) noexcept
{
   std::swap(pImpl, QList.pImpl);
}

namespace
{
   NiftyCounter::nifty_counter<SymmetryListImpl::InitializeInstances> SymmetryListInitCounter;
} // namespace

} // namespace QuantumNumbers

#include "symmetrylist.cc"

#endif
