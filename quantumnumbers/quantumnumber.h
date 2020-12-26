// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/quantumnumber.h
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
  an optimized quantum number type.

  Created 2001-05-17 Ian McCulloch
*/

#if !defined(MPTOOLKIT_QUANTUMNUMBERS_QUANTUMNUMBER_H)
#define MPTOOLKIT_QUANTUMNUMBERS_QUANTUMNUMBER_H

/*
  rationale:
  The original method was to store quantum numbers as strings.  This very flexible,
  but unfortunately too slow, and with some disadvantages (for example, the quantum numbers
  don't know which symmetry group they are representing).
  This seeks to remedy this by storing quantum numbers directly without need for
  string conversions.

  The proper OO way of doing this would be to store each quantum number as a vector
  of pointers to an abstract base class,
  with virtual functions for obtaining the coupling coefficients etc.  But copying
  these objects would still be pretty slow.

  The compromise position is to store quantum numbers as a raw block of memory, such
  that the copy constructor can do a bit-for-bit copy of the memory.  This places
  some restrictions on how quantum numbers are stored, but in practice quantum numbers
  are pretty ordinary objects anyway.  In principle, they can be labelled by a single
  integer, although it is more convenient in general to use multiple integers.
  A multiple-integer representation is the intent.

  The first sizeof(SymmetryListImpl*) bytes of the block of memory for the
  quantum number is a pointer to the associated SymmetryListImpl object.

  A reasonable alternative would be to have the quantum number type as some kind of template
  type, with a polymorphic wrapper.  (Purely as a template type would be bad, it would mean
  that virtually every class in the program would become a template parameterized by the
  quantum number type.)  This has the advantage that constructing the quantum number
  (say, in GetTargetState() function) could be done directly, without converting it to
  a string and other such hackery.  In any event, this need not change the existing
  interface, so could be done later.

  Persistent streaming of quantum number types is messy.

  Storage allocation of the QuantumNumberType:

             +-----------------------+
  Storage -> | pointer to            |
             | SymmetryListImpl      |
             +-----------------------+
             |      Data             |
             |  (variable size)      |
             .                       .
             .                       .
             +-----------------------+

  2003-01-04: changed the implementation slightly to be a list of int's rather than a list of char's.

  FIXME: the storage has changed completely, the above docs are out of date.

*/

#include "symmetrylist.h"
#include "pstream/pstream.h"
#include <string>
#include <algorithm>

#define QUANTUM_NUMBER_FIXED_SIZE 4

#if defined(QN_DEBUG)
#define QN_TRACE(Msg) TRACE(Msg)
#else
#define QN_TRACE(Msg) DUMMY_TRACE(Msg)
#endif

namespace QuantumNumbers
{

//
// QuantumNumber is a generic class for labelling a representation
// associated with a particular symmetry group.
// This is enabled by a theorem of Peter and Weyl (1927):
// For a compact Lie group, the number of inequivalent irreducible
// representations is infinite but countable.
// Thus, in principle, a single integer would be sufficient.
// For convenience it is better to have an array of integers.
// The array size is fixed by the symmetry group.
// For finite groups, it is obvious that the set of all representations
// can be labelled by a finite set of integers.
//
// RepLabelBase abstracts out common components of QuantumNumber and Projection
//

template <typename Tag>
class RepLabelBase
{
   public:
      typedef int value_type;
      typedef int* iterator;
      typedef int const* const_iterator;

      // returns true if this object doesnt have an associated RepLabelBaseList
      bool is_null() const { return this->GetSymmetryListImpl() == NULL; }

      // returns the number of integers used to represent this quantum number
      size_t size() const { return static_cast<Tag const&>(*this).size(); }

      iterator begin() { return &Storage.NumberArray[0]; }
      const_iterator begin() const { return &Storage.NumberArray[0]; }

      iterator end() { return &Storage.NumberArray[0] + size(); }
      const_iterator end() const { return &Storage.NumberArray[0] + size(); }

      // returns the quantum number list object associated with this.
      // Precondition: !IsNull()
      SymmetryList GetSymmetryList() const { return SymmetryList(Storage.SList); }

      void WriteRaw(PStream::opstream& out) const;
      void ReadRaw(PStream::ipstream& in);

   protected:
      RepLabelBase() noexcept;
      RepLabelBase(RepLabelBase const& q);
      RepLabelBase(RepLabelBase&& Other) noexcept = default;
      ~RepLabelBase() noexcept;
      RepLabelBase& operator=(RepLabelBase const& q);
      RepLabelBase& operator=(RepLabelBase&& q) noexcept;

      // assigns the quantum number list, and allocates the storage array, does not initialize it
      RepLabelBase(SymmetryListImpl const* q, int Size);

      // assigns the quantum number list, and allocates the storage array and
      // initializes it from the provided input iterator
      template <typename InputIter>
      RepLabelBase(SymmetryListImpl const* q, int Size, InputIter InitIter);

      bool is_equal_to(RepLabelBase const& Q) const;

      bool is_less_than(RepLabelBase const& Q) const;

      size_t StorageSize() const
         { return this->CalculatePrivateSize(this->size()); }

      SymmetryListImpl const* GetSymmetryListImpl() const { return Storage.SList; }

      // implements swap(*this, Other)
      void DoSwap(RepLabelBase& Other)
         { std::swap(Storage, Other.Storage); }

   private:
      struct RepLabelBaseStorage
      {
         SymmetryListImpl const* SList;
         int NumberArray[QUANTUM_NUMBER_FIXED_SIZE];
      };

      static size_t CalculatePrivateSize(int Size)
         { return sizeof(RepLabelBaseStorage) + (Size - QUANTUM_NUMBER_FIXED_SIZE) * sizeof(int); }

      RepLabelBaseStorage Storage;

   friend class SymmetryListImpl;
   friend class SymmetryList;
};

class QuantumNumber : public RepLabelBase<QuantumNumber>
{
   public:
      struct NoInitialization {};  // tag class for the constructor that
                                   // allocates space but does not initialize

      QuantumNumber() noexcept;

      // Compiler-generated copy ctor/assignment are OK

      // assigns the identity quantum number to this
      explicit QuantumNumber(SymmetryList const& q);

      // constructs a quantum number and initializes it from the provided input iterator
      template <typename InputIter>
      QuantumNumber(SymmetryList const& q, InputIter InitIter);

      // construction from string
      QuantumNumber(SymmetryList const& q, std::string const& s);
      QuantumNumber(SymmetryList const& q, char const* s);
      QuantumNumber(SymmetryList const& q, char* s);  // needed on intel 6.0 - is it a broken compiler?

      // Occasionally, we want to construct a quantum number and initialize it later.
      // This allocates memory but does not initialize it.
      QuantumNumber(SymmetryList const& q, NoInitialization);

      size_t size() const { return this->GetSymmetryListImpl()->QuantumNumberSize(); }

      bool operator==(QuantumNumber const& Q) const;
      bool operator!=(QuantumNumber const& Q) const;

      bool operator<(QuantumNumber const& Q) const { return this->is_less_than(Q); }

      // returns a string representation of the quantum number
      std::string ToString() const;

      // returns the degree of the representation
      int degree() const;

      // Expands *this to have symmetry list SList.
      // Missing quantum numbers become the scalar quantum number.
      void CoerceSymmetryList(SymmetryList const& SList);

      void swap(QuantumNumber& Other) { this->DoSwap(Other); }

      // Get the quantum number component of the given Name,
      // which is of type T.  For example, given a symmetry list "N:U(1),S:SU(2)",
      // and a quantum number q = "3,0.5",
      // q.get<U1>("N") will return a U1 object with value 3,
      // and q.get<SU2>("S") will return a SU2 object with value 0.5.
      template <typename T>
      T get(std::string Name) const;

      // counterpart to get: set the given component of the quantum number
      template <typename T>
      void set(std::string Name, T const& q);
};

QuantumNumber CoerceSymmetryList(QuantumNumber const& q, SymmetryList const& SList)
   __attribute__((warn_unused_result));

inline
void
CoerceSymmetryListInPlace(QuantumNumber& q, SymmetryList const& SList)
{
   q.CoerceSymmetryList(SList);
}

inline
void
QuantumNumber::CoerceSymmetryList(SymmetryList const& SList)
{
   (*this) = ::QuantumNumbers::CoerceSymmetryList(*this, SList);
}

// so we can call CoerceSymmetryList on nested components that may or may not
// contain quantum number components
template <typename T>
inline
T CoerceSymmetryList(T const& x, SymmetryList const&)
{
   return x;
}

template <typename T>
inline
void CoerceSymmetryListInPlace(T&, SymmetryList const&)
{
}


std::ostream& operator<<(std::ostream& out, QuantumNumber const& Q);

PStream::opstream& operator<<(PStream::opstream& out, QuantumNumber const& L);
PStream::ipstream& operator>>(PStream::ipstream& in, QuantumNumber& L);

//
// QuantumNumberList
//
// An array of quantum numbers of the same symmetry type.  This greatly accelerates
// serializing the array because the SymmetryList needs to be written only once
// per container, rather than once per quantum number.
//

class QuantumNumberList
{
   private:
      typedef std::vector<QuantumNumber> ImplType;

   public:
      typedef QuantumNumber             value_type;
      typedef ImplType::const_iterator  const_iterator;
      typedef ImplType::iterator        iterator;
      typedef ImplType::reference       reference;
      typedef ImplType::const_reference const_reference;
      typedef ImplType::pointer         pointer;
      typedef ImplType::size_type       size_type;

      QuantumNumberList() {}

      explicit QuantumNumberList(size_t Size, QuantumNumber const& q = QuantumNumber())
         : Impl(Size, q) {}

      // compiler generated copy and assignment are OK

      template <typename FwdIter>
      QuantumNumberList(FwdIter first, FwdIter last) : Impl(first, last) {}

      void push_back(QuantumNumber const& q)
      { DEBUG_CHECK(Impl.empty() || q.GetSymmetryList() == Impl.back().GetSymmetryList());
      Impl.push_back(q); }

      // adds q to the array, and also coerces the symmetry lists.  (not yet implemented)
      // void push_back_coerce(QuantumNumber const& q);

      size_type size() const { return Impl.size(); }

      const_iterator begin() const { return Impl.begin(); }
      const_iterator end() const { return Impl.end(); }

      iterator begin() { return Impl.begin(); }
      iterator end() { return Impl.end(); }

      value_type const& operator[](size_type x) const { return Impl[x]; }
      value_type& operator[](size_type x) { return Impl[x]; }

      value_type const& back() const { return Impl.back(); }
      value_type const& front() const { return Impl.front(); }

      void clear() { Impl.clear(); }

      SymmetryList GetSymmetryList() const
         { if (Impl.empty()) return SymmetryList();
           return Impl.begin()->GetSymmetryList(); }

      bool operator==(QuantumNumberList const& Other) const { return Impl == Other.Impl; }
      bool operator!=(QuantumNumberList const& Other) const { return Impl != Other.Impl; }

      void delta_shift(QuantumNumber const& q);

      void CoerceSymmetryList(SymmetryList const& sl);

   private:
      ImplType Impl;

   friend PStream::opstream& operator<<(PStream::opstream& out, QuantumNumberList const& q);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, QuantumNumberList& q);
};

inline
std::ostream& operator<<(std::ostream& out, QuantumNumberList const& QL)
{
   std::copy(QL.begin(), QL.end(), std::ostream_iterator<QuantumNumber>(out, " "));
   return out;
}

inline
void
QuantumNumberList::CoerceSymmetryList(SymmetryList const& sl)
{
   for (QuantumNumberList::iterator I = this->begin(); I != this->end(); ++I)
   {
      I->CoerceSymmetryList(sl);
   }
}

inline
QuantumNumberList
CoerceSymmetryList(QuantumNumberList const& b, SymmetryList const& sl)
{
   QuantumNumberList Result;
   for (QuantumNumberList::const_iterator I = b.begin(); I != b.end(); ++I)
   {
      Result.push_back(CoerceSymmetryList(*I, sl));
   }
   return Result;
}

QuantumNumberList adjoint(QuantumNumberList ql);

std::set<QuantumNumber> adjoint(std::set<QuantumNumber> const& ql);

#if 0
// returns a quantum number with the name numbers, but a different symmetry list.
// The only allowed differences in the symmetry list are different names,
// the symmetry types and order must be identical.
inline
QuantumNumber RenameSymmetry(QuantumNumber const& q, SymmetryList const& NewSL)
{
   // TODO: insert some assert here
   return QuantumNumber(NewSL, q.begin());
}
#endif

//
// MakeQN
//
// Helper function to convert from concrete symmetries to QuantumNumber.
//

template <typename T1>
inline
QuantumNumber MakeQN(SymmetryList SList, T1 n1)
{
   PRECONDITION_EQUAL(SList.NumSymmetries(), 1);
   PRECONDITION_EQUAL(SList.SymmetryType(0), n1.Type());
   QuantumNumber q(SList, QuantumNumber::NoInitialization());
   n1.Convert(q.begin());
   return q;
}

template <typename T1, typename T2>
inline
QuantumNumber MakeQN(SymmetryList SList, T1 n1, T2 n2)
{
   PRECONDITION_EQUAL(SList.NumSymmetries(), 2);
   PRECONDITION_EQUAL(SList.SymmetryType(0), n1.Type());
   PRECONDITION_EQUAL(SList.SymmetryType(1), n2.Type());
   QuantumNumber q(SList, QuantumNumber::NoInitialization());
   n2.Convert(n1.Convert(q.begin()));
   return q;
}

template <typename T1, typename T2, typename T3>
inline
QuantumNumber MakeQN(SymmetryList SList, T1 n1, T2 n2, T3 n3)
{
   PRECONDITION(SList.NumSymmetries() == 3);
   PRECONDITION(SList.SymmetryType(0) == n1.Type());
   PRECONDITION(SList.SymmetryType(1) == n2.Type());
   PRECONDITION(SList.SymmetryType(2) == n3.Type());
   QuantumNumber q(SList, QuantumNumber::NoInitialization());
   n3.Convert(n2.Convert(n1.Convert(q.begin())));
   return q;
}

template <typename T1, typename T2, typename T3, typename T4>
inline
QuantumNumber MakeQN(SymmetryList SList, T1 n1, T2 n2, T3 n3, T4 n4)
{
   PRECONDITION(SList.NumSymmetries() == 4);
   PRECONDITION(SList.SymmetryType(0) == n1.Type());
   PRECONDITION(SList.SymmetryType(1) == n2.Type());
   PRECONDITION(SList.SymmetryType(2) == n3.Type());
   PRECONDITION(SList.SymmetryType(3) == n4.Type());
   QuantumNumber q(SList, QuantumNumber::NoInitialization());
   n4.Convert(n3.Convert(n2.Convert(n1.Convert(q.begin()))));
   return q;
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
inline
QuantumNumber MakeQN(SymmetryList SList, T1 n1, T2 n2, T3 n3, T4 n4, T5 n5)
{
   PRECONDITION(SList.NumSymmetries() == 5);
   PRECONDITION(SList.SymmetryType(0) == n1.Type());
   PRECONDITION(SList.SymmetryType(1) == n2.Type());
   PRECONDITION(SList.SymmetryType(2) == n3.Type());
   PRECONDITION(SList.SymmetryType(3) == n4.Type());
   PRECONDITION(SList.SymmetryType(4) == n5.Type());
   QuantumNumber q(SList, QuantumNumber::NoInitialization());
   n5.Convert(n4.Convert(n3.Convert(n2.Convert(n1.Convert(q.begin())))));
   return q;
}

//
// QNConstructor
//
// A functor that can be used to quickly convert from concrete quantum number types
// to QuantumNumber.
//
// Typical usage (basis states for the Hubard model with SO(4) symmetry):
//    QNConstructor<SU2, SU2> MakeSO4("S:SU(2),Q:SU(2)");
//    SiteBasisType Basis(MakeSO4.GetSymmetryList());
//    Basis.Append("Spinon", MakeSO4(0.5, 0));
//    Basis.Append("Holon", MakeSO4(0, 0.5));
//
// Or, we could do this the other way around:
//    SiteBasisType Basis(SymmetryList("S:SU(2),Q:SU(2)"));
//    QNConstructor<SU2, SU2> MakeSO4(Basis.GetSymmetryList());
//    Basis.Append("Spinon", MakeSO4(0.5, 0));
//    Basis.Append("Holon", MakeSO4(0, 0.5));

struct Dummy { typedef Dummy ProjectionType; };

template <class T1, class T2 = Dummy, class T3 = Dummy, class T4 = Dummy, class T5 = Dummy>
class QNConstructor
{
   public:
      typedef typename T1::ProjectionType P1;
      typedef typename T2::ProjectionType P2;
      typedef typename T3::ProjectionType P3;
      typedef typename T4::ProjectionType P4;
      typedef typename T5::ProjectionType P5;

      QNConstructor(std::string const& Name) : SList(Name) {}

      QNConstructor(char const* c) : SList(c) {}

      QNConstructor(SymmetryList SList_) : SList(SList_) {}

      QuantumNumber operator()(T1 n1) const { return MakeQN(SList, n1); }
      QuantumNumber operator()(T1 n1, T2 n2) const { return MakeQN(SList, n1, n2); }
      QuantumNumber operator()(T1 n1, T2 n2, T3 n3) const { return MakeQN(SList, n1, n2, n3); }
      QuantumNumber operator()(T1 n1, T2 n2, T3 n3, T4 n4) const { return MakeQN(SList, n1, n2, n3, n4); }
      QuantumNumber operator()(T1 n1, T2 n2, T3 n3, T4 n4, T5 n5) const { return MakeQN(SList, n1, n2, n3, n4, n5); }

      SymmetryList GetSymmetryList() const { return SList; }

   private:
      SymmetryList SList;
};

// returns true if q is the identity
inline
bool is_scalar(QuantumNumber const& q)
{
   return q == QuantumNumber(q.GetSymmetryList());
}

// returns the degee (dimension) of the representation q
// this can be fractional, in the case of anyonic quantum numbers
int degree(QuantumNumber const& q);

inline
bool is_pointed(QuantumNumber const& q)
{
   // TODO: optimize this implementation
   return degree(q) == 1;
}

// The multiplicity of the representation.  Currently, we do not handle
// non-multiplicity-free algebras at all, so this must always be 0 or 1.  But
// one day we will maybe do SU(3) ?  And this will be useful ;)
int multiplicity(QuantumNumber const& q1, QuantumNumber const& q2, QuantumNumber const& q);

// true if the cross product between operators transforming as these quantum numbers exists.
bool cross_product_exists(QuantumNumber const& q1, QuantumNumber const& q2);

// returns the quantum number of the cross product of q1 x q2
QuantumNumber cross_product_transforms_as(QuantumNumber const& q1, QuantumNumber const& q2);

// The cross product of operators transforming as q1 x a2 is defined to be:
// cross_product_factor(q1,q2) * prod(q1,q2,cross_product_transforms_as(q1,q2))
std::complex<double> cross_product_factor(QuantumNumber const& q1, QuantumNumber const& q2);

// coupling coefficent c such that
// <q' | AB(k) | q > = sum_{q''} c * < q' | A(k1) | q'' > < q'' | B(k2) | q >
double product_coefficient(QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
                           QuantumNumber const& qp, QuantumNumber const& q, QuantumNumber const& qpp);
//{
//   return coupling_3j_phase(q, qp, k) * sqrt(degree(qpp) / degree(k))
//      * coupling_6j(qp, k1, qpp, k2, q, k);
//}

// coupling coefficent c such that a product can be decomposed as
// < q' | A(k1) | q'' > < q'' | B(k2) | q > = sum_k c * <q' | AB(k) | q >
double
inverse_product_coefficient(QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
                            QuantumNumber const& qp, QuantumNumber const& q, QuantumNumber const& qpp);

// coupling coefficient <q' | (A \otimes B)(k) | q> = c * <q1' | A(k1) | q1> <q2' | B(k2) | q2>
// precondition: (k1,k2,k), (q1p, q2p, qp), (q1, q2, q) are all valid transform triplets.
double tensor_coefficient(QuantumNumber const& q1, QuantumNumber const& q2, QuantumNumber const& q,
                          QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
                          QuantumNumber const& q1p, QuantumNumber const& q2p, QuantumNumber const& qp);

// coupling coefficient <q1' | A(k1) | q1> <q2' | B(k2) | q2> = c * <q' | (A \otimes B)(k) | q>
// precondition: (k1,k2,k), (q1p, q2p, qp), (q1, q2, q) are all valid transform triplets.
double inverse_tensor_coefficient(QuantumNumber const& q1, QuantumNumber const& q2,
                                  QuantumNumber const& q,
                                  QuantumNumber const& k1, QuantumNumber const& k2,
                                  QuantumNumber const& k,
                                  QuantumNumber const& q1p, QuantumNumber const& q2p,
                                  QuantumNumber const& qp);

// the coupling coefficient between < q1q2(q12), q3 q | q1, q2q3(q23) q >
double recoupling(QuantumNumber const& q1, QuantumNumber const& q2, QuantumNumber const& q12,
                  QuantumNumber const& q3, QuantumNumber const& q, QuantumNumber const& q23);

// the coupling coefficient between < q1q2(q12) q3 q | q1q3(q13) q2 q >
double recoupling_12_3__13_2(QuantumNumber const& q1, QuantumNumber const& q2,
                             QuantumNumber const& q12,
                             QuantumNumber const& q3, QuantumNumber const& q,
                             QuantumNumber const& q13);

// returns the conjugate of q (from the metric interpretation,
// adjoint(q) * q contains the identity representation)
QuantumNumber adjoint(QuantumNumber const& q);

// returns the coefficient c that gives <qp || adjoint(T[k]) || q> = c conj( <q || T[k] || qp> )
// precondition: multiplicity(qp,k,q) == 1
double adjoint_coefficient(QuantumNumber const& qp,
                           QuantumNumber const& k,
                           QuantumNumber const& q);

double conj_phase(QuantumNumber const& qp,
                  QuantumNumber const& k,
                  QuantumNumber const& q);

// returns the coefficient c that gives <qp || T[k] || q> = c conj( <q || adjoint(T[k]) || qp> )
// precondition: multiplicity(q,k,qp) == 1
inline
double inverse_adjoint_coefficient(QuantumNumber const& qp,
                                   QuantumNumber const& k, QuantumNumber const& q)
{
   return 1.0 / adjoint_coefficient(q,adjoint(k),qp);
}

// returns true if q is a member of the product basis q1 * q2
bool is_transform_target(QuantumNumber const& q1, QuantumNumber const& q2, QuantumNumber const& q);

// returns the number of quantum numbers in the Clebsch-Gordan expansion of q1 * q2
int num_transform_targets(QuantumNumber const& q1, QuantumNumber const& q2);

// enumerates the Clebsch-Gordan expansion of q1 * q2
template <typename OutIter>
void transform_targets(QuantumNumber const& q1, QuantumNumber const& q2, OutIter Out);

// enumerates the Clebsch-Gordan expansion of q1 * q2
inline
QuantumNumberList transform_targets(QuantumNumber const& q1, QuantumNumber const& q2)
{
   QuantumNumberList Q;
   transform_targets(q1, q2, std::back_inserter(Q));
   return Q;
}

inline
QuantumNumberList operator*(QuantumNumber const& q1, QuantumNumber const& q2)
{
   return transform_targets(q1, q2);
}

inline
QuantumNumber operator+(QuantumNumber const& q1, QuantumNumber const& q2)
{
   PRECONDITION(is_pointed(q1));
   PRECONDITION(is_pointed(q2));
   return transform_targets(q1, q2)[0];
}

// 'addition' of quantum numbers.  This is only possible if tghe quantum numbers
// are pointed (possibly through being in the abelian sector of a non-abelian symmetry).
// This is an alias for operator+
inline
QuantumNumber delta_shift(QuantumNumber const& q1, QuantumNumber const& q2)
{
   QuantumNumberList QL = transform_targets(q1, q2);
   CHECK_EQUAL(QL.size(), 1);
   return QL[0];
}

// returns the number of quantum numbers q2 such that q is in the C-G expansion of q1 * q2
int num_inverse_transform_targets(QuantumNumber const& q1, QuantumNumber const& q);

// enumerates the quantum numbers q2 such that q is in the C-G expansion of q1 * q2
template <typename OutIter>
void inverse_transform_targets(QuantumNumber const& q1, QuantumNumber const& q, OutIter Out);

inline
QuantumNumberList inverse_transform_targets(QuantumNumber const& q1, QuantumNumber const& q)
{
   QuantumNumberList Q;
   inverse_transform_targets(q1, q, std::back_inserter(Q));
   return Q;
}

double casimir(QuantumNumber const& q, int n);

} // namespace QuantumNumbers

template <typename T, typename Enable = void>
struct Adjoint
{
};

template <>
struct Adjoint<QuantumNumbers::QuantumNumber>
{
   typedef QuantumNumbers::QuantumNumber argument_type;
   typedef QuantumNumbers::QuantumNumber result_type;
   typedef QuantumNumbers::QuantumNumber value_type;
   result_type operator()(argument_type const& x) const { return QuantumNumbers::adjoint(x); }
};

#include "quantumnumber.cc"

#endif
