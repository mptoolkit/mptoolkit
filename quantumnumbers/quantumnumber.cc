// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/quantumnumber.cc
//
// Copyright (C) 2004-2022 Ian McCulloch <ian@qusim.net>
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

#include "common/poolallocator.h"

namespace QuantumNumbers
{

//
// RepLabelBase
//

template <typename Tag>
inline
RepLabelBase<Tag>::RepLabelBase() noexcept
{
   Storage.SList = NULL;
}

template <typename Tag>
inline
RepLabelBase<Tag>::RepLabelBase(RepLabelBase const& q)
   : Storage(q.Storage)
{
}

template <typename Tag>
inline
RepLabelBase<Tag>::RepLabelBase(SymmetryListImpl const* q, int Size)
{
   DEBUG_PRECONDITION(q != NULL);
   CHECK(Size <= QUANTUM_NUMBER_FIXED_SIZE)(Size)(QUANTUM_NUMBER_FIXED_SIZE)
      ("Increase the constant QUANTUM_NUMBER_FIXED_SIZE!");
   Storage.SList = q;
}

template <typename Tag>
template <typename InputIter>
inline
RepLabelBase<Tag>::RepLabelBase(SymmetryListImpl const* q, int Size, InputIter InitIter)
{
   static_assert(std::is_same<typename std::iterator_traits<InputIter>::value_type, int>::value,
                 "quantum number iterator value must be int");
   DEBUG_PRECONDITION(q != NULL);
   CHECK(Size <= QUANTUM_NUMBER_FIXED_SIZE)(Size)
      ("Increase the constant QUANTUM_NUMBER_FIXED_SIZE!");
   Storage.SList = q;
   iterator I = begin();
   for (int i = 0; i < Size; ++i)
   {
      *I = *InitIter;
      ++I; ++InitIter;
   }
}

template <typename Tag>
inline
RepLabelBase<Tag>::~RepLabelBase() noexcept
{
}

template <typename Tag>
RepLabelBase<Tag>&
RepLabelBase<Tag>::operator=(RepLabelBase<Tag> const& q)
{
   Storage = q.Storage;
   return *this;
}

template <typename Tag>
RepLabelBase<Tag>&
RepLabelBase<Tag>::operator=(RepLabelBase<Tag>&& q) noexcept
{
   using std::swap;
   swap(Storage, q.Storage);
   return *this;
}

template <typename Tag>
inline
bool
RepLabelBase<Tag>::is_equal_to(RepLabelBase<Tag> const& Q) const
{
   DEBUG_PRECONDITION(GetSymmetryListImpl() == Q.GetSymmetryListImpl());
   return memcmp(&Storage.NumberArray[0], &Q.Storage.NumberArray[0],
                 this->size()*sizeof(int)) == 0;
}

template <typename Tag>
inline
bool
RepLabelBase<Tag>::is_less_than(RepLabelBase<Tag> const& Q) const
{
   DEBUG_PRECONDITION(GetSymmetryListImpl() == Q.GetSymmetryListImpl());
   return memcmp(&Storage.NumberArray[0], &Q.Storage.NumberArray[0],
                 this->size()*sizeof(int)) < 0;
}

template <typename Tag>
inline
void
RepLabelBase<Tag>::WriteRaw(PStream::opstream& out) const
{
   PStream::copy_n(this->begin(), this->size(), PStream::opstream_iterator<value_type>(out));
}

template <typename Tag>
inline
void
RepLabelBase<Tag>::ReadRaw(PStream::ipstream& in)
{
   PStream::copy_n(PStream::ipstream_iterator<value_type>(in), this->size(), this->begin());
}

//
// QuantumNumber
//

inline
QuantumNumberList
adjoint(QuantumNumberList ql)
{
   for (auto& q : ql)
   {
      q = adjoint(q);
   }
   return ql;
}

inline
std::set<QuantumNumber>
adjoint(std::set<QuantumNumber> const& ql)
{
   std::set<QuantumNumber> Result;
   for (const auto& q : ql)
   {
      Result.insert(adjoint(q));
   }
   return Result;
}

inline
QuantumNumber::QuantumNumber() noexcept
{
}

inline
bool
QuantumNumber::operator==(QuantumNumber const& Q) const
{
   return is_equal_to(Q);
}

inline
bool
QuantumNumber::operator!=(QuantumNumber const& Q) const
{
   return !is_equal_to(Q);
}

inline
QuantumNumber::QuantumNumber(SymmetryList const& qL)
  : RepLabelBase<QuantumNumber>(qL.GetSymmetryListImpl(), qL.QuantumNumberSize())
{
   // set the quantum number to the scalar representation
   for (int i = 0; i < qL.NumSymmetries(); ++i)
   {
      qL.GetSymmetryBase(i)->scalar_transforms_as(begin()+qL.QuantumNumberOffset(i));
   }
}

template <typename InputIter>
inline
QuantumNumber::QuantumNumber(SymmetryList const& qList, InputIter InitIter)
  : RepLabelBase<QuantumNumber>(qList.GetSymmetryListImpl(), qList.QuantumNumberSize(), InitIter)
{
   QN_TRACE("Constructing from ")(typeid(InitIter).name());
}

inline
QuantumNumber::QuantumNumber(SymmetryList const& qList, NoInitialization)
  : RepLabelBase<QuantumNumber>(qList.GetSymmetryListImpl(), qList.QuantumNumberSize())
{
}

template <typename T>
T
QuantumNumber::get(std::string Name) const
{
   int SymmetryNumber = this->GetSymmetryList().WhichSymmetry(Name);
   if (SymmetryNumber == -1)
   {
      PANIC("Unknown quantum number")(Name);
   }

   SymmetryBase const* Sb = this->GetSymmetryList().GetSymmetryBase(SymmetryNumber);
   BasicSymmetry<T> const* SbT = dynamic_cast<BasicSymmetry<T> const*>(Sb);
   if (!SbT)
   {
      PANIC("Quantum number is not convertible to the given type")
         (Name)(typeid(T).name());
   }

   return SbT->MakeQN(this->begin()
                      + this->GetSymmetryList().QuantumNumberOffset(SymmetryNumber));
}

template <typename T>
void
QuantumNumber::set(std::string Name, T const& q)
{
   int SymmetryNumber = this->GetSymmetryList().WhichSymmetry(Name);
   if (SymmetryNumber == -1)
   {
      PANIC("Unknown quantum number")(Name);
   }

   SymmetryBase const* Sb = this->GetSymmetryList().GetSymmetryBase(SymmetryNumber);
   BasicSymmetry<T> const* SbT = dynamic_cast<BasicSymmetry<T> const*>(Sb);
   if (!SbT)
   {
      PANIC("Quantum number is not convertible to the given type")
         (Name)(typeid(T).name());
   }

   q.Convert(this->begin() + this->GetSymmetryList().QuantumNumberOffset(SymmetryNumber));
}

//
// Projection
//

template <typename InputIter>
inline
Projection::Projection(SymmetryList const& qList, InputIter InitIter)
  : RepLabelBase<Projection>(qList.GetSymmetryListImpl(), qList.ProjectionSize(), InitIter)
{
}

inline
Projection::Projection(SymmetryList const& qList, NoInitialization)
  : RepLabelBase<Projection>(qList.GetSymmetryListImpl(), qList.ProjectionSize())
{
}

inline
bool
Projection::operator==(Projection const& Q) const
{
   return is_equal_to(Q);
}

inline
bool
Projection::operator!=(Projection const& Q) const
{
   return !is_equal_to(Q);
}

//
// free functions for coupling coefficients etc
//

inline
int degree(QuantumNumber const& q)
{
   return q.degree();
}

inline
double trace(QuantumNumber const& q)
{
   return q.GetSymmetryList().trace(q.begin());
}

inline
double identity(QuantumNumber const& q)
{
   return degree(q) / trace(q);
}

inline
int multiplicity(QuantumNumber const& q1, QuantumNumber const& q2, QuantumNumber const& q)
{
   return is_transform_target(q1, q2, q);
}

inline
bool cross_product_exists(QuantumNumber const& q1, QuantumNumber const& q2)
{
   DEBUG_PRECONDITION(q1.GetSymmetryList() == q2.GetSymmetryList());
   return q1.GetSymmetryList().cross_product_exists(q1.begin(), q2.begin());
}

inline
QuantumNumber cross_product_transforms_as(QuantumNumber const& q1, QuantumNumber const& q2)
{
   DEBUG_PRECONDITION(q1.GetSymmetryList() == q2.GetSymmetryList());
   QuantumNumber Ret(q1.GetSymmetryList());
   q1.GetSymmetryList().cross_product_transforms_as(q1.begin(), q2.begin(), Ret.begin());
   return Ret;
}

inline
std::complex<double> cross_product_factor(QuantumNumber const& q1, QuantumNumber const& q2)
{
   DEBUG_PRECONDITION(q1.GetSymmetryList() == q2.GetSymmetryList());
   return q1.GetSymmetryList().cross_product_factor(q1.begin(), q2.begin());
}

inline
double clebsch_gordan(QuantumNumber const& q1, QuantumNumber const& q2, QuantumNumber const& q,
                      Projection const&    m1, Projection const&    m2, Projection const&    m)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == q1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q2.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == m1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == m2.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == m.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(multiplicity(q1, q2, q), 1)(q1)(q2)(q);
   return q.GetSymmetryList().
      clebsch_gordan(q1.begin(), q2.begin(), q.begin(), m1.begin(), m2.begin(), m.begin());
}

inline
double product_coefficient(QuantumNumber const& k1, QuantumNumber const& k2,
                           QuantumNumber const& k,
                          QuantumNumber const& qp, QuantumNumber const& q,
                           QuantumNumber const& qpp)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == qp.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == qpp.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k2.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k.GetSymmetryList());
   return q.GetSymmetryList().
      product_coefficient(k1.begin(), k2.begin(), k.begin(), qp.begin(), q.begin(), qpp.begin());
}

inline
double inverse_product_coefficient(QuantumNumber const& k1, QuantumNumber const& k2,
                                   QuantumNumber const& k,
                                   QuantumNumber const& qp, QuantumNumber const& q,
                                   QuantumNumber const& qpp)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == qp.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == qpp.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k2.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k.GetSymmetryList());
   return q.GetSymmetryList().
      inverse_product_coefficient(k1.begin(), k2.begin(), k.begin(), qp.begin(), q.begin(), qpp.begin());
}

inline
double tensor_coefficient(QuantumNumber const& k1,  QuantumNumber const& k2,
                          QuantumNumber const& k,
                         QuantumNumber const& q1p, QuantumNumber const& q2p,
                          QuantumNumber const& qp,
                         QuantumNumber const& q1,  QuantumNumber const& q2,
                          QuantumNumber const& q)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == k1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k2.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q1p.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q2p.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == qp.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q2.GetSymmetryList());
   return q.GetSymmetryList().tensor_coefficient(k1.begin()  , k2.begin() , k.begin(),
                                                q1p.begin(), q2p.begin(), qp.begin(),
                                                q1.begin() , q2.begin() , q.begin());
}

inline
double inverse_tensor_coefficient(QuantumNumber const& k1,  QuantumNumber const& k2,
                          QuantumNumber const& k,
                         QuantumNumber const& q1p, QuantumNumber const& q2p,
                          QuantumNumber const& qp,
                         QuantumNumber const& q1,  QuantumNumber const& q2,
                          QuantumNumber const& q)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == k1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k2.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == k.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q1p.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q2p.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == qp.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q2.GetSymmetryList());
   return q.GetSymmetryList().inverse_tensor_coefficient(k1.begin()  , k2.begin() , k.begin(),
                                                         q1p.begin(), q2p.begin(), qp.begin(),
                                                         q1.begin() , q2.begin() , q.begin());
}

inline
double recoupling(QuantumNumber const& q1, QuantumNumber const& q3, QuantumNumber const& q13,
                  QuantumNumber const& q2, QuantumNumber const& q,  QuantumNumber const& q23)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == q1.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q2.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q3.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q13.GetSymmetryList());
   DEBUG_PRECONDITION(q.GetSymmetryList() == q23.GetSymmetryList());
   return q.GetSymmetryList().
      recoupling(q1.begin(), q3.begin(), q13.begin(), q2.begin(), q.begin(), q23.begin());
}

inline
double recoupling_12_3__13_2(QuantumNumber const& q1, QuantumNumber const& q2,
                             QuantumNumber const& q12,
                             QuantumNumber const& q3, QuantumNumber const& q,
                             QuantumNumber const& q13)
{
   return q.GetSymmetryList().recoupling_12_3__13_2(q1.begin(), q2.begin(), q12.begin(),
                                                    q3.begin(), q.begin(), q13.begin());
}

inline
QuantumNumber adjoint(QuantumNumber const& q)
{
   QuantumNumber Ret(q.GetSymmetryList());
   q.GetSymmetryList().adjoint(q.begin(), Ret.begin());
   return Ret;
}

inline
double adjoint_coefficient(QuantumNumber const& qp, QuantumNumber const& k, QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), qp.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), k.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(multiplicity(qp,k,q), 1)(qp)(k)(q);
   return q.GetSymmetryList().adjoint_coefficient(qp.begin(), k.begin(), q.begin());
}

inline
double conj_phase(QuantumNumber const& qp, QuantumNumber const& k, QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), qp.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), k.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(multiplicity(qp,k,q), 1)(qp)(k)(q);
   return q.GetSymmetryList().conj_phase(qp.begin(), k.begin(), q.begin());
}

inline
bool
is_transform_target(QuantumNumber const& q1, QuantumNumber const& q2, QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), q1.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), q2.GetSymmetryList());
   return q.GetSymmetryList().is_transform_target(q1.begin(), q2.begin(), q.begin());
}

inline
int
num_transform_targets(QuantumNumber const& q1, QuantumNumber const& q2)
{
   DEBUG_PRECONDITION_EQUAL(q1.GetSymmetryList(), q2.GetSymmetryList());
   return q1.GetSymmetryList().num_transform_targets(q1.begin(), q2.begin());
}

template <typename OutIter>
inline
void
transform_targets(QuantumNumber const& q1, QuantumNumber const& q2, OutIter Out)
{
   DEBUG_PRECONDITION_EQUAL(q1.GetSymmetryList(), q2.GetSymmetryList());
   std::vector<int> RetList;
   q1.GetSymmetryList().transform_targets(q1.begin(), q2.begin(), RetList);
   int Size = q1.GetSymmetryList().QuantumNumberSize();
   for (std::vector<int>::const_iterator I = RetList.begin(); I != RetList.end(); I += Size)
   {
      *Out++ = QuantumNumber(q1.GetSymmetryList(), I);
   }
}

inline
int
num_inverse_transform_targets(QuantumNumber const& q1, QuantumNumber const& q2)
{
   DEBUG_PRECONDITION(q1.GetSymmetryList() == q2.GetSymmetryList());
   return q1.GetSymmetryList().num_inverse_transform_targets(q1.begin(), q2.begin());
}

template <typename OutIter>
inline
void
inverse_transform_targets(QuantumNumber const& q1, QuantumNumber const& q, OutIter Out)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == q1.GetSymmetryList());
   std::vector<int> RetList;
   q.GetSymmetryList().inverse_transform_targets(q1.begin(), q.begin(), RetList);
   int Size = q1.GetSymmetryList().QuantumNumberSize();
   for (std::vector<int>::const_iterator I = RetList.begin(); I != RetList.end(); I += Size)
   {
      *Out++ = QuantumNumber(q1.GetSymmetryList(), I);
   }
}

template <typename OutIter>
inline
void
enumerate_projections(QuantumNumber const& q, OutIter Out)
{
  //  std::cout << "quantumnumber.cc: enumerating projections of " << q << std::endl;
   std::vector<int> RetList;
   q.GetSymmetryList().enumerate_projections(q.begin(), RetList);
   //   std::cout << "RetList size is " << RetList.size() << std::endl;
   int Size = q.GetSymmetryList().ProjectionSize();
   for (std::vector<int>::const_iterator I = RetList.begin(); I != RetList.end(); I += Size)
   {
      *Out++ = Projection(q.GetSymmetryList(), I);
      //      std::cout << "Projection is " << Projection(q.GetSymmetryList(), I) << std::endl;
   }
}

inline
bool is_projection(QuantumNumber const& q, Projection const& p)
{
   DEBUG_PRECONDITION(q.GetSymmetryList() == p.GetSymmetryList());
   return q.GetSymmetryList().is_projection(q.begin(), p.begin());
}

inline
bool is_delta(QuantumNumber const& q1, QuantumNumber const& Q,
              Projection const& P, QuantumNumber const& q2)
{
   DEBUG_PRECONDITION(q1.GetSymmetryList() == Q.GetSymmetryList());
   DEBUG_PRECONDITION(q2.GetSymmetryList() == Q.GetSymmetryList());
   DEBUG_PRECONDITION(P.GetSymmetryList() == Q.GetSymmetryList());
   return Q.GetSymmetryList().is_delta(q1.begin(), Q.begin(), P.begin(), q2.begin());
}

inline
Projection difference(QuantumNumber const& q1, QuantumNumber const& q2)
{
   DEBUG_PRECONDITION_EQUAL(q1.GetSymmetryList(), q2.GetSymmetryList());
   Projection p(q1.GetSymmetryList());
   q1.GetSymmetryList().difference(q1.begin(), q2.begin(), p.begin());
   return p;
}

inline
Projection negate(Projection const& p)
{
   Projection r = p;
   r.GetSymmetryList().negate(r.begin());
   return r;
}

inline
Projection sum(Projection const& p1, Projection const& p2)
{
   Projection r(p1.GetSymmetryList());
   r.GetSymmetryList().sum(p1.begin(), p2.begin(), r.begin());
   return r;
}

inline
bool is_possible(QuantumNumber const& q, Projection const& p)
{
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), p.GetSymmetryList());
   return q.GetSymmetryList().is_possible(q.begin(), p.begin());
}

inline
QuantumNumber change(QuantumNumber const& q, Projection const& p)
{
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), p.GetSymmetryList());
   QuantumNumber Q(q.GetSymmetryList());
   q.GetSymmetryList().change(q.begin(), p.begin(), Q.begin());
   return Q;
}

inline
QuantumNumber heighest_weight(Projection const& p)
{
   QuantumNumber Ret(p.GetSymmetryList());
   p.GetSymmetryList().heighest_weight(p.begin(), Ret.begin());
   return Ret;
}

inline
double weight(Projection const& p)
{
   return p.GetSymmetryList().weight(p.begin());
}

inline
double delta_shift_coefficient(QuantumNumber const& qp, QuantumNumber const& k,
                               QuantumNumber const& q,
                               QuantumNumber const& Delta)
{
   DEBUG_PRECONDITION(k.GetSymmetryList() == q.GetSymmetryList());
   DEBUG_PRECONDITION(k.GetSymmetryList() == qp.GetSymmetryList());
   DEBUG_PRECONDITION(k.GetSymmetryList() == Delta.GetSymmetryList());
   return q.GetSymmetryList().
      delta_shift_coefficient(qp.begin(), k.begin(), q.begin(), Delta.begin());
}

inline
double casimir(QuantumNumber const& q, int n)
{
   DEBUG_PRECONDITION(n < q.GetSymmetryList().NumCasimirOperators())(n)
      (q.GetSymmetryList().NumCasimirOperators());
   return q.GetSymmetryList().casimir(q.begin(), n);
}

} // namespace QuantumNumbers
