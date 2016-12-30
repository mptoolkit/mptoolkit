// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/symmetrybase.cc
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// $Id
#include "common/trace.h"

#include <iostream>
#include <iomanip>

namespace QuantumNumbers
{

//
// BasicSymmetry is templated on a concrete symmetry type,
// and handles unwrapping the packed integer representation
// of the quantum number to the concrete type, and delegates
// the coupling coefficent functions.
//

template <typename T>
class BasicSymmetry : public SymmetryBase
{
   public:
      typedef T                          QuantumNumberType;
      typedef typename T::ProjectionType ProjectionType;
      typedef typename T::FactoryType    FactoryType;

      explicit BasicSymmetry(FactoryType Factory_ = FactoryType());

      inline QuantumNumberType MakeQN(int const* q) const { return Factory.MakeQuantumNumber(q); }
      inline ProjectionType    MakeP(int const* p) const { return Factory.MakeProjection(p); }

      virtual std::string Type() const;

      virtual std::string ProjectionSuffix() const;

      virtual std::string ToString(int const* q) const;

      virtual void FromString(std::string const& s, int* q) const;

      virtual std::string ProjectionToString(int const* q) const;

      virtual void ProjectionFromString(std::string const& s, int* q) const;

      virtual int multiplicity(int const* q1, int const* q2, int const* q) const;

      virtual bool cross_product_exists(int const* q1, int const* q2) const;

      virtual void cross_product_transforms_as(int const* q1, int const* q2, int* q) const;

      virtual std::complex<double> cross_product_factor(int const* q1, int const* q2) const;

      virtual double recoupling(int const* q1, int const* q3, int const* q13,
                                int const* q2, int const* q,  int const* q23) const;

      virtual double recoupling_12_3__13_2(int const* q1, int const* q3, int const* q13,
                                          int const* q2, int const* q,  int const* q23) const;

      virtual double product_coefficient(int const* k1, int const* k2, int const* k,
                                        int const* qp, int const* q,  int const* qpp) const;

      virtual double inverse_product_coefficient(int const* k1, int const* k2, int const* k,
                                                 int const* qp, int const* q,  int const* qpp) const;

      virtual double tensor_coefficient(int const* k1,  int const* k2,  int const* k,
                                           int const* q1p, int const* q2p, int const* qp,
                                           int const* q1,  int const* q2,  int const* qd) const;

      virtual double inverse_tensor_coefficient(int const* k1,  int const* k2,  int const* k,
                                                int const* q1p, int const* q2p, int const* qp,
                                                int const* q1,  int const* q2,  int const* qd) const;

      virtual int num_transform_targets(int const* q1, int const* q2) const;

      virtual void transform_targets(int const* q1,
                                    int const* q2,
                                    int* Targets) const;

      virtual int num_inverse_transform_targets(int const* q1, int const* q) const;

      virtual bool is_transform_target(int const* q1, int const* q2, int const* q) const;

      virtual void inverse_transform_targets(int const* q1,
                                             int const* q,
                                             int* InverseTargets) const;

      virtual void adjoint(int const* q, int* adjointq) const;

      virtual double adjoint_coefficient(int const* qp, int const* k, int const* q) const;

      virtual double conj_phase(int const* qp, int const* k, int const* q) const;

      virtual void enumerate_projections(int const* q,
                                        int* Projections) const;

      virtual bool is_delta(int const* q1, int const* Q, int const* P, int const* q2) const;

      virtual double clebsch_gordan(int const* qp,  int const* k,  int const* q,
                                    int const* qpm, int const* km, int const* qm) const;

      virtual int degree(int const* q) const;

      virtual double trace(int const* q) const;

      virtual void scalar_transforms_as(int* q) const;

      virtual void difference(int const* q1, int const* q2, int* p) const;

      virtual void negate(int* p) const;

      virtual void sum(int const* p1, int const* p2, int* r) const;

      virtual bool is_projection(int const* q, int const* p) const;

      virtual bool is_possible(int const* q, int const* p) const;

      virtual void change(int const* q, int const* p, int* Q) const;

      virtual void heighest_weight(int const* p, int* q) const;

      virtual double weight(int const* p) const;

      virtual double delta_shift_coefficient(int const* qp, int const* k,
                                             int const* q, int const* Delta) const;

      virtual int num_casimir() const;

      virtual std::string casimir_name(std::string const& QName, int n) const;

      virtual double casimir(int const* q, int n) const;

   private:
      FactoryType Factory;
};

//
// QuantumNumberConvertIter is a proxy output iterator that converts T (where T is either a
// quantum number type or a projection type) into its packed representation as an array of int.
// The array must be preallocated to have the correct length.
//

template <typename T>
class QuantumNumberConvertIter
{
   public:
      explicit QuantumNumberConvertIter(int* p_) : p(p_) {}

      QuantumNumberConvertIter& operator++() { return *this; }
      QuantumNumberConvertIter& operator++(int) { return *this; }
      QuantumNumberConvertIter& operator*() { return *this; }
      void operator=(T const& x) { p = x.Convert(p); }

   private:
      int* p;
};

template <typename T>
BasicSymmetry<T>::BasicSymmetry(FactoryType Factory_)
  : SymmetryBase(QuantumNumberType::Size(), ProjectionType::Size()),
    Factory(Factory_)
{
}

template <typename T>
std::string
BasicSymmetry<T>::Type() const
{
   return Factory.Type();
}

template <typename T>
std::string
BasicSymmetry<T>::ProjectionSuffix() const
{
   return Factory.ProjectionSuffix();
}

template <typename T>
std::string
BasicSymmetry<T>::ToString(int const* q) const
{
   return this->MakeQN(q).ToString();
}

template <typename T>
void
BasicSymmetry<T>::FromString(std::string const& s, int* q) const
{
   QuantumNumberType qn(Factory.MakeQuantumNumber(s));
   qn.Convert(q);
}

template <typename T>
std::string
BasicSymmetry<T>::ProjectionToString(int const* p) const
{
   return this->MakeP(p).ToString();
}

template <typename T>
void
BasicSymmetry<T>::ProjectionFromString(std::string const& s, int* p) const
{
  ProjectionType P(Factory.MakeProjection(s));
  P.Convert(p);
}

//
// To survive 2-phase lookup and forward to the correct function
// at the point of instantiation of BasicSymmetry<T>, we need
// to use an *unqualified* name for the free function.
// Unfortunately, the names of the free functions are the same as the
// names of the BasicSymmetry<T> members, so we cannot simply
// forward from the member function to the free function (even though
// the signatures are completely different).  Instead, we forward
// to the ADL_xxxx functions, which in turn, forward to the free
// functions.  ...And 2-phase name binding was supposed to make life easier!
//

template <typename T>
inline
int ADL_multiplicity(T const* S, int const* q1, int const* q2, int const* q)
{
   return multiplicity(S->MakeQN(q1), S->MakeQN(q2), S->MakeQN(q));
}

template <typename T>
inline
bool
ADL_cross_product_exists(T const* S, int const* q1, int const* q2)
{
   return cross_product_exists(S->MakeQN(q1), S->MakeQN(q2));
}

template <typename T>
void
ADL_cross_product_transforms_as(T const* S, int const* q1, int const* q2, int* q)
{
   cross_product_transforms_as(S->MakeQN(q1), S->MakeQN(q2)).Convert(q);
}

template <typename T>
std::complex<double>
ADL_cross_product_factor(T const* S, int const* q1, int const* q2)
{
   return cross_product_factor(S->MakeQN(q1), S->MakeQN(q2));
}


template <typename T>
inline
double ADL_recoupling(T const* S, int const* q1, int const* q3, int const* q13,
                      int const* q2, int const* q,  int const* q23)
{
   return recoupling(S->MakeQN(q1), S->MakeQN(q3), S->MakeQN(q13),
                     S->MakeQN(q2), S->MakeQN(q),  S->MakeQN(q23));
}

template <typename T>
inline
double ADL_recoupling_12_3__13_2(T const* S, int const* q1, int const* q3, int const* q13,
                                int const* q2, int const* q,  int const* q23)
{
   return recoupling_12_3__13_2(S->MakeQN(q1), S->MakeQN(q3), S->MakeQN(q13),
                               S->MakeQN(q2), S->MakeQN(q),  S->MakeQN(q23));
}

template <typename T>
inline
double ADL_product_coefficient(T const* S, int const* k1, int const* k2, int const* k,
                              int const* qp, int const* q,  int const* qpp)
{
   return product_coefficient(S->MakeQN(k1), S->MakeQN(k2), S->MakeQN(k),
                             S->MakeQN(qp), S->MakeQN(q),  S->MakeQN(qpp));
}

template <typename T>
inline
double ADL_inverse_product_coefficient(T const* S, int const* k1, int const* k2, int const* k,
                                       int const* qp, int const* q,  int const* qpp)
{
   return inverse_product_coefficient(S->MakeQN(k1), S->MakeQN(k2), S->MakeQN(k),
                                      S->MakeQN(qp), S->MakeQN(q),  S->MakeQN(qpp));
}

template <typename T>
inline
double ADL_tensor_coefficient(T const* S, int const* k1,  int const* k2,  int const* k,
                             int const* q1p, int const* q2p, int const* qp,
                             int const* q1,  int const* q2,  int const* q)
{
   return tensor_coefficient(S->MakeQN(k1),  S->MakeQN(k2),  S->MakeQN(k),
                            S->MakeQN(q1p), S->MakeQN(q2p), S->MakeQN(qp),
                            S->MakeQN(q1),  S->MakeQN(q2),  S->MakeQN(q));
}

template <typename T>
inline
double ADL_inverse_tensor_coefficient(T const* S, int const* k1,  int const* k2,  int const* k,
                                      int const* q1p, int const* q2p, int const* qp,
                                      int const* q1,  int const* q2,  int const* q)
{
   return inverse_tensor_coefficient(S->MakeQN(k1),  S->MakeQN(k2),  S->MakeQN(k),
                                     S->MakeQN(q1p), S->MakeQN(q2p), S->MakeQN(qp),
                                     S->MakeQN(q1),  S->MakeQN(q2),  S->MakeQN(q));
}

template <typename T>
inline
int ADL_num_transform_targets(T const* S, int const* q1, int const* q2)
{
   return num_transform_targets(S->MakeQN(q1), S->MakeQN(q2));
}

template <typename T>
inline
void ADL_transform_targets(T const* S, int const* q1,
                                        int const* q2,
                                        int* Targets)
{
   transform_targets(S->MakeQN(q1), S->MakeQN(q2),
                    QuantumNumberConvertIter<typename T::QuantumNumberType>(Targets));
}

template <typename T>
inline
int ADL_num_inverse_transform_targets(T const* S, int const* q1, int const* q)
{
   return num_inverse_transform_targets(S->MakeQN(q1), S->MakeQN(q));
}

template <typename T>
inline
void ADL_inverse_transform_targets(T const* S, int const* q1,
                                 int const* q,
                                 int* InverseTargets)
{
   inverse_transform_targets(S->MakeQN(q1), S->MakeQN(q),
                             QuantumNumberConvertIter<typename T::QuantumNumberType>(InverseTargets));
}
template <typename T>
inline
bool ADL_is_transform_target(T const* S, int const* q1, int const* q2, int const* q)
{
  return is_transform_target(S->MakeQN(q1), S->MakeQN(q2), S->MakeQN(q));
}

template <typename T>
inline
void ADL_adjoint(T const* S, int const* q, int* adjointq)
{
   adjoint(S->MakeQN(q)).Convert(adjointq);
}

template <typename T>
inline
double ADL_adjoint_coefficient(T const* S, int const* qp, int const* k, int const* q)
{
   return adjoint_coefficient(S->MakeQN(qp), S->MakeQN(k), S->MakeQN(q));
}

template <typename T>
inline
double ADL_conj_phase(T const* S, int const* qp, int const* k, int const* q)
{
   return conj_phase(S->MakeQN(qp), S->MakeQN(k), S->MakeQN(q));
}

template <typename T>
inline
void ADL_enumerate_projections(T const* S, int const* q,
                              int* Projections)
{
   enumerate_projections(S->MakeQN(q), QuantumNumberConvertIter<typename T::ProjectionType>(Projections));
}

template <typename T>
inline
bool ADL_is_delta(T const* S, int const* q1, int const* Q, int const* P, int const* q2)
{
   return is_delta(S->MakeQN(q1), S->MakeQN(Q), S->MakeP(P), S->MakeQN(q2));
}

template <typename T>
inline
double ADL_clebsch_gordan(T const* S, int const* qp,  int const* k,  int const* q,
                                 int const* qpm, int const* km, int const* qm)
{
   return clebsch_gordan(S->MakeQN(qp), S->MakeQN(k), S->MakeQN(q),
                         S->MakeP(qpm), S->MakeP(km), S->MakeP(qm));
}

template <typename T>
inline
int ADL_degree(T const* S, int const* q)
{
   return degree(S->MakeQN(q));
}

template <typename T>
inline
double ADL_trace(T const* S, int const* q)
{
   return trace(S->MakeQN(q));
}

template <typename T>
inline
void ADL_difference(T const* S, int const* q1, int const* q2, int* p)
{
   difference(S->MakeQN(q1), S->MakeQN(q2)).Convert(p);
}

template <typename T>
inline
void ADL_negate(T const* S, int* p)
{
   negate(S->MakeP(p)).Convert(p);
}

template <typename T>
inline
void ADL_sum(T const* S, int const* p1, int const* p2, int* r)
{
   sum(S->MakeP(p1), S->MakeP(p2)).Convert(r);
}

template <typename T>
inline
bool ADL_is_projection(T const* S, int const* q, int const* p)
{
   return is_projection(S->MakeQN(q), S->MakeP(p));
}

template <typename T>
inline
bool ADL_is_possible(T const* S, int const* q, int const* p)
{
   return is_possible(S->MakeQN(q), S->MakeP(p));
}

template <typename T>
inline
void ADL_change(T const* S, int const* q, int const* p, int* Q)
{
   change(S->MakeQN(q), S->MakeP(p)).Convert(Q);
}

template <typename T>
inline
void ADL_heighest_weight(T const* S, int const* p, int* q)
{
   heighest_weight(S->MakeP(p)).Convert(q);
}

template <typename T>
inline
double ADL_weight(T const* S, int const* p)
{
   return weight(S->MakeP(p));
}

template <typename T>
inline
double ADL_delta_shift_coefficient(T const* S, int const* qp, int const* k, int const* q, int const* Delta)
{
   return delta_shift_coefficient(S->MakeQN(qp), S->MakeQN(k), S->MakeQN(q), S->MakeQN(Delta));
}

template <typename T>
inline
double ADL_casimir(T const* S, int const* q, int n)
{
   return casimir(S->MakeQN(q), n);
}

//
// BasicSymmetry<T>
//

template <typename T>
int BasicSymmetry<T>::multiplicity(int const* q1, int const* q2, int const* q) const
{
   return QuantumNumbers::ADL_multiplicity(this, q1, q2, q);
}


template <typename T>
bool
BasicSymmetry<T>::cross_product_exists(int const* q1, int const* q2) const
{
   return QuantumNumbers::ADL_cross_product_exists(this, q1, q2);
}

template <typename T>
void
BasicSymmetry<T>::cross_product_transforms_as(int const* q1, int const* q2, int* q) const
{
   QuantumNumbers::ADL_cross_product_transforms_as(this, q1, q2, q);
}

template <typename T>
std::complex<double>
BasicSymmetry<T>::cross_product_factor(int const* q1, int const* q2) const
{
   return QuantumNumbers::ADL_cross_product_factor(this, q1, q2);
}

template <typename T>
double BasicSymmetry<T>::recoupling(int const* q1, int const* q3, int const* q13,
                                    int const* q2, int const* q,  int const* q23) const
{
   return QuantumNumbers::ADL_recoupling(this, q1, q3, q13, q2, q,  q23);
}

template <typename T>
double BasicSymmetry<T>::recoupling_12_3__13_2(int const* q1, int const* q3, int const* q13,
                                              int const* q2, int const* q,  int const* q23) const
{
   return QuantumNumbers::ADL_recoupling_12_3__13_2(this, q1, q3, q13, q2, q,  q23);
}

template <typename T>
double BasicSymmetry<T>::product_coefficient(int const* k1, int const* k2, int const* k,
                                            int const* qp, int const* q,  int const* qpp) const
{
   return QuantumNumbers::ADL_product_coefficient(this, k1, k2, k, qp, q, qpp);
}

template <typename T>
double BasicSymmetry<T>::inverse_product_coefficient(int const* k1, int const* k2, int const* k,
                                                     int const* qp, int const* q,  int const* qpp) const
{
   return QuantumNumbers::ADL_inverse_product_coefficient(this, k1, k2, k, qp, q, qpp);
}

template <typename T>
double BasicSymmetry<T>::tensor_coefficient(int const* j1,  int const* j2,  int const* j12,
                                           int const* j3, int const* j4, int const* j34,
                                           int const* j13,  int const* j24,  int const* j) const
{
   return QuantumNumbers::ADL_tensor_coefficient(this, j1, j2, j12, j3, j4, j34, j13, j24, j);
}

template <typename T>
double BasicSymmetry<T>::inverse_tensor_coefficient(int const* j1,  int const* j2,  int const* j12,
                                                    int const* j3, int const* j4, int const* j34,
                                                    int const* j13,  int const* j24,  int const* j) const
{
   return QuantumNumbers::ADL_inverse_tensor_coefficient(this, j1, j2, j12, j3, j4, j34, j13, j24, j);
}

template <typename T>
int BasicSymmetry<T>::num_transform_targets(int const* q1, int const* q2) const
{
   return QuantumNumbers::ADL_num_transform_targets(this, q1, q2);
}

template <typename T>
void BasicSymmetry<T>::transform_targets(int const* q1,
                                        int const* q2,
                                        int* Targets) const
{
   return QuantumNumbers::ADL_transform_targets(this, q1, q2, Targets);
}

template <typename T>
int BasicSymmetry<T>::num_inverse_transform_targets(int const* q1, int const* q) const
{
   return QuantumNumbers::ADL_num_inverse_transform_targets(this, q1, q);
}

template <typename T>
void BasicSymmetry<T>::inverse_transform_targets(int const* q1,
                                               int const* q,
                                               int* InverseTargets) const
{
   return QuantumNumbers::ADL_inverse_transform_targets(this, q1, q, InverseTargets);
}
template <typename T>
bool BasicSymmetry<T>::is_transform_target(int const* q1, int const* q2, int const* q) const
{
  return QuantumNumbers::ADL_is_transform_target(this, q1, q2, q);
}

template <typename T>
void BasicSymmetry<T>::adjoint(int const* q, int* adjointq) const
{
   return QuantumNumbers::ADL_adjoint(this, q, adjointq);
}

template <typename T>
double BasicSymmetry<T>::adjoint_coefficient(int const* qp, int const* k, int const* q) const
{
   return QuantumNumbers::ADL_adjoint_coefficient(this, qp, k, q);
}

template <typename T>
double BasicSymmetry<T>::conj_phase(int const* qp, int const* k, int const* q) const
{
   return QuantumNumbers::ADL_conj_phase(this, qp, k, q);
}

template <typename T>
void BasicSymmetry<T>::enumerate_projections(int const* q,
                                            int* Projections) const
{
   QuantumNumbers::ADL_enumerate_projections(this, q, Projections);
}

template <typename T>
bool BasicSymmetry<T>::is_delta(int const* q1, int const* Q, int const* P, int const* q2) const
{
   return QuantumNumbers::ADL_is_delta(this, q1, Q, P, q2);
}

template <typename T>
double BasicSymmetry<T>::clebsch_gordan(int const* qp,  int const* k,  int const* q,
                                        int const* qpm, int const* km, int const* qm) const
{
   return QuantumNumbers::ADL_clebsch_gordan(this, qp, k, q, qpm, km, qm);
}

template <typename T>
int BasicSymmetry<T>::degree(int const* q) const
{
   return QuantumNumbers::ADL_degree(this, q);
}

template <typename T>
double BasicSymmetry<T>::trace(int const* q) const
{
   return QuantumNumbers::ADL_trace(this, q);
}


template <typename T>
void BasicSymmetry<T>::scalar_transforms_as(int* q) const
{
   Factory.MakeQuantumNumber().Convert(q);
}

template <typename T>
void BasicSymmetry<T>::difference(int const* q1, int const* q2, int* p) const
{
   return QuantumNumbers::ADL_difference(this, q1, q2, p);
}

template <typename T>
void BasicSymmetry<T>::negate(int* p) const
{
   return QuantumNumbers::ADL_negate(this, p);
}

template <typename T>
void BasicSymmetry<T>::sum(int const* p1, int const* p2, int* r) const
{
   return QuantumNumbers::ADL_sum(this, p1, p2, r);
}

template <typename T>
bool BasicSymmetry<T>::is_projection(int const* q, int const* p) const
{
   return QuantumNumbers::ADL_is_projection(this, q, p);
}

template <typename T>
bool BasicSymmetry<T>::is_possible(int const* q, int const* p) const
{
   return QuantumNumbers::ADL_is_possible(this, q, p);
}

template <typename T>
void BasicSymmetry<T>::change(int const* q, int const* p, int* Q) const
{
   return QuantumNumbers::ADL_change(this, q, p, Q);
}

template <typename T>
void BasicSymmetry<T>::heighest_weight(int const* p, int* q) const
{
   return QuantumNumbers::ADL_heighest_weight(this, p, q);
}

template <typename T>
double BasicSymmetry<T>::weight(int const* p) const
{
   return QuantumNumbers::ADL_weight(this, p);
}

template <typename T>
double BasicSymmetry<T>::delta_shift_coefficient(int const* qp, int const* k,
                                                 int const* q, int const* Delta) const
{
   return QuantumNumbers::ADL_delta_shift_coefficient(this, qp, k, q, Delta);
}

template <typename T>
int BasicSymmetry<T>::num_casimir() const
{
   return T::num_casimir();
}

template <typename T>
std::string BasicSymmetry<T>::casimir_name(std::string const& QName, int n) const
{
   return T::casimir_name(QName, n);
}

template <typename T>
double BasicSymmetry<T>::casimir(int const* q, int n) const
{
   return QuantumNumbers::ADL_casimir(this, q, n);
}

//
// Implementation of RegisterStaticSymmetry
//

template <typename T>
class StaticSymmetryFactory : public SymmetryFactory
{
   public:
      virtual SymmetryBase* AttemptCreate(std::string const& N)
      {
         return (N == T::Type()) ? new BasicSymmetry<T>() : NULL;
      }
};

template <typename T>
void RegisterStaticSymmetry()
{
   SymmetryFactory::Register(new StaticSymmetryFactory<T>());
}

} // namespace QuantumNumbers
