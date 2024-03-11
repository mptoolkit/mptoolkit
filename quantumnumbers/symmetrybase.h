// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/symmetrybase.h
//
// Copyright (C) 2001-2016 Ian McCulloch <ian@qusim.net>
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
  Defines the abstract base class SymmetryBase and the
  template <class T, class Proj> class BasicQuantumNumber
  which are used to implement the symmetry groups.

  Created 2001-11-19 Ian McCulloch, from a file reorganization
*/

#if !defined(QUANTUMNUMBERBASE_H_DHFJKH48589UFEIH54OIHFRKU43R89)
#define QUANTUMNUMBERBASE_H_DHFJKH48589UFEIH54OIHFRKU43R89

#include "pstream/pstream.h"
#include "common/stringutil.h"
#include "common/convertstring.h"
#include "common/niftycounter.h"

namespace QuantumNumbers
{

//
// SymmetryBase
//
// This is an abstract base class that defines an interface for obtaining various quantities
// from the internal representation of a quantum number.  The string representation can be
// anything, but it CANNOT include comma ',' or colon ':' characters as these are used to
// delimit multiple quantum numbers and quantum number names.  The recommended separator for
// compound quantum numbers (eg with more than one casimir operator) is semicolon ';'.
// For example, the string representation of a quantum number of type SU(3) would
// be of the form "1;2", being the D(1,2) representation of SU(3).
//

class SymmetryBase
{
   public:
      SymmetryBase(int Count_, int ProjectionCount_)
         : count(Count_), projectionCount(ProjectionCount_) {}

      // returns the number of integers required to store a quantum number label
      size_t QuantumNumberSize() const { return count; }

      // returns the number of integers required to store a projection label
      size_t ProjectionSize() const { return projectionCount; }

      // returns the type of this quantum number.
      virtual std::string Type() const = 0;

      // returns the suffix used to denote the projection quantum numbers.
      virtual std::string ProjectionSuffix() const = 0;

      // returns a string representation of the quantum number
      virtual std::string ToString(int const* q) const = 0;

      // given a string, and an output buffer, converts the string to the raw reprentation.
      virtual void FromString(std::string const& s, int* q) const = 0;

      // given a buffer of PojectionCount() integers, returns a string
      // representation of the projection
      virtual std::string ProjectionToString(int const* q) const = 0;

      // given a string, and an output buffer, converts the string to the raw reprentation.
      virtual void ProjectionFromString(std::string const& s, int* q) const = 0;

      // returns the multiplicity of the representation q in the Clebsch-Gordan series expansion
      // of q1 \otimes q2
      virtual int multiplicity(int const* q1, int const* q2, int const* q) const = 0;

      // returns true if the cross product between operators transforming as these quantum numbers exists.
      virtual bool cross_product_exists(int const* q1, int const* q2) const = 0;

      // returns the quantum number of the cross product of q1 x q2
      virtual void cross_product_transforms_as(int const* q1, int const* q2, int* q) const = 0;

      // The cross product of operators transforming as q1 x a2 is defined to be:
      // cross_product_factor(q1,q2) * prod(q1,q2,cross_product_transforms_as(q1,q2))
      virtual std::complex<double> cross_product_factor(int const* q1, int const* q2) const = 0;

      // the coupling coefficient between < q1, q2q3(q23) q | q1q3(q13) q2 q >
      virtual double recoupling(int const* q1, int const* q3, int const* q13,
                                int const* q2, int const* q,  int const* q23) const = 0;

      // the coupling coefficient between < q1, q2q3(q23) q | q1q3(q13) q2 q >
      virtual double recoupling_12_3__13_2(int const* q1, int const* q3, int const* q13,
                                           int const* q2, int const* q,  int const* q23) const = 0;

      // coupling coefficent <q' | AB(k) | q > =
      // sum_{q''} c * < q' | A(k1) | q'' > < q'' | B(k2) | q >
      virtual double product_coefficient(int const* k1, int const* k2, int const* k,
                                        int const* qp, int const* q,  int const* qpp) const = 0;

      // coupling coefficent c such that a product can be decomposed as
      // < q' | A(k1) | q'' > < q'' | B(k2) | q > = sum_k c * <q' | AB(k) | q >
      virtual double inverse_product_coefficient(int const* k1, int const* k2, int const* k,
                                                 int const* qp, int const* q,  int const* qpp) const = 0;

      // coupling coefficient <q' | (A \otimes B)(k) | q> = c * <q1' | A(k1) | q1> <q2' | B(k2) | q2>
      // precondition: (k1,k2,k), (q1p, q2p, qp), (q1, q2, q) are all valid transform triplets.
      virtual double tensor_coefficient(int const* k1,  int const* k2,  int const* k,
                                       int const* q1p, int const* q2p, int const* qp,
                                       int const* q1,  int const* q2,  int const* qd) const = 0;

      // precondition: (k1,k2,k), (q1p, q2p, qp), (q1, q2, q) are all valid transform triplets.
      virtual double inverse_tensor_coefficient(int const* k1,  int const* k2,  int const* k,
                                                int const* q1p, int const* q2p, int const* qp,
                                                int const* q1,  int const* q2,  int const* qd) const = 0;

      // returns the number of quantum numbers in the Clebsch-Gordan expansion of q1 * q2
      virtual int num_transform_targets(int const* q1, int const* q2) const = 0;

      // returns true if q is a member of the product basis q1 * q2
      virtual bool is_transform_target(int const* q1, int const* q2, int const* q) const = 0;

      // appends the quantum numbers in the product basis q1 * q2 to the vector Targets.
      virtual void transform_targets(int const* q1,
                                    int const* q2,
                                    int* Targets) const = 0;

      virtual int num_inverse_transform_targets(int const* q1, int const* q2) const = 0;

      // appends to InverseTargets the possible quantum numbers q2 such that (q1,q2) -> q
      // is a valid transform target
      virtual void inverse_transform_targets(int const* q1,
                                             int const* q,
                                             int* InverseTargets) const = 0;

      // returns the quantum number of the Hermitian adjoint.  The projections of the hermitian
      // adjoint should be the same as the projections of the original quantum number.
      virtual void adjoint(int const* q, int* adjointq) const = 0;

      // returns the coefficient c in <q' | adjoint(T(k)) | q> = c <q | T(k) | q'>
      // precondition: q', k, q are valid for the symmetry group.
      virtual double adjoint_coefficient(int const* qp, int const* k, int const* q) const = 0;

      virtual double conj_phase(int const* qp, int const* k, int const* q) const = 0;

      // enumerates all the possible projections of this quantum number and
      // appends them to Projections
      virtual void enumerate_projections(int const* q,
                                        int* Projections) const = 0;

      virtual bool is_delta(int const* q1, int const* Q, int const* P, int const* q2) const = 0;

      // gives the matrix element of the projection < q';q'_m | k;k_m | q;q_m >
      virtual double clebsch_gordan(int const* qp,  int const* k,  int const* q,
                                    int const* qpm, int const* km, int const* qm) const = 0;

      //  returns the degree of the irreducible rep
      virtual int degree(int const* q) const = 0;

      // returns the trace of the matrix element |q><q|
      virtual double trace(int const* q) const = 0;

      // returns the quantum number of a scalar operator, by inserting it into q
      virtual void scalar_transforms_as(int* q) const = 0;

      // returns the projection p such that q1 = q2+p
      virtual void difference(int const* q1, int const* q2, int* p) const = 0;

      // negate a projection
      virtual void negate(int* p) const = 0;

      virtual void sum(int const* p1, int const* p2, int* r) const = 0;

      // returns true iff q+p is a well-defined quantum number
      virtual bool is_possible(int const* q, int const* p) const = 0;

      // returns true if p is a valid projection of q
      virtual bool is_projection(int const* q, int const* p) const = 0;

      // reutrns the quantum number Q = q+p
      virtual void change(int const* q, int const* p, int* Q) const = 0;

      // returns the heightest weight state allowed by the projection p
      virtual void heighest_weight(int const* p, int* q) const = 0;

      // returns a 'weight' associated with the projection; this should be a
      // non-negative number that is zero for 'identity', and non-zero away
      // from the identity.  This is somewhat heuristic.
      virtual double weight(int const* p) const = 0;

      virtual double delta_shift_coefficient(int const* qp, int const* k, int const* q,
                                             int const* Delta) const = 0;

      virtual int num_casimir() const = 0;

      virtual std::string casimir_name(std::string const& QName, int n) const = 0;

      virtual double casimir(int const* q, int n) const = 0;

      virtual ~SymmetryBase() = 0;

      // returns an instance of SymmetryBase that will handle a quantum number of the given Name
      static SymmetryBase* Create(std::string const& Type);

      // initialize the CreatedInstances and global data, called from nifty counter
      static void InitializeInstances();

   private:
      int count, projectionCount;

      // So we don't construct an instance more than once for each symmetry type, we keep
      // a list of all objects created via Create().  Initialized via NiftyCounter
      static std::map<std::string, SymmetryBase*>* CreatedInstances;
};

template <typename T>
struct StaticQuantumNumberFactory
{
   typedef T                          QuantumNumberType;
   typedef typename T::ProjectionType ProjectionType;

   inline static QuantumNumberType MakeQuantumNumber() { return QuantumNumberType(); }
   inline static QuantumNumberType MakeQuantumNumber(int const* q) { return QuantumNumberType(q); }
   inline static QuantumNumberType MakeQuantumNumber(std::string const& s) { return QuantumNumberType(s); }

   inline static ProjectionType MakeProjection(int const* p) { return ProjectionType(p); }
   inline static ProjectionType MakeProjection(std::string const& s) { return ProjectionType(s); }

   inline static char const* Type() { return QuantumNumberType::Type(); }
   inline static char const* ProjectionSuffix() { return ProjectionType::Suffix(); }
};

//
// SymmetryFactory
//
// The system maintains a list of factory objects which are used to create
// instances of SymmetryBase*, given the name of the symmetry group.
//

class SymmetryFactory
{
   public:
      virtual ~SymmetryFactory();

      virtual SymmetryBase* AttemptCreate(std::string const& Type) = 0;

      // registers a factory.  Thereafter, F is owned by the implementation (ie. allocate it on the heap!)
      static void Register(SymmetryFactory* F);
};

// shortcut function for registering class T as a concrete quantum number.
// Useful only if the class T is
// capable of handling only one symmetry group.
template <typename T>
void RegisterStaticSymmetry();

inline
void RegisterSymmetry(SymmetryFactory* F)
{
  SymmetryFactory::Register(F);
}

namespace
{
    NiftyCounter::nifty_counter<SymmetryBase::InitializeInstances> SymmetryBaseInitCounter;
} // namespace

} // namespace QuantumNumbers

#include "symmetrybase.cc"

#endif
