// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/z2.h
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
// Quantum numbers for Z_2 multiplicative group

#if !defined(Z2_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY)
#define Z2_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY

#include "common/niftycounter.h"
#include "quantumnumber.h"
#include "common/halfint.h"
#include <boost/lexical_cast.hpp>

namespace QuantumNumbers
{

class Z2
{
   public:
      typedef Z2                             ProjectionType;
      typedef Z2                             QuantumNumberType;
      typedef StaticQuantumNumberFactory<Z2> FactoryType;

      Z2() : x(1) {}
      Z2(int x_) : x(x_) { DEBUG_PRECONDITION(x_ == -1 || x_ == 1)(x_); }
      explicit Z2(int const* InIter) : x(*InIter) { DEBUG_PRECONDITION(x == -1 || x == 1)(x); }
      explicit Z2( std::string const& s);

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { *OutIter = x; return OutIter+1; }

      static char const* Type() { return "Z2"; }

      static int Size() { return 1; }

      static int num_casimir() { return 1; }
      static std::string casimir_name(std::string const& QName, int)
      { return QName; }

      // Registation is automatic via a nifty counter
      static void Register();

      // suffix for treating Z2 as a projection
      static char const* Suffix() { return ""; }

      bool operator==(Z2 const& y) const { return x == y.x; }
      bool operator!=(Z2 const& y) const { return x != y.x; }

      int x;
};

inline
std::ostream& operator<<(std::ostream& out, Z2 const& s)
{
   return out << s.x;
}

inline
Z2 operator*(Z2 const& a, Z2 const& b)
{
   return Z2(a.x*b.x);
}

typedef Z2 Z2Projection;

//
// inlines
//

inline
int degree(Z2 const&)
{
   return 1;
}

inline
double trace(Z2 const&)
{
   return 1;
}

inline
int multiplicity(Z2 const& q1, Z2 const& q2, Z2 const& q)
{
   DEBUG_PRECONDITION(q1.x*q2.x == q.x);
   return 1;
}

inline
bool cross_product_exists(Z2 const&, Z2 const&)
{
   return true;
}

inline
Z2 cross_product_transforms_as(Z2 const& q1, Z2 const& q2)
{
   return Z2(q1.x*q2.x);
}

inline
std::complex<double> cross_product_factor(Z2 const&, Z2 const&)
{
   return 1.0;
}

inline
double clebsch_gordan(Z2 const& q1, Z2 const& q2, Z2 const& q,
	  Z2Projection const& m1,  Z2Projection const& m2,  Z2Projection const& m)
{
   DEBUG_CHECK_EQUAL(q1,m1);
   DEBUG_CHECK_EQUAL(q2,m2);
   DEBUG_CHECK_EQUAL(q,m);
   return (q1*q2 == q) ? 1 : 0;
}

inline
double product_coefficient(Z2 const& k1, Z2 const& k2, Z2 const& k,
			  Z2 const& qp, Z2 const& q, Z2 const& qpp)
{
   return ((k1 *k2 == k) && (qp == k1*qpp) && (qpp == k2*q)) ? 1 : 0;
}

inline
   double inverse_product_coefficient(Z2 const& k1, Z2 const& k2, Z2 const& k,
                                      Z2 const& qp, Z2 const& q, Z2 const& qpp)
{
   return ((k1*k2 == k) && (qp == k1*qpp) && (qpp == k2*q)) ? 1 : 0;
}


inline
double tensor_coefficient(Z2 const& k1,  Z2 const& k2,  Z2 const& k,
			 Z2 const& q1p, Z2 const& q2p, Z2 const& qp,
			 Z2 const& q1,  Z2 const& q2,  Z2 const& q)
{
   PRECONDITION(k1*k2 == k)(k1)(k2)(k);
   PRECONDITION(q1p*q2p == qp)(q1p)(q2p)(qp);
   PRECONDITION(q1*q2 == q)(q1)(q2)(q);
   PRECONDITION(q1p == k1*q1)(q1p)(k1)(q1);
   PRECONDITION(q2p == k2*q2)(q2p)(k2)(q2);
   PRECONDITION(qp == k*q)(qp)(k)(q);
   return 1;
}

inline
double inverse_tensor_coefficient(Z2 const& k1,  Z2 const& k2,  Z2 const& k,
                                  Z2 const& q1p, Z2 const& q2p, Z2 const& qp,
                                  Z2 const& q1,  Z2 const& q2,  Z2 const& q)
{
   PRECONDITION(k1*k2 == k)(k1)(k2)(k);
   PRECONDITION(q1p*q2p == qp)(q1p)(q2p)(qp);
   PRECONDITION(q1*q2 == q)(q1)(q2)(q);
   PRECONDITION(q1p == k1*q1)(q1p)(k1)(q1);
   PRECONDITION(q2p == k2*q2)(q2p)(k2)(q2);
   PRECONDITION(qp == k*q)(qp)(k)(q);
   return 1;
}
inline
double recoupling(Z2 const& q1, Z2 const& q3, Z2 const& q13,
		  Z2 const& q2, Z2 const& q, Z2 const& q23)
{
   DEBUG_PRECONDITION_EQUAL(q13, q1*q3);
   DEBUG_PRECONDITION_EQUAL(q23, q2*q3);
   DEBUG_PRECONDITION_EQUAL(q, q1*q23);
   return 1;
}

inline
double recoupling_12_3__13_2(Z2 const& q1, Z2 const& q2, Z2 const& q12,
                             Z2 const& q3, Z2 const& q, Z2 const& q13)
{
   DEBUG_PRECONDITION_EQUAL(q12, q1*q2);
   DEBUG_PRECONDITION_EQUAL(q13, q1*q3);
   DEBUG_PRECONDITION_EQUAL(q, q12*q3);
   return 1;
}

inline
Z2 adjoint(Z2 const& q)
{
   return q;
}

inline
double adjoint_coefficient(Z2 const& qp, Z2 const& k, Z2 const& q)
{
   PRECONDITION_EQUAL(qp, k*q);
   return 1;
}

inline
double conj_phase(Z2 const& qp, Z2 const& k, Z2 const& q)
{
   PRECONDITION_EQUAL(qp, k*q);
   return 1;
}

inline
bool is_transform_target(Z2 const& q1, Z2 const& q2, Z2 const& q)
{
   return q.x == q1.x*q2.x;
}

inline
int num_transform_targets(Z2 const& q1, Z2 const& q2)
{
   return 1;
}

template <typename OutIter>
inline
void transform_targets(Z2 const& q1, Z2 const& q2, OutIter Out)
{
   *Out++ = q1*q2;
}

inline
int num_inverse_transform_targets(Z2 const& q1, Z2 const& q)
{
   return 1;
}

template <typename OutIter>
inline
void inverse_transform_targets(Z2 const& q1, Z2 const& q, OutIter Out)
{
   *Out++ = q1*q;
}

template <typename OutIter>
inline
void enumerate_projections(Z2 const& q, OutIter Out)
{
   *Out++ = Z2Projection(q);
}

inline
bool is_delta(Z2 const& q1, Z2 const& Q, Z2Projection const& P, Z2 const& q2)
{
   DEBUG_PRECONDITION_EQUAL(Q, P);
   return q1 == P*q2;
}

inline
Z2Projection difference(Z2 const& q1, Z2 const& q2)
{
   return q1*q2;
}

inline
Z2Projection negate(Z2Projection const& p)
{
   return p;
}

inline
Z2Projection sum(Z2Projection const& q1, Z2Projection const& q2)
{
   return q1*q2;
}

inline
bool is_possible(Z2 const& q, Z2Projection const& p)
{
   return true;
}

inline
bool is_projection(Z2 const& q, Z2Projection const& p)
{
   return q == p;
}

inline
Z2 change(Z2 const& q, Z2Projection const& p)
{
   return q*p;
}

inline
Z2 heighest_weight(Z2Projection const& p)
{
   return p;
}

inline
double weight(Z2Projection const& p)
{
   return p.x == 1 ? 0 : 1;
}

inline
double delta_shift_coefficient(Z2 const& qp, Z2 const& k, Z2 const& q, Z2 const& Delta)
{
   PRECONDITION_EQUAL(qp, q*k);
   return 1;
}

inline
double casimir(Z2 const& s, int n)
{
   CHECK_EQUAL(n, 0);
   return s.x;
}

namespace
{
   NiftyCounter::nifty_counter<Z2::Register> Z2Registrator;
}

} // namespace QuantumNumbers

#endif
