// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/u1.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  Created 2000-09-07 Ian McCulloch

  The U(1) quantum number registers an object that allows the name "U(1)"
  in the SymmetryList.

  This header simply has a nifty counter that ensures that the static object is initialized.

  Modifed 2003-10-16 Ian McCulloch: It does no harm to use the half_int version always; renamed
  U(1)_half_int to U(1) and got rid of U(1)_integer.  Updated to new quantum number interface.

*/

#if !defined(U1HALFINT_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY)
#define U1HALFINT_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY

#include "common/niftycounter.h"
#include "common/halfint.h"
#include "quantumnumber.h"

namespace QuantumNumbers
{

   //class U1Projection;

class U1
{
   public:
   //      typedef U1Projection                   ProjectionType;
         typedef U1                   ProjectionType;
         typedef U1                   QuantumNumberType;
      typedef StaticQuantumNumberFactory<U1> FactoryType;

      U1() : x(0) {}
      U1(half_int x_) : x(x_) {}
      U1(int x_) : x(x_) {}
      U1(double x_) : x(x_) {}
      explicit U1(int const* InIter);
      explicit U1(std::string const& s);

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { *OutIter = x.twice(); return OutIter+1; }

      static char const* Type() { return "U(1)"; }

      static int Size() { return 1; }

      static int num_casimir() { return 1; }
      static std::string casimir_name(std::string const& QName, int) { return QName; }

      // Registation is automatic via a nifty counter
      static void Register();

      static char const* Suffix() { return ""; }

      half_int x;

      bool operator==(U1 const& y) const { return x == y.x; }
      bool operator!=(U1 const& y) const { return x != y.x; }
};

inline
std::ostream& operator<<(std::ostream& out, U1 const& s)
{
   return out << s.x;
}

typedef U1 U1Projection;

#if 0
class U1Projection
{
   public:
     typedef U1 QuantumNumberType;

      U1Projection() {}
      explicit U1Projection(int const* InIter) {};
      explicit U1Projection(std::string const& s) {};

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { return OutIter; }

      static char const* Type() { return "U(1)"; }
      static char const* Suffix() { return ""; }

      static int Size() { return 0; }
};
inline
std::ostream& operator<<(std::ostream& out, U1Projection const&)
{
   return out;
}
#endif


//
// inlines
//

inline
U1::U1(int const* InIter)
  : x(from_twice(*InIter))
{
}

inline
int degree(U1 const&)
{
   return 1;
}

inline
double trace(U1 const&)
{
   return 1;
}

inline
int multiplicity(U1 const&, U1 const&, U1 const&)
{
   return 1;
}

inline
bool cross_product_exists(U1 const&, U1 const&)
{
   return true;
}

inline
U1 cross_product_transforms_as(U1 const& q1, U1 const& q2)
{
   return U1(q1.x+q2.x);
}

inline
std::complex<double> cross_product_factor(U1 const&, U1 const&)
{
   return 1.0;
}

inline
double clebsch_gordan(U1 const& q1, U1 const& q2, U1 const& q,
                      U1Projection const& m1,  U1Projection const& m2,  U1Projection const& m)
{
   DEBUG_CHECK_EQUAL(q1,m1);
   DEBUG_CHECK_EQUAL(q2,m2);
   DEBUG_CHECK_EQUAL(q,m);
   return (q1.x + q2.x == q.x) ? 1 : 0;
}

inline
double product_coefficient(U1 const& k1, U1 const& k2, U1 const& k,
                          U1 const& qp, U1 const& q, U1 const& qpp)
{
   return ((k1.x + k2.x == k.x) && (qp.x == k1.x + qpp.x) && (qpp.x == k2.x + q.x)) ? 1 : 0;
}

inline
double inverse_product_coefficient(U1 const& k1, U1 const& k2, U1 const& k,
                          U1 const& qp, U1 const& q, U1 const& qpp)
{
   return ((k1.x + k2.x == k.x) && (qp.x == k1.x + qpp.x) && (qpp.x == k2.x + q.x)) ? 1 : 0;
}

inline
double tensor_coefficient(U1 const& q1,  U1 const& q2,  U1 const& q,
                          U1 const& k1,  U1 const& k2,  U1 const& k,
                          U1 const& q1p, U1 const& q2p, U1 const& qp)
{
   DEBUG_PRECONDITION(k1.x + k2.x == k.x)(k1.x)(k2.x)(k.x);
   DEBUG_PRECONDITION(q1p.x + q2p.x == qp.x)(q1p.x)(q2p.x)(qp.x);
   DEBUG_PRECONDITION(q1.x + q2.x == q.x)(q1.x)(q2.x)(q.x);
   DEBUG_PRECONDITION(q1p.x == k1.x + q1.x)(q1p.x)(k1.x)(q1.x);
   DEBUG_PRECONDITION(q2p.x == k2.x + q2.x)(q2p.x)(k2.x)(q2.x);
   DEBUG_PRECONDITION(qp.x == k.x + q.x)(qp.x)(k.x)(q.x);
   return 1;
}

inline
double inverse_tensor_coefficient(U1 const& q1,  U1 const& q2,  U1 const& q,
                                  U1 const& k1,  U1 const& k2,  U1 const& k,
                                  U1 const& q1p, U1 const& q2p, U1 const& qp)
{
   DEBUG_PRECONDITION(k1.x + k2.x == k.x)(k1.x)(k2.x)(k.x);
   DEBUG_PRECONDITION(q1p.x + q2p.x == qp.x)(q1p.x)(q2p.x)(qp.x);
   DEBUG_PRECONDITION(q1.x + q2.x == q.x)(q1.x)(q2.x)(q.x);
   DEBUG_PRECONDITION(q1p.x == k1.x + q1.x)(q1p.x)(k1.x)(q1.x);
   DEBUG_PRECONDITION(q2p.x == k2.x + q2.x)(q2p.x)(k2.x)(q2.x);
   DEBUG_PRECONDITION(qp.x == k.x + q.x)(qp.x)(k.x)(q.x);
   return 1;
}

inline
double recoupling(U1 const& q1, U1 const& q2, U1 const& q12,
                  U1 const& q3, U1 const& q, U1 const& q23)
{
   DEBUG_CHECK(q12.x == q1.x + q2.x)(q12.x)(q1.x)(q2.x);
   DEBUG_CHECK(q23.x == q2.x + q3.x)(q23.x)(q2.x)(q3.x)(q1.x)(q12.x)(q.x);
   DEBUG_CHECK(q.x == q1.x + q2.x + q3.x)(q.x)(q1.x)(q2.x)(q3.x);
   return 1;
}

inline
double recoupling_12_3__13_2(U1 const& q1, U1 const& q2, U1 const& q12,
                             U1 const& q3, U1 const& q, U1 const& q13)
{
   double x = (q13.x == q1.x + q3.x) && (q12.x == q1.x + q2.x) && (q.x == q1.x+q2.x+q3.x) && (q12.x + q3.x == q.x) &&
     (q13.x + q2.x == q.x) ? 1 : 0;
   DEBUG_CHECK(x==1);
   return x;
}

inline
U1 adjoint(U1 const& q)
{
  return U1(-q.x);
}

inline
double adjoint_coefficient(U1 const& qp, U1 const& k, U1 const& q)
{
   PRECONDITION(qp.x == q.x - k.x);
   return 1;
}

inline
double conj_phase(U1 const& qp, U1 const& k, U1 const& q)
{
   PRECONDITION(qp.x == q.x - k.x);
   return 1;
}

inline
bool is_transform_target(U1 const& q1, U1 const& q2, U1 const& q)
{
   return q1.x + q2.x == q.x;
}

inline
int num_transform_targets(U1 const& q1, U1 const& q2)
{
   return 1;
}

template <typename OutIter>
inline
void transform_targets(U1 const& q1, U1 const& q2, OutIter Out)
{
   *Out++ = q1.x + q2.x;
}

inline
int num_inverse_transform_targets(U1 const& q1, U1 const& q)
{
   return 1;
}

template <typename OutIter>
inline
void inverse_transform_targets(U1 const& q1, U1 const& q, OutIter Out)
{
   *Out++ = q.x - q1.x;
}

template <typename OutIter>
inline
void enumerate_projections(U1 const& q, OutIter Out)
{
   *Out++ = U1Projection(q);
}

inline
bool is_delta(U1 const& q1, U1 const& Q, U1Projection const& P, U1 const& q2)
{
   DEBUG_CHECK_EQUAL(Q.x, P.x);
   return q1.x == P.x + q2.x;
}

inline
U1Projection difference(U1 const& q1, U1 const& q2)
{
   return U1Projection(q1.x-q2.x);
}

inline
U1Projection negate(U1Projection const& p)
{
   return U1Projection(-p.x);
}

inline
U1Projection sum(U1Projection const& q1, U1Projection const& q2)
{
   return U1Projection(q1.x+q2.x);
}

inline
bool is_possible(U1 const& q, U1Projection const& p)
{
   return true;
}

inline
bool is_projection(U1 const& q, U1Projection const& p)
{
   return q.x == p.x;
}

inline
U1 change(U1 const& q, U1Projection const& p)
{
   return U1(q.x+p.x);
}

inline
U1 heighest_weight(U1Projection const& p)
{
   return U1(p.x);
}

inline
double weight(U1Projection const& p)
{
   return fabs(p.x.to_double());
}

inline
double delta_shift_coefficient(U1 const& qp, U1 const& k, U1 const& q, U1 const& Delta)
{
   PRECONDITION(qp.x == q.x + k.x)(qp.x)(q.x)(k.x);
   return 1;
}

inline
double casimir(U1 const& s, int n)
{
   CHECK_EQUAL(n, 0);
   return s.x.to_double();
}

//
// U1Projection
//

namespace
{
   NiftyCounter::nifty_counter<U1::Register> U1HalfIntCounter;
}

} // namespace QuantumNumbers

#endif
