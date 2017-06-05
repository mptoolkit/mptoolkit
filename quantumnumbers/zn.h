// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/zn.h
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
// Quantum numbers for Z_N multiplicative group

#if !defined(MPTOOLKIT_QUANTUMNUMBERS_ZN_H)
#define MPTOOLKIT_QUANTUMNUMBERS_ZN_H

#include "common/niftycounter.h"
#include "quantumnumber.h"
#include "common/halfint.h"
#include <boost/lexical_cast.hpp>

namespace QuantumNumbers
{

extern char const* TypeStr[];

template <int n>
class Zn
{
   public:
      typedef Zn                                ProjectionType;
      typedef Zn                                QuantumNumberType;
      typedef StaticQuantumNumberFactory<Zn<n>> FactoryType;

      Zn() : x(0) {}
      Zn(half_int x_) : x((x_+n)%n) {}
      Zn(int x_) : x((x_ + n)%n) {}
      explicit Zn(int const* InIter) : x(from_twice(*InIter)) {}
      explicit Zn( std::string const& s);

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { *OutIter = x.twice(); return OutIter+1; }

      static char const* Type() { return TypeStr[n]; }

      static int Size() { return 1; }

      static int num_casimir() { return 1; }
      static std::string casimir_name(std::string const& QName, int)
      { return QName; }

      // Registation is automatic via a nifty counter
      static void Register();

      // suffix for treating Zn as a projection
      static char const* Suffix() { return ""; }

      bool operator==(Zn const& y) const { return x == y.x; }
      bool operator!=(Zn const& y) const { return x != y.x; }

      half_int x;
};

template <int n>
inline
Zn<n>
operator+(Zn<n> x, Zn<n> y)
{
   return Zn<n>((x.x + y.x)%n);
}

template <int n>
inline
Zn<n>
operator-(Zn<n> x, Zn<n> y)
{
   return Zn<n>((x.x - y.x + n)%n);
}

template <int n>
inline
Zn<n>
operator-(Zn<n> x)
{
   return Zn<n>(n - x.x);
}

template <int n>
inline
std::ostream& operator<<(std::ostream& out, Zn<n> const& s)
{
   return out << s.x;
}

//
// inlines
//

template <int n>
inline
int degree(Zn<n> const&)
{
   return 1;
}

template <int n>
inline
double trace(Zn<n> const&)
{
   return 1;
}

template <int n>
inline
bool cross_product_exists(Zn<n> const&, Zn<n> const&)
{
   return true;
}

template <int n>
inline
Zn<n>
cross_product_transforms_as(Zn<n> const& q1, Zn<n> const& q2)
{
   return Zn<n>(q1+q2);
}

template <int n>
inline
std::complex<double>
cross_product_factor(Zn<n> const&, Zn<n> const&)
{
   return 1.0;
}

template <int n>
inline
double clebsch_gordan(Zn<n> const& q1, Zn<n> const& q2, Zn<n> const& q,
                      Zn<n> const& m1, Zn<n> const& m2, Zn<n> const& m)
{
   DEBUG_CHECK_EQUAL(q1,m1);
   DEBUG_CHECK_EQUAL(q2,m2);
   DEBUG_CHECK_EQUAL(q,m);
   return (q1+q2 == q) ? 1 : 0;
}

template <int n>
inline
double product_coefficient(Zn<n> const& k1, Zn<n> const& k2, Zn<n> const& k,
                           Zn<n> const& qp, Zn<n> const& q, Zn<n> const& qpp)
{
   return ((k1 + k2 == k) && (qp == k1 + qpp) && (qpp == k2 + q)) ? 1 : 0;
}


template <int n>
inline
double inverse_product_coefficient(Zn<n> const& k1, Zn<n> const& k2, Zn<n> const& k,
                                      Zn<n> const& qp, Zn<n> const& q, Zn<n> const& qpp)
{
   return ((k1+k2 == k) && (qp == k1+qpp) && (qpp == k2+q)) ? 1 : 0;
}

template <int n>
inline
double tensor_coefficient(Zn<n> const& k1,  Zn<n> const& k2,  Zn<n> const& k,
                         Zn<n> const& q1p, Zn<n> const& q2p, Zn<n> const& qp,
                         Zn<n> const& q1,  Zn<n> const& q2,  Zn<n> const& q)
{
   PRECONDITION(k1 + k2 == k)(k1)(k2)(k);
   PRECONDITION(q1p + q2p == qp)(q1p)(q2p)(qp);
   PRECONDITION(q1 + q2 == q)(q1)(q2)(q);
   PRECONDITION(k1 + q1p == q1)(q1p)(k1)(q1);
   PRECONDITION(k2 + q2p == q2)(q2p)(k2)(q2);
   PRECONDITION(k + qp == q)(qp)(k)(q);
   return 1;
}

template <int n>
inline
double inverse_tensor_coefficient(Zn<n> const& k1,  Zn<n> const& k2,  Zn<n> const& k,
                                  Zn<n> const& q1p, Zn<n> const& q2p, Zn<n> const& qp,
                                  Zn<n> const& q1,  Zn<n> const& q2,  Zn<n> const& q)
{
   PRECONDITION(k1 + k2 == k)(k1)(k2)(k);
   PRECONDITION(q1p + q2p == qp)(q1p)(q2p)(qp);
   PRECONDITION(q1 + q2 == q)(q1)(q2)(q);
   PRECONDITION(k1 + q1p == q1)(q1p)(k1)(q1);
   PRECONDITION(k2 + q2p == q2)(q2p)(k2)(q2);
   PRECONDITION(k + qp == q)(qp)(k)(q);
   return 1;
}

template <int n>
inline
double recoupling(Zn<n> const& q1, Zn<n> const& q3, Zn<n> const& q13,
                  Zn<n> const& q2, Zn<n> const& q, Zn<n> const& q23)
{
   DEBUG_PRECONDITION_EQUAL(q13, q1+q3);
   DEBUG_PRECONDITION_EQUAL(q23, q2+q3);
   DEBUG_PRECONDITION_EQUAL(q, q1+q23);
   return 1;
}

template <int n>
inline
double recoupling_12_3__13_2(Zn<n> const& q1, Zn<n> const& q2, Zn<n> const& q12,
                             Zn<n> const& q3, Zn<n> const& q, Zn<n> const& q13)
{
   DEBUG_PRECONDITION_EQUAL(q12, q1+q2);
   DEBUG_PRECONDITION_EQUAL(q13, q1+q3);
   DEBUG_PRECONDITION_EQUAL(q, q12+q3);
   return 1;
}

template <int n>
inline
Zn<n> adjoint(Zn<n> const& q)
{
   return Zn<n>(n-q.x);
}

template <int n>
inline
double adjoint_coefficient(Zn<n> const& qp, Zn<n> const& k, Zn<n> const& q)
{
   PRECONDITION_EQUAL(qp, q-k);
   return 1;
}

template <int n>
inline
double conj_phase(Zn<n> const& qp, Zn<n> const& k, Zn<n> const& q)
{
   PRECONDITION_EQUAL(qp, q-k);
   return 1;
}

template <int n>
inline
bool is_transform_target(Zn<n> const& q1, Zn<n> const& q2, Zn<n> const& q)
{
   return q == q1+q2;
}

template <int n>
inline
int num_transform_targets(Zn<n> const& q1, Zn<n> const& q2)
{
   return 1;
}

template <int n, typename OutIter>
inline
void transform_targets(Zn<n> const& q1, Zn<n> const& q2, OutIter Out)
{
   *Out++ = q1+q2;
}

template <int n>
inline
int num_inverse_transform_targets(Zn<n> const& q1, Zn<n> const& q)
{
   return 1;
}

template <int n, typename OutIter>
inline
void inverse_transform_targets(Zn<n> const& q1, Zn<n> const& q, OutIter Out)
{
   *Out++ = q1+q;
}

template <int n, typename OutIter>
inline
void enumerate_projections(Zn<n> const& q, OutIter Out)
{
   *Out++ = q;
}

template <int n>
inline
bool is_delta(Zn<n> const& q1, Zn<n> const& Q, Zn<n> const& P, Zn<n> const& q2)
{
   DEBUG_PRECONDITION_EQUAL(Q, P);
   return q1 == P+q2;
}

template <int n>
inline
Zn<n>
difference(Zn<n> const& q1, Zn<n> const& q2)
{
   return q1-q2;
}

template <int n>
inline
Zn<n> negate(Zn<n> const& p)
{
   return p;
}

template <int n>
inline
Zn<n> sum(Zn<n> const& q1, Zn<n> const& q2)
{
   return q1+q2;
}

template <int n>
inline
bool is_possible(Zn<n> const& q, Zn<n> const& p)
{
   return true;
}

template <int n>
inline
bool is_projection(Zn<n> const& q, Zn<n> const& p)
{
   return q == p;
}

template <int n>
inline
Zn<n> change(Zn<n> const& q, Zn<n> const& p)
{
   return q+p;
}

template <int n>
inline
Zn<n> heighest_weight(Zn<n> const& p)
{
   return p;
}

template <int n>
inline
double weight(Zn<n> const& p)
{
   return std::min(p.x.twice(), 2*n-p.x.twice());
}

template <int n>
inline
double delta_shift_coefficient(Zn<n> const& qp, Zn<n> const& k, Zn<n> const& q, Zn<n> const& Delta)
{
   PRECONDITION_EQUAL(qp, q+k);
   return 1;
}

template <int n>
inline
double casimir(Zn<n> const& s, int i)
{
   CHECK_EQUAL(i, 0);
   return s.x.to_double();
}

template <int n>
Zn<n>::Zn(std::string const& s)
  : x(convert_string<int>(s.begin(), s.end()))
{
}

template <int n>
std::string
Zn<n>::ToString() const
{
   return ConvertToString(x);
}

template <int n>
void Zn<n>::Register()
{
   RegisterStaticSymmetry<Zn<n>>();
}

namespace
{
   NiftyCounter::nifty_counter<Zn<2>::Register> ZnRegistrator2;
   NiftyCounter::nifty_counter<Zn<3>::Register> ZnRegistrator3;
   NiftyCounter::nifty_counter<Zn<4>::Register> ZnRegistrator4;
   NiftyCounter::nifty_counter<Zn<5>::Register> ZnRegistrator5;
   NiftyCounter::nifty_counter<Zn<6>::Register> ZnRegistrator6;
   NiftyCounter::nifty_counter<Zn<7>::Register> ZnRegistrator7;
   NiftyCounter::nifty_counter<Zn<8>::Register> ZnRegistrator8;
   NiftyCounter::nifty_counter<Zn<9>::Register> ZnRegistrator9;
   NiftyCounter::nifty_counter<Zn<10>::Register> ZnRegistrator10;
   NiftyCounter::nifty_counter<Zn<11>::Register> ZnRegistrator11;
   NiftyCounter::nifty_counter<Zn<12>::Register> ZnRegistrator12;
}

} // namespace QuantumNumbers

#endif
