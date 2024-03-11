// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/su2.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ian@qusim.net>
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

  The SU(2) quantum number registers an object that allows the name "SU(2)"
  in the SymmetryList.  eg, MyQuantumNumbers = SymmetryList("S:SU(2)")

  This header simply has a nifty counter that ensures that the static object is initialized.
*/

#if !defined(SU2_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY)
#define SU2_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY

#include "common/niftycounter.h"
#include "common/halfint.h"
#include "quantumnumber.h"
#include "coupling.h"
#include <iostream>

//#define PRESERVE_NORM

namespace QuantumNumbers
{

class Sz;

class SU2
{
   public:
      typedef Sz                              ProjectionType;
      typedef StaticQuantumNumberFactory<SU2> FactoryType;

      SU2() : j(0) {}
      SU2(half_int j_) : j(j_) {}
      SU2(int j_) : j(j_) {}
      SU2(double j_) : j(j_) {}
      explicit SU2(int const* InIter);
      explicit SU2(std::string const& s);

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { *OutIter = j.twice(); return OutIter+1; }

      static char const* Type() { return "SU(2)"; }

      static int Size() { return 1; }

      static int num_casimir() { return 1; }
      static std::string casimir_name(std::string const& QName, int)
      { return QName + "^2"; }

      half_int j;

      // Registation is automatic via a nifty counter
      static void Register();
};

inline
std::ostream& operator<<(std::ostream& out, SU2 const& s)
{
   return out << s.j;
}

class Sz
{
   public:
      typedef SU2 QuantumNumberType;

      Sz(half_int m_) : m(m_) {}
      Sz(int m_) : m(m_) {}
      Sz(double m_) : m(m_) {}
      explicit Sz(int const* InIter);
      explicit Sz(std::string const& s);

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { *OutIter = m.twice(); return ++OutIter; }

      static char const* Type() { return "SU(2)"; }

      static char const* Suffix() { return "_m"; }

      static int Size() { return 1; }

      half_int m;
};

inline
std::ostream& operator<<(std::ostream& out, Sz const& s)
{
   return out << s.m;
}
//
// inlines
//

inline
SU2::SU2(int const* InIter)
  : j(from_twice(*InIter))
{
}

inline
std::string SU2::ToString() const
{
   return ConvertToString(j);
}

inline
int degree(SU2 const& q)
{
   return q.j.twice() + 1;
}

inline
double trace(SU2 const& q)
{
#if defined(PRESERVE_NORM)
   return 1;
#else
   return q.j.twice() + 1;
#endif
}

inline
double identity(SU2 const& q)
{
#if defined(PRESERVE_NORM)
   return q.j.twice() + 1;
#else
   return 1;
#endif
}

inline
int multiplicity(SU2 const& q1, SU2 const& q2, SU2 const& q)
{
   return 1;
}

// we only know how to define the cross product for vectors
inline
bool cross_product_exists(SU2 const& q1, SU2 const& q2)
{
   return q1.j == 1 && q2.j == 1;
}

inline
SU2 cross_product_transforms_as(SU2 const& q1, SU2 const& q2)
{
   DEBUG_CHECK_EQUAL(q1.j, 1);
   DEBUG_CHECK_EQUAL(q2.j, 1);
   return q1;
}

inline
std::complex<double> cross_product_factor(SU2 const& q1, SU2 const& q2)
{
   DEBUG_CHECK_EQUAL(q1.j, 1);
   DEBUG_CHECK_EQUAL(q2.j, 1);
   return std::complex<double>(0.0, std::sqrt(2.0));
}

inline
double clebsch_gordan(SU2 const& q1, SU2 const& q2, SU2 const& q,
                      Sz const& m1,  Sz const& m2,  Sz const& m)
{
   PRECONDITION(is_triangle(q1.j, q2.j, q.j));

   // to get the convention of Varshalovich et al insert factor
   return
#if defined(PRESERVE_NORM)
      (minus1pow(q2.j.twice()) / sqrt(q1.j.twice()+1)) *
#endif
      ClebschGordan(q1.j, m1.m, q2.j, m2.m, q.j, m.m);
}

inline
double product_coefficient(SU2 const& k1, SU2 const& k2, SU2 const& k,
                          SU2 const& qp, SU2 const& q, SU2 const& qpp)
{
   DEBUG_PRECONDITION(is_triangle(k1.j, k2.j, k.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k1.j, qpp.j));
   DEBUG_PRECONDITION(is_triangle(qpp.j, k2.j, q.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k.j, q.j));

   return minus1pow(to_int(k1.j + k2.j - k.j)) * sqrt(degree(qpp)) * sqrt(degree(k)) *
     Racah(qp.j, k1.j, q.j, k2.j, qpp.j, k.j);
}

inline
double inverse_product_coefficient(SU2 const& k1, SU2 const& k2, SU2 const& k,
                                   SU2 const& qp, SU2 const& q, SU2 const& qpp)
{
   DEBUG_PRECONDITION(is_triangle(k1.j, k2.j, k.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k1.j, qpp.j));
   DEBUG_PRECONDITION(is_triangle(qpp.j, k2.j, q.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k.j, q.j));

   return minus1pow(to_int(k1.j + k2.j - k.j)) * sqrt(degree(qpp)) * sqrt(degree(k)) *
      Racah(qp.j, k1.j, q.j, k2.j, qpp.j, k.j);
}

inline
double tensor_coefficient(SU2 const& j1,  SU2 const& j2,  SU2 const& j12,
                          SU2 const& j3,  SU2 const& j4,  SU2 const& j34,
                          SU2 const& j13, SU2 const& j24, SU2 const& j)
{
   DEBUG_PRECONDITION(is_triangle(j1.j,  j2.j, j12.j))(j1.j)(j2.j)(j12.j);
   DEBUG_PRECONDITION(is_triangle(j3.j, j4.j, j34.j))(j3.j)(j4.j)(j34.j);
   DEBUG_PRECONDITION(is_triangle(j13.j, j24.j, j.j))(j13.j)(j24.j)(j.j);
   DEBUG_PRECONDITION(is_triangle(j1.j,  j3.j, j13.j))(j1.j)(j3.j)(j13.j);
   DEBUG_PRECONDITION(is_triangle(j2.j, j4.j, j24.j))(j2.j)(j4.j)(j24.j);
   DEBUG_PRECONDITION(is_triangle(j12.j, j34.j, j.j))(j12.j)(j34.j)(j.j);

   double f1 = sqrt(degree(j12));
   double f2 = sqrt(degree(j13));
   double f3 = sqrt(degree(j24));
   double f4 = sqrt(degree(j34));

   double Result =  f1*f2*f3*f4*Coupling9j(j1.j,  j2.j,  j12.j,
                                           j3.j,  j4.j,  j34.j,
                                           j13.j, j24.j, j.j);
   return Result;
}

inline
double inverse_tensor_coefficient(SU2 const& j1,  SU2 const& j2,  SU2 const& j12,
                                  SU2 const& j3,  SU2 const& j4,  SU2 const& j34,
                                  SU2 const& j13, SU2 const& j24, SU2 const& j)
{
   DEBUG_PRECONDITION(is_triangle(j1.j,  j2.j, j12.j))(j1.j)(j2.j)(j12.j);
   DEBUG_PRECONDITION(is_triangle(j3.j, j4.j, j34.j))(j3.j)(j4.j)(j34.j);
   DEBUG_PRECONDITION(is_triangle(j13.j, j24.j, j.j))(j13.j)(j24.j)(j.j);
   DEBUG_PRECONDITION(is_triangle(j1.j,  j3.j, j13.j))(j1.j)(j3.j)(j13.j);
   DEBUG_PRECONDITION(is_triangle(j2.j, j4.j, j24.j))(j2.j)(j4.j)(j24.j);
   DEBUG_PRECONDITION(is_triangle(j12.j, j34.j, j.j))(j12.j)(j34.j)(j.j);

   double f1 = sqrt(degree(j12));
   double f2 = sqrt(degree(j13));
   double f3 = sqrt(degree(j24));
   double f4 = sqrt(degree(j34));

   double denominator = f1*f2*f3*f4;
   double numerator = degree(j3)*degree(j4)*degree(j12)*degree(j);

   double Result =  (numerator/denominator)*Coupling9j(j1.j,  j2.j,  j12.j,
                                                       j3.j,  j4.j,  j34.j,
                                                       j13.j, j24.j, j.j);
   return Result;
}

inline
double recoupling(SU2 const& q1, SU2 const& q2, SU2 const& q12,
                  SU2 const& q3, SU2 const& q, SU2 const& q23)
{
   return sqrt(degree(q12)) * sqrt(degree(q23)) *
      Racah(q1.j, q2.j, q.j, q3.j, q12.j, q23.j);
}

inline
double recoupling_12_3__13_2(SU2 const& q1, SU2 const& q2,
                             SU2 const& q12,
                             SU2 const& q3, SU2 const& q,
                             SU2 const& q13)
{
   return minus1pow(to_int(q1.j+q.j+q12.j+q13.j)) * sqrt(degree(q13)) * sqrt(degree(q12))
      * Racah(q1.j, q2.j, q3.j, q.j, q12.j, q13.j);
}

inline
SU2 adjoint(SU2 const& q)
{
  return q;
}

inline
double adjoint_coefficient(SU2 const& qp, SU2 const& k, SU2 const& q)
{
   return minus1pow(to_int(q.j + k.j - qp.j)) * sqrt(degree(q)/double(degree(qp)));
}

inline
double conj_phase(SU2 const& qp, SU2 const& k, SU2 const& q)
{
   return minus1pow(to_int(q.j + k.j - qp.j));
}

inline
bool is_transform_target(SU2 const& q1, SU2 const& q2, SU2 const& q)
{
   return is_triangle(q1.j, q2.j, q.j);
}

inline
int num_transform_targets(SU2 const& q1, SU2 const& q2)
{
   return std::min(q1.j.twice(), q2.j.twice()) + 1;
}

template <typename OutIter>
inline
void transform_targets(SU2 const& q1, SU2 const& q2, OutIter Out)
{
   for (half_int j = abs(q1.j - q2.j); j <= q1.j + q2.j; ++j)
   {
      *Out++ = j;
   }
}

inline
int num_inverse_transform_targets(SU2 const& q1, SU2 const& q)
{
   return std::min(q1.j.twice(), q.j.twice()) + 1;
}

template <typename OutIter>
inline
void inverse_transform_targets(SU2 const& q1, SU2 const& q, OutIter Out)
{
   for (half_int j = abs(q1.j - q.j); j <= q1.j + q.j; ++j)
   {
      *Out++ = j;
   }
}

template <typename OutIter>
inline
void enumerate_projections(SU2 const& q, OutIter Out)
{
  //  std::cout << "Enumerating projections of " << q << std::endl;
   for (half_int m = -q.j; m <= q.j; ++m)
   {
      *Out++ = m;
   }
}

inline
bool is_delta(SU2 const& q1, SU2 const& Q, Sz const& P, SU2 const& q2)
{
   // note: Q is unused here.
   return q1.j == P.m + q2.j;
}

inline
Sz difference(SU2 const& q1, SU2 const& q2)
{
   return Sz(q1.j-q2.j);
}

inline
Sz negate(Sz const& p)
{
   return Sz(-p.m);
}

inline
Sz sum(Sz const& a, Sz const& b)
{
   return Sz(a.m + b.m);
}

inline
bool is_possible(SU2 const& q, Sz const& p)
{
   return q.j+p.m >= 0;
}

inline
bool is_projection(SU2 const& q, Sz const& p)
{
   return -q.j <= p.m && p.m <= q.j;
}

inline
SU2 change(SU2 const& q, Sz const& p)
{
   DEBUG_CHECK(q.j+p.m >= 0)(q.j)(p.m);
   return SU2(q.j+p.m);
}

inline
SU2 heighest_weight(Sz const& z)
{
   return SU2(abs(z.m));
}

inline
double weight(Sz const& p)
{
   return fabs(p.m.to_double());
}

inline
double delta_shift_coefficient(SU2 const& qp, SU2 const& k, SU2 const& q, SU2 const& Delta)
{
   double Result = minus1pow(to_int(q.j - k.j + qp.j))
      * sqrt(double(q.j.twice()+Delta.j.twice()+1)*double(qp.j.twice()+1))
      * Racah(q.j, qp.j, q.j+Delta.j, qp.j+Delta.j, k.j, Delta.j);

   double ResultCheck = tensor_coefficient(Delta, q, SU2(Delta.j+q.j),
                                          SU2(0), k, k,
                                           Delta, qp, SU2(Delta.j+qp.j));

   CHECK(fabs(Result-ResultCheck) < 1E-12)(Result)(ResultCheck);

   return Result;
}

inline
double casimir(SU2 const& s, int n)
{
   CHECK_EQUAL(n, 0);
   return s.j.to_double() * (s.j.to_double() + 1);
}

//
// Sz
//

inline
Sz::Sz(int const* InIter)
  : m(from_twice(*InIter))
{
}

inline
std::string Sz::ToString() const
{
   return ConvertToString(m);
}

namespace
{
   NiftyCounter::nifty_counter<SU2::Register> SU2Counter;
}

} // namespace QuantumNumbers

#endif
