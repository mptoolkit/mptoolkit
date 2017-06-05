// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/su2.h
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

  The SU(2) quantum number registers an object that allows the name "SU(2)"
  in the SymmetryList.  eg, MyQuantumNumbers = SymmetryList("S:SU(2)")

  We have two braid group representations defined, bosonic and fermionic.
*/

// definitions:
// integral: all dimensions are integers
// pointed: all dimensions are 1


#if !defined(MPTOOKLIT_QUANTUMNUMBERS_SU2_H)
#define MPTOOKLIT_QUANTUMNUMBERS_SU2_H

// definitions:
// integral: all dimensions are integers
// pointed: all dimensions are 1

#include "common/types.h"
#include "common/niftycounter.h"
#include "common/halfint.h"
#include "quantumnumber.h"
#include "coupling_su2.h"
#include <iostream>
#include <type_traits>

namespace QuantumNumbers
{

// detect if class T has a .degree(T::value_type) member function
template <typename T>
struct has_degree_function
{
   template <typename C> static constexpr
   decltype(std::declval<C>().degree(std::declval<typename T::value_type>()), bool())
      test(int)
   {
      return true;
   }

   template <typename C> static constexpr bool test(...)
   {
      return false;
   }
   static constexpr bool value = test<T>(int());
   using type = std::integral_constant<bool, value>;
};

// helper class for a multiplicity-free Lie group

template <typename Derived>
class StaticLieGroup_MF
{
   public:
      using value_type = Derived::value_type;
      using is_multiplicity_free = std::true_type;
      using is_integral = std::true_type;
      using is_pointed = std::true_type;
      using is_finite = std::false_type;
      using is_real = std::true_type;

      static real qdim(value_type v)
      {
         return Derived::degree(v);
      }
};

// helper class for a multiplicity-free finite group
template <typename Derived>
class StaticFiniteGroup_MF

{
   public:
      using value_type = Derived::value_type;
      using is_multiplicity_free = std::true_type;
      using is_integral = std::true_type;
      using is_pointed = std::true_type;
      using is_finite = std::true_type;
      using is_real = std::true_type;

      static real qdim(value_type v)
      {
         return Derived::degree(v);
      }
};


>>>>>>> Work in progress on new QN library
// the identity rep of the braid group
template <typename Symmetry>
class Bosonic
{
   public:
      using value_type = Symmetry::value_type;
      using is_real = std::true_type;
      using is_symmetric = std::true_type;     // alias for is_real

      explicit Bosonic(Symmetry const& s_, std::string const& Name) : symmetry(s_)
      {
         CHECK_EQUAL(Name, "bosonic");
      }

      static std::string name() { return "bosonic"; }

      static bool is_valid(Symmetry const& s, std::string const& Name)
      {
         return Name == "bosonic";
      }

      static real r_matrix(value_type a, value_type b, value_type c)
      {
         DEBUG_CHECK(symmetry.is_transform_target(a,b,c));
         return 1;
      }

   private:
      Symmetry symmetry;

};

template <typename Symmetry>
class FermionicSpin
{
   public:
      using is_real = std::true_type;
      using is_symmetric = std::true_type;     // alias for is_real
      using value_type = Symmetry::value_type;

      static_assert(std::is_same<value_type, half_int>::value);

      explicit FermionicSpin(Symmetry const& s_, std::string const& Name) : symmetry(s_)
      {
         CHECK_EQUAL(Name, "FermionicSpin");
      }

      static std::string name() { return "FermionicSpin"; }

      static bool is_valid(Symmetry const& s, std::string const& Name)
      {
         return Name == "FermionicSpin";
      }

      static real r_matrix(half_int a, half_int b, half_int c)
      {
         DEBUG_CHECK(symmetry.is_transform_target(a,b,c));
         return (is_integral(a) || is_integral(b)) ? 1.0 : -1.0;
      }

   private:
      Symmetry symmetry;

};

template <typename T>
class FermionicParticleImpl;

template <>
class FermionicParticleImpl<int>
{
   public:
      static real r_matrix(int a, int b, int)
      {
         return ((a%2 == 1) && (b%2 == 1)) ? -1.0 : 1.0;
      }
};

template <>
class FermionicParticleImpl<half_int>
{
   public:
      static real r_matrix(half_int a, half_int b, half_int)
      {
         DEBUG_CHECK(is_integral(a));
         DEBUG_CHECK(is_integral(b));
         return ((a.to_int_assert()%2 == 1) && (b.to_int_assert()%2 == 1)) ? -1.0 : 1.0;
      }
};

template <typename Symmetry>
class FermionicParticle : public FermionicParticleImpl<typename Symmetry::value_type>
{
   public:
      using is_real = std::true_type;
      using is_symmetric = std::true_type;     // alias for is_real
      using value_type = Symmetry::value_type;

      static_assert(std::is_same<value_type, half_int>::value
                    || std::is_same<value_type, int>::value);

      explicit FermionicParticle(Symmetry const& s_, std::string const& Name) : symmetry(s_)
      {
         CHECK_EQUAL(Name, "FermionicParticle");
      }

      static std::string name() { return "FermionicParticle"; }

      static bool is_valid(Symmetry const& s, std::string const& Name)
      {
         return Name == "FermionicParticle";
      }

      static real r_matrix(value_type a, value_type b, value_type c)
      {
         DEBUG_CHECK(symmetry.is_transform_target(a,b,c));
         return FermionicParticleImpl<value_type>::r_matrix(a,b,c);
      }

   private:
      Symmetry symmetry;

};

class SU2 : public StaticLieGroup_MF<SU2>
>>>>>>> sketch for how braid reps would work
{
   public:
      using value_type = Derived::value_type;
      using is_multiplicity_free = std::true_type;
      using is_integral = std::true_type;
      using is_pointed = std::true_type;
      using is_finite = std::true_type;
      using is_real = std::true_type;

      static real qdim(value_type v)
      {
         return Derived::degree(v);
      }
};


// the identity rep of the braid group
template <typename Symmetry>
class Bosonic
{
   public:
      using value_type = Symmetry::value_type;
      using is_real = std::true_type;
      using is_symmetric = std::true_type;     // alias for is_real

      explicit Bosonic(Symmetry const& s_, std::string const& Name) : symmetry(s_)
      {
         CHECK_EQUAL(Name, "bosonic");
      }

      static std::string name() { return "bosonic"; }

      static bool is_valid(Symmetry const& s, std::string const& Name)
      {
         return Name == "bosonic";
      }

      static real r_matrix(value_type a, value_type b, value_type c)
      {
         DEBUG_CHECK(symmetry.is_transform_target(a,b,c));
         return 1;
      }

   private:
      Symmetry symmetry;

};

template <typename Symmetry>
class FermionicSpin
{
   public:
<<<<<<< HEAD
      using is_real = std::true_type;
      using is_symmetric = std::true_type;     // alias for is_real
      using value_type = Symmetry::value_type;

      static_assert(std::is_same<value_type, half_int>::value);

      explicit FermionicSpin(Symmetry const& s_, std::string const& Name) : symmetry(s_)
      {
         CHECK_EQUAL(Name, "FermionicSpin");
      }

      static std::string name() { return "FermionicSpin"; }

      static bool is_valid(Symmetry const& s, std::string const& Name)
      {
         return Name == "FermionicSpin";
      }

      static real r_matrix(half_int a, half_int b, half_int c)
      {
         DEBUG_CHECK(symmetry.is_transform_target(a,b,c));
         return (is_integral(a) || is_integral(b)) ? 1.0 : -1.0;
      }

   private:
      Symmetry symmetry;

};

template <typename T>
class FermionicParticleImpl;

template <>
class FermionicParticleImpl<int>
{
   public:
      static real r_matrix(int a, int b, int)
      {
         return ((a%2 == 1) && (b%2 == 1)) ? -1.0 : 1.0;
      }
};

template <>
class FermionicParticleImpl<half_int>
{
   public:
      static real r_matrix(half_int a, half_int b, half_int)
      {
         DEBUG_CHECK(is_integral(a));
         DEBUG_CHECK(is_integral(b));
         return ((a.to_int_assert()%2 == 1) && (b.to_int_assert()%2 == 1)) ? -1.0 : 1.0;
      }
};

template <typename Symmetry>
class FermionicParticle : public FermionicParticleImpl<typename Symmetry::value_type>
{
   public:
      using is_real = std::true_type;
      using is_symmetric = std::true_type;     // alias for is_real
      using value_type = Symmetry::value_type;

      static_assert(std::is_same<value_type, half_int>::value
                    || std::is_same<value_type, int>::value);

      explicit FermionicParticle(Symmetry const& s_, std::string const& Name) : symmetry(s_)
      {
         CHECK_EQUAL(Name, "FermionicParticle");
      }

      static std::string name() { return "FermionicParticle"; }

      static bool is_valid(Symmetry const& s, std::string const& Name)
      {
         return Name == "FermionicParticle";
      }

      static real r_matrix(value_type a, value_type b, value_type c)
      {
         DEBUG_CHECK(symmetry.is_transform_target(a,b,c));
         return FermionicParticleImpl<value_type>::r_matrix(a,b,c);
      }

   private:
      Symmetry symmetry;

};

class SU2 : public StaticLieGroup_MF<SU2>
{
   public:
      using value_type      = half_int;
      using difference_type = half_int;

      // conversion to/from string representation
      static std::string as_string(half_int j);
      static half_int from_string(std::string const& s);

      // converstion to/from int storage (not to be confused with the integer enumeration!)
      static half_int from_int(int n)
      {
         return half_int::from_twice(n);
      }

      static int to_int(half_int j)
      {
         return j.twice();
      }

      static std::string static_name() { return "SU(2)"; }

      // 'friendly' ordering of quantum numbers.  The identity must
      // compare as less than all other non-identity reps.
      static bool less(half_int j1, half_int j2)
      {
         // order is 0,1/2,-1/2,1,-1,3/2,-3/2,...
         return abs(j1) < abs(j2) || (abs(j1) == abs(j2) && j2 < j1);
      }

      // 'size' is the number of members of the category.  For an infinite category,
      // return 0
      static int size()
      {
         return 0;
      }

      // enumerate the possible quantum numbers from a non-negative integer.
      // This must be consistent with the ordering from the 'less' function:
      // for integers i,j, we have:
      // i < j iff less(enumerate(i), enumerate(j))
      // (this implies that the scalar (identity) quantum number is n=0
      static half_int enumerate(int n)
      {
         if (n == 0)
            return 0;
         if (n % 2 == 1)
         {
            return half_int::from_twice(n/2);
         }
         return -half_int::from_twice(n/2);
      }

      using std::tuple<FermionicSpin<SU2>, Bosonic<SU2>> braid_rep_types;

      using default_braid_type = Bosonic;

      static half_int from_int(int n)
      {
         return half_int::from_twice(n);
      }

      static int to_int(half_int j)
      {
         return j.twice();
      }

      static std::string static_name() { return "SU(2)"; }

      // 'friendly' ordering of quantum numbers.  The identity must
      // compare as less than all other non-identity reps.
      static bool less(half_int j1, half_int j2)
      {
         // order is 0,1/2,-1/2,1,-1,3/2,-3/2,...
         return abs(j1) < abs(j2) || (abs(j1) == abs(j2) && j2 < j1);
      }

      // 'size' is the number of members of the category.  For an infinite category,
      // return 0
      static int size()
      {
         return 0;
      }

      // enumerate the possible quantum numbers from a non-negative integer.
      // This must be consistent with the ordering from the 'less' function:
      // for integers i,j, we have:
      // i < j iff less(enumerate(i), enumerate(j))
      // (this implies that the scalar (identity) quantum number is n=0
      static half_int enumerate(int n)
      {
         if (n == 0)
            return 0;
         if (n % 2 == 1)
         {
            return half_int::from_twice(n/2);
         }
         return -half_int::from_twice(n/2);
      }

      using std::tuple<FermionicSpin<SU2>, Bosonic<SU2>> braid_rep_types;

      using default_braid_type = Bosonic;

      static default_braid_type default_braid_rep()
      {
         return default_braid_type();
      }

      // Casmir invariant operators (Lie algebra terminology).
      // For finite groups, these are the expectation values of the centre
      // of the enveloping algebra (a linearly independent set that spans)
      static int num_casimir() { return 1; }

      static std::string casimir_name(std::string const& QName, int)
      { return QName + "^2"; }

      static real casimir(half_int j, int n)
      {
         DEBUG_CHECK_EQUAL(n, 0);
         return j.to_real() * (j.to_real() + 1);
      }

      // true if this is the scalar (identity) quantum number
      static bool is_scalar(half_int j)
      {
         return j == 0;
      }

      // degree of the representation
      static int degree(half_int j)
      {
         return j.twice() + 1;
      }

      static half_int adjoint(half_int j)
      {
         return j;
      }

      // cross product is defined only for vector operators
      static bool cross_product_exists(half_int j1, half_int j2)
      {
         return j1 == 1 && j2 == 1;
      }

      static half_int cross_product_transforms_as(half_int j1, half_int j2)
      {
         DEBUG_CHECK_EQUAL(j1, 1);
         DEBUG_CHECK_EQUAL(j2, 1);
         return j1;
      }

      static complex cross_product_factor(half_int j1, half_int j2)
      {
         DEBUG_CHECK_EQUAL(j1, 1);
         DEBUG_CHECK_EQUAL(j2, 1);
         return complex(0.0, math_const::sqrt_2);
      }

      static real coupling_3j_phase(half_int j1, half_int j2, half_int j)
      {
         return minus1pow(j1+j2-j);
      }

      static real coupling_6j(half_int j1, half_int j2, half_int j3,
                              half_int j4, half_int j5, half_int j6)
      {
         return CouplingSU2::Coupling6j(j1, j2, j3, j4, j5, j6);
      }

      static real coupling_9j(half_int j1 , half_int j2 , half_int j12,
                              half_int j3 , half_int j4 , half_int j34,
                              half_int j13, half_int j24, half_int j  )
      {
         return CouplingSU2::Coupling9j(j1, j2, j12, j3, j4, j34, j12, j24, j);
      }

      static bool is_transform_target(half_int j1, half_int j2, half_int j)
      {
         return (j.twice() <= j1.twice() + j2.twice())
            && (j.twice() >= std::abs(j1.twice() - j2.twice()));
      }

      static int num_transform_targets(half_int j1, half_int j2)
      {
         return std::min(j1.twice(), j2.twice()) + 1;
      }

      template <typename OutIter>
      static OutIter transform_targets(half_int j1, half_int j2, OutIter Out)
      {
         for (half_int j = abs(j1 - j2); j <= j1 + j2; ++j)
         {
            *Out++ = j;
         }
         return Out;
      }

      static difference_type difference(half_int j1, half_int j2)
      {
         return j1 - j2;
      }

      static bool is_possible_change(half_int j, difference_type delta)
      {
         return j + delta >= 0;
      }

      static half_int change_by(half_int j, difference_type delta)
      {
         DEBUG_CHECK(J+delta >= 0);
         return j + delta;
      }

      static real weight(difference_type x)
      {
         return std::abs(x.to_real());
      }

      SU2() : j(0) {}
      SU2(half_int j_) : j(j_) {}
      SU2(int j_) : j(j_) {}
      SU2(real j_) : j(j_) {}
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
real trace(SU2 const& q)
{
   return q.j.twice() + 1;
}

inline
real identity(SU2 const& q)
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
complex cross_product_factor(SU2 const& q1, SU2 const& q2)
{
   DEBUG_CHECK_EQUAL(q1.j, 1);
   DEBUG_CHECK_EQUAL(q2.j, 1);
   return complex(0.0, math_const::sqrt_2);
}

inline
real clebsch_gordan(SU2 const& q1, SU2 const& q2, SU2 const& q,
                      Sz const& m1,  Sz const& m2,  Sz const& m)
{
   PRECONDITION(is_triangle(q1.j, q2.j, q.j));
   using std::sqrt;

   // to get the convention of Varshalovich et al insert factor
   return
#if defined(PRESERVE_NORM)
      (minus1pow(q2.j.twice()) / sqrt(q1.j.twice()+1)) *
#endif
      ClebschGordan(q1.j, m1.m, q2.j, m2.m, q.j, m.m);
}

inline
real product_coefficient(SU2 const& k1, SU2 const& k2, SU2 const& k,
                         SU2 const& qp, SU2 const& q, SU2 const& qpp)
{
   using std::sqrt;
   DEBUG_PRECONDITION(is_triangle(k1.j, k2.j, k.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k1.j, qpp.j));
   DEBUG_PRECONDITION(is_triangle(qpp.j, k2.j, q.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k.j, q.j));

   return minus1pow(to_int(qp.j + q.j + k.j)) * sqrt(degree(qpp)) * sqrt(degree(k)) *
      Coupling6j(qp.j, k1.j, qpp.j, k2.j, q.j, k.j);
}



inline
real inverse_product_coefficient(SU2 const& k1, SU2 const& k2, SU2 const& k,
                                 SU2 const& qp, SU2 const& q, SU2 const& qpp)
{
   using std::sqrt;
   DEBUG_PRECONDITION(is_triangle(k1.j, k2.j, k.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k1.j, qpp.j));
   DEBUG_PRECONDITION(is_triangle(qpp.j, k2.j, q.j));
   DEBUG_PRECONDITION(is_triangle(qp.j, k.j, q.j));

   return minus1pow(to_int(k1.j + k2.j - k.j)) * sqrt(degree(qpp)) * sqrt(degree(k)) *
      Racah(qp.j, k1.j, q.j, k2.j, qpp.j, k.j);
}

inline
real tensor_coefficient(SU2 const& j1,  SU2 const& j2,  SU2 const& j12,
                        SU2 const& j3,  SU2 const& j4,  SU2 const& j34,
                        SU2 const& j13, SU2 const& j24, SU2 const& j)
{
   using std::sqrt;
   DEBUG_PRECONDITION(is_triangle(j1.j,  j2.j, j12.j))(j1.j)(j2.j)(j12.j);
   DEBUG_PRECONDITION(is_triangle(j3.j, j4.j, j34.j))(j3.j)(j4.j)(j34.j);
   DEBUG_PRECONDITION(is_triangle(j13.j, j24.j, j.j))(j13.j)(j24.j)(j.j);
   DEBUG_PRECONDITION(is_triangle(j1.j,  j3.j, j13.j))(j1.j)(j3.j)(j13.j);
   DEBUG_PRECONDITION(is_triangle(j2.j, j4.j, j24.j))(j2.j)(j4.j)(j24.j);
   DEBUG_PRECONDITION(is_triangle(j12.j, j34.j, j.j))(j12.j)(j34.j)(j.j);

   real f1 = sqrt(degree(j12));
   real f2 = sqrt(degree(j13));
   real f3 = sqrt(degree(j24));
   real f4 = sqrt(degree(j34));

   real Result =  f1*f2*f3*f4*Coupling9j(j1.j,  j2.j,  j12.j,
                                           j3.j,  j4.j,  j34.j,
                                           j13.j, j24.j, j.j);
   return Result;
}

inline
real inverse_tensor_coefficient(SU2 const& j1,  SU2 const& j2,  SU2 const& j12,
                                  SU2 const& j3,  SU2 const& j4,  SU2 const& j34,
                                  SU2 const& j13, SU2 const& j24, SU2 const& j)
{
   using std::sqrt;
   DEBUG_PRECONDITION(is_triangle(j1.j,  j2.j, j12.j))(j1.j)(j2.j)(j12.j);
   DEBUG_PRECONDITION(is_triangle(j3.j, j4.j, j34.j))(j3.j)(j4.j)(j34.j);
   DEBUG_PRECONDITION(is_triangle(j13.j, j24.j, j.j))(j13.j)(j24.j)(j.j);
   DEBUG_PRECONDITION(is_triangle(j1.j,  j3.j, j13.j))(j1.j)(j3.j)(j13.j);
   DEBUG_PRECONDITION(is_triangle(j2.j, j4.j, j24.j))(j2.j)(j4.j)(j24.j);
   DEBUG_PRECONDITION(is_triangle(j12.j, j34.j, j.j))(j12.j)(j34.j)(j.j);

   real f1 = sqrt(degree(j12));
   real f2 = sqrt(degree(j13));
   real f3 = sqrt(degree(j24));
   real f4 = sqrt(degree(j34));

   real denominator = f1*f2*f3*f4;
   real numerator = degree(j3)*degree(j4)*degree(j12)*degree(j);

   real Result =  (numerator/denominator)*Coupling9j(j1.j,  j2.j,  j12.j,
                                                       j3.j,  j4.j,  j34.j,
                                                       j13.j, j24.j, j.j);
   return Result;
}

inline
real recoupling(SU2 const& q1, SU2 const& q2, SU2 const& q12,
                SU2 const& q3, SU2 const& q, SU2 const& q23)
{
   using std::sqrt;
   return sqrt(degree(q12)) * sqrt(degree(q23)) *
      Racah(q1.j, q2.j, q.j, q3.j, q12.j, q23.j);
}

inline
real recoupling_12_3__13_2(SU2 const& q1, SU2 const& q2,
                             SU2 const& q12,
                             SU2 const& q3, SU2 const& q,
                             SU2 const& q13)
{
   using std::sqrt;
   return minus1pow(to_int(q1.j+q.j+q12.j+q13.j)) * sqrt(degree(q13)) * sqrt(degree(q12))
      * Racah(q1.j, q2.j, q3.j, q.j, q12.j, q13.j);
}

inline
SU2 adjoint(SU2 const& q)
{
   return q;
}

inline
real adjoint_coefficient(SU2 const& qp, SU2 const& k, SU2 const& q)
{
   using std::sqrt;
   return minus1pow(to_int(q.j + k.j - qp.j)) * sqrt(degree(q)/real(degree(qp)));
}

inline
real conj_phase(SU2 const& qp, SU2 const& k, SU2 const& q)
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
real weight(Sz const& p)
{
   using std::abs;
   return abs(p.m.to_real());
}

inline
real delta_shift_coefficient(SU2 const& qp, SU2 const& k, SU2 const& q, SU2 const& Delta)
{
   using std::sqrt;
   using std::abs;
   real Result = minus1pow(to_int(q.j - k.j + qp.j))
      * sqrt(real(q.j.twice()+Delta.j.twice()+1)*real(qp.j.twice()+1))
      * Racah(q.j, qp.j, q.j+Delta.j, qp.j+Delta.j, k.j, Delta.j);

   real ResultCheck = tensor_coefficient(Delta, q, SU2(Delta.j+q.j),
                                           SU2(0), k, k,
                                           Delta, qp, SU2(Delta.j+qp.j));

   CHECK(abs(Result-ResultCheck) < 1E-12)(Result)(ResultCheck);

   return Result;
}

inline
real casimir(SU2 const& s, int n)
{
   CHECK_EQUAL(n, 0);
   return s.j.to_real() * (s.j.to_real() + 1);
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
