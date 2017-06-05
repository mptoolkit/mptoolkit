// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/fibonacci.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  Fibonacci anyons, equivalent to SU(2)_3.  Quantum number labels are '1' and 't'
*/

#if !defined(MPTOOKLIT_QUANTUMNUMBERS_FIBONACCI_H)
#define MPTOOKLIT_QUANTUMNUMBERS_FIBONACCI_H

#include "common/types.h"
#include "common/niftycounter.h"
#include "common/math_const.h"
#include "quantumnumber.h"
#include <iostream>
#include <type_traits>

namespace QuantumNumbers
{

class Fibonacci;

class FibonacciBraid
{
   public:
      using is_real    = std::false_type;
      using value_type = int;

      FibonacciBraid() {}

      explicit FibonacciBraid(Fibonacci const&, std::string const& Name)
      {
         CHECK_EQUAL(Name, "fibonacci");
      }

      static std::string name() { return "fibonacci"; }

      static bool is_valid(Fibonacci const&, std::string const& Name)
      {
         retrun Name == "fibonacci";
      }

      using std::exp;
      static constexpr complex exp_3_5 = exp(complex(0,3*math_const::pi/5));
      static constexpr complex exp_4_5 = exp(complex(0,-4*math_const::pi/5));

      static complex r_matrix(int a, int b, int c)
      {
         DEBUG_CHECK((a+b == c) || (a+b > 1));
         if (a == 1 && b == 1)
         {
            return c ? exp_3_5 : exp_4_5;
         }
         return 1;
      }
};

class Fibonacci
{
   public:
      using value_type      = int;
      using difference_type = int;

      Fibonacci() {}

      static std::string static_name() { return "fibonacci"; }

      // conversion to/from string representation
      static std::string as_string(int j);
      static int from_string(std::string const& s);

      // converstion to/from int storage (not to be confused with the integer enumeration!)
      static int from_int(int n)
      {
         return n;
      }

      static int to_int(int j)
      {
         return j;
      }

      // 'friendly' ordering of quantum numbers.  The identity must
      // compare as less than all other non-identity reps.
      static int less(int j1, int j2)
      {
         return j1 < j2;
      }

      static int size() { return 2; }

      // enumerate the possible quantum numbers from a non-negative integer.
      // This must be consistent with the ordering from the 'less' function:
      // for integers i,j, we have:
      // i < j iff less(enumerate(i), enumerate(j))
      // (this implies that the scalar (identity) quantum number is n=0
      static int enumerate(int n)
      {
         return n;
      }

      using std::tuple<FibonacciBraid> braid_rep_types;

      using default_braid_type = FibonacciBraid;

      static FibonacciBraid default_braid_rep()
      {
         return FibonacciBraid();
      }

      // Casmir invariant operators (Lie algebra terminology).
      // For finite groups, these are the expectation values of the centre
      // of the enveloping algebra (a linearly independent set that spans)
      static int num_casimir() { return 1; }

      static std::string casimir_name(std::string const& QName, int)
      { return QName; }

      static real casimir(int j, int n)
      {
         DEBUG_CHECK_EQUAL(n, 0);
         return j;
      }

      // true if this is the scalar (identity) quantum number
      static bool is_scalar(int j)
      {
         return !j;
      }

      // quantum dimensions
      static real qdim(int j)
      {
         return j ? math_const::sqrt_2 : real(1.0);
      }

      static int adjoint(int j)
      {
         return j;
      }

      // cross product is just ordinary product
      static int cross_product_exists(int j1, int j2)
      {
         return true;
      }

      static int cross_product_transforms_as(int j1, int j2)
      {
         return j1 ^ j2;
      }

      static real cross_product_factor(int j1, int j2)
      {
         return 1;
      }

      static real coupling_3j_phase(int j1, int j2, int j)
      {
         return 1;
      }

      static real fusion_31(int a, int b, int c, int d, int ab, int bc)
      {
         static constexpr real phi1 = math_const::phi;
         static constexpr real phi12 = math_const::sqrt_phi;

         if (a == 1 && b == 1 && c == 1 && d == 1)
         {
            if (ab == 0)
            {
               return (bc == 0) ? phi1 : phi12;
            }
            else
            {
               return (bc == 0) ? phi12 : -phi1;
            }
         }
         return 1;
      }

      static real fusion_22(int a, int b, int c, int d, int e, int f);

      static real coupling_6j(int j1, int j2, int j3,
                              int j4, int j5, int j6)
      {
      }

      static real coupling_9j(int j1 , int j2 , int j12,
                              int j3 , int j4 , int j34,
                              int j13, int j24, int j  )
      {
      }

      static bool is_transform_target(int j1, int j2, int j)
      {
         return (j == j1 + j2) || (j1 == 1 && j2 == 1);
      }

      static int num_transform_targets(int j1, int j2)
      {
         return (j1 == 1 && j2 == 1) ? 2 : 1;
      }

      template <typename OutIter>
      static OutIter transform_targets(int j1, int j2, OutIter Out)
      {
         if (j1 == 1 && j2 == 1)
         {
            *Out++ = 0;
            *Out++ = 1;
         }
         else
         {
            *Out++ = j1+j2;
         }
         return Out;
      }

      static difference_type difference(int j1, int j2)
      {
         return j1 - j2;
      }

      static bool is_possible_change(int j, difference_type delta)
      {
         return j+delta >= 0 && j+delta <= 1;
      }

      static int change_by(int j, difference_type delta)
      {
         return j + delta;
      }

      static real weight(difference_type x)
      {
         return abs(x);
      }
};

std::string
Fibonacci::as_string(int j)
{
   return j ? "t" : "1";
}

int
Fibonacci::from_string(std::string const& s)
{
   if (s == "1")
      return 0;
   if (s == "t")
      return 1;
   throw invalid_quantum_number("Fibonacci", s);
}

namespace
{
NiftyCounter::nifty_counter<RegisterSymmetry<Fibonacci>> FibonacciCounter;
}

} // namespace QuantumNumbers

#endif
