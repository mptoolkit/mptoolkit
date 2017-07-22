// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/semion.h
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

#if !defined(MPTOOKLIT_QUANTUMNUMBERS_SEMION_H)
#define MPTOOKLIT_QUANTUMNUMBERS_SEMION_H

#include "common/types.h"
#include "common/niftycounter.h"
#include "common/halfint.h"
#include "quantumnumber.h"
#include "coupling_su2.h"
#include <iostream>
#include <type_traits>

namespace QuantumNumbers
{

class Semion;

class SemionBraid
{
   public:
      using is_real    = std::false_type;
      using value_type = bool;

      explicit SemionBraid(Semion const&) {}

      static std::string static_name() { return "semion"; }

      static complex r_matrix(bool a, bool b, bool c)
      {
         DEBUG_CHECK(a^b == c);
         return (a && b) ? complex(0,1) : complex(1,0);
      }
};

class AntiSemionBraid
{
   public:
      using is_real    = std::false_type;
      using value_type = bool;

      explicit AntiSemionBraid(Semion const&) {}

      static std::string static_name() { return "antisemion"; }

      static complex r_matrix(bool a, bool b, bool c)
      {
         DEBUG_CHECK(a^b == c);
         return (a && b) ? complex(0,-1) : complex(1,0);
      }
};

class Semion
{
   public:
      using value_type      = bool;
      using difference_type = bool;

      using is_pointed      = std::true_type;

      Semion() {}

      // conversion to/from string representation
      static std::string as_string(bool j);
      static bool from_string(std::string const& s);

      // converstion to/from int storage (not to be confused with the
      // integer enumeration!)
      static bool from_int(int n)
      {
         return n;
      }

      static int to_int(bool j)
      {
         return j;
      }

      static std::string static_name() { return "semion"; }

      // 'friendly' ordering of quantum numbers.  The identity must
      // compare as less than all other non-identity reps.
      static bool less(bool j1, bool j2)
      {
         return j1 < j2;
      }

      static int size() { return 2; }

      // enumerate the possible quantum numbers from a non-negative integer.
      // This must be consistent with the ordering from the 'less' function:
      // for integers i,j, we have:
      // i < j iff less(enumerate(i), enumerate(j))
      // (this implies that the scalar (identity) quantum number is n=0
      static bool enumerate(int n)
      {
         return n;
      }

      using std::tuple<SemionBraid, AntiSemionBraid> braid_rep_types;

      // no default braid type

      // Casmir invariant operators (Lie algebra terminology).
      // For finite groups, these are the expectation values of the center
      // of the enveloping algebra,
      static int num_casimir() { return 1; }

      static std::string casimir_name(std::string const& QName, int)
      { return QName; }

      static real casimir(bool j, int n)
      {
         DEBUG_CHECK_EQUAL(n, 0);
         return j;
      }

      // true if this is the scalar (identity) quantum number
      static bool is_scalar(bool j)
      {
         return !j;
      }

      // degree of the representation
      // this is redundant for a pointed category
      static int degree(bool j)
      {
         return 1;
      }

      static bool adjoint(bool j)
      {
         return j;
      }

      // Frobenius-Schur indicator
      static real fb_indicator(bool j)
      {
         return j ? -1 : 1;
      }

      // cross product is just ordinary product.
      // In this case we shouldn't need to define it explicitly
      static bool cross_product_exists(bool j1, bool j2)
      {
         return true;
      }

      static bool cross_product_transforms_as(bool j1, bool j2)
      {
         return j1 ^ j2;
      }

      static real cross_product_factor(bool j1, bool j2)
      {
         return 1;
      }

      static real coupling_3j_phase(bool j1, bool j2, bool j)
      {
         return 1;
      }

      static real fusion_31(bool a, bool b, bool c, bool d, bool ab, bool bc)
      {
         DEBUG_CHECK_EQUAL(a^b, ab);
         DEBUG_CHECK_EQUAL(b^c, bc);
         DEBUG_CHECK_EQUAL(ab^c, d);
         DEBUG_CHECK_EQUAL(a^bc, d);

         return (a && b && c && d) ? -1 : 1;
      }

      static real fusion_22(bool a, bool b, bool c, bool d, bool e, bool f);

      static real coupling_6j(bool j1, bool j2, bool j3,
                              bool j4, bool j5, bool j6)
      {
      }

      static real coupling_9j(bool j1 , bool j2 , bool j12,
                              bool j3 , bool j4 , bool j34,
                              bool j13, bool j24, bool j  )
      {
      }

      static bool product(bool j1, bool j2)
      {
         return j1 ^ j2;
      }

      // following functions are redundant for a pointed category
      static bool is_transform_target(bool j1, bool j2, bool j)
      {
         return j == (j1 ^ j2);
      }

      static int num_transform_targets(bool j1, bool j2)
      {
         return 1;
      }

      template <typename OutIter>
      static OutIter transform_targets(bool j1, bool j2, OutIter Out)
      {
         *Out++ = j1 ^ j2;
         }
         return Out;
      }

      static difference_type difference(bool j1, bool j2)
      {
         return j1 ^ j2;
      }

      static bool is_possible_change(bool j, difference_type delta)
      {
         return true;
      }

      static bool change_by(bool j, difference_type delta)
      {
         return j ^ delta;
      }

      static real weight(difference_type x)
      {
         return int(x);
      }
};

namespace
{
NiftyCounter::nifty_counter<RegisterSymmetry<Semion>> SemionCounter;
}

} // namespace QuantumNumbers

#endif
