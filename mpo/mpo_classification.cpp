// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/mpo_classification.cpp
//
// Copyright (C) 2013-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mpo_classification.h"

OperatorClassification::OperatorClassification()
   : Factor_(0.0), Product_(false), Unitary_(false),
     Identity_(false), PropUnitary_(false), PropIdentity_(false), Null_(false)
{
}

bool OperatorClassification::is_null() const
{
   return Null_;
}

bool OperatorClassification::is_product() const
{
   return Product_;
}

bool OperatorClassification::is_unitary() const
{
   return Unitary_;
}

bool OperatorClassification::is_prop_unitary() const
{
   return PropUnitary_;
}

bool OperatorClassification::is_prop_identity() const
{
   return PropIdentity_;
}

bool OperatorClassification::is_identity() const
{
   return Identity_;
}

bool OperatorClassification::is_unclassified() const
{
   return !Product_ && !Null_;
}

std::complex<double> OperatorClassification::factor() const
{
   return Factor_;
}

std::ostream& operator<<(std::ostream& out, OperatorClassification const& Class)
{
   out << "null: " << Class.is_null() << '\n';
   out << "product: " << Class.is_product() << '\n';
   out << "unitary: " << Class.is_unitary() << '\n';
   out << "prop_unitary: " << Class.is_prop_unitary() << '\n';
   out << "prop_identity: " << Class.is_prop_identity() << '\n';
   out << "complex_identity: " << Class.is_complex_identity() << '\n';
   out << "identity: " << Class.is_identity() << '\n';
   out << "factor: " << Class.factor() << '\n';
   return out;
}

OperatorClassification classify(OperatorComponent c, double UnityEpsilon)
{
   OperatorClassification Result;

   // Early return if the component is null
   if (c.is_null())
   {
      Result.Null_ = true;
      return Result;
   }

   bool IsPropIdentity = true;  // true if the operator is proportional to identity
   bool IsPropUnitary = true;   // true if the operator is proportional to a unitary operator
   std::complex<double> Factor  = 1.0;

   // firstly, check to see if it is 1x1
   if (c.Basis1().size() != 1 || c.Basis2().size() != 1)
      return Result;  // default constructed return is unclassified

   if (IsPropUnitary)
   {
      SimpleRedOperator X = c(0,0);

      if (IsPropIdentity)
      {
         if (X.Basis1() != X.Basis2() || !is_pure_scalar(X))
            IsPropIdentity = false;
         else
         {
            std::complex<double> x = PropIdent(X.scalar(), UnityEpsilon);
            if (x == 0.0)
               IsPropIdentity = false;
            else
               Factor *= x;
         }
      }

      if (!IsPropIdentity)
      {
         // is it unitary?
         std::complex<double> x = PropIdent(scalar_prod(X, herm(X)), UnityEpsilon);
         std::complex<double> y = PropIdent(scalar_prod(herm(X), X), UnityEpsilon);

         if (x == 0.0 || y == 0.0)
         {
            IsPropUnitary = false;
         }
         else
         {
            Factor *= std::sqrt(x);
         }
      }
   }

   Result.Product_ = true;
   if (IsPropUnitary)
   {
      Result.PropUnitary_ = true;

      if (IsPropIdentity)
      {
         Result.PropIdentity_ = true;
         Result.Identity_ = LinearAlgebra::norm_frob_sq(Factor - std::complex<double>(1.0, 0))
            < UnityEpsilon*UnityEpsilon;

         // if we claim to be an identity operator, we might as well make it exact
         if (Result.Identity_)
            Factor = 1.0;

         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
      else
      {
         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
   }
   return Result;
}

OperatorClassification classify(GenericMPO const& Op, double UnityEpsilon)
{
   OperatorClassification Result;

   // Early return if the operator is null
   if (Op.is_null())
   {
      Result.Null_ = true;
      return Result;
   }

   bool IsPropIdentity = true;  // true if the operator is proportional to identity
   bool IsPropUnitary = true;   // true if the operator is proportional to a unitary operator
   std::complex<double> Factor  = 1.0;

   for (unsigned i = 0; i < Op.size(); ++i)
   {
      // firstly, check to see if it is 1x1
      if (Op[i].Basis1().size() != 1 || Op[i].Basis2().size() != 1)
         return Result;  // default constructed return is unclassified

      if (IsPropUnitary)
      {
         SimpleRedOperator X = Op[i](0,0);

         if (IsPropIdentity)
         {
            if (X.Basis1() != X.Basis2() || !is_pure_scalar(X))
               IsPropIdentity = false;
            else
            {
               std::complex<double> x = PropIdent(X.scalar(), UnityEpsilon);
               if (x == 0.0)
                  IsPropIdentity = false;
               else
                  Factor *= x;
            }
         }

         if (!IsPropIdentity)
         {
            // is it unitary?
            std::complex<double> x = PropIdent(scalar_prod(X, herm(X)), UnityEpsilon);
            std::complex<double> y = PropIdent(scalar_prod(herm(X), X), UnityEpsilon);

            if (x == 0.0 || y == 0.0)
            {
               IsPropUnitary = false;
            }
            else
            {
               Factor *= std::sqrt(x);
            }
         }

      }
   }

   Result.Product_ = true;
   if (IsPropUnitary)
   {
      Result.PropUnitary_ = true;

      if (IsPropIdentity)
      {
         Result.PropIdentity_ = true;
         Result.Identity_ = LinearAlgebra::norm_frob_sq(Factor - std::complex<double>(1.0, 0))
            < UnityEpsilon*UnityEpsilon;

         // if we claim to be an identity operator, we might as well make it exact
         if (Result.Identity_)
            Factor = 1.0;

         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
      else
      {
         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
   }

   return Result;
}
