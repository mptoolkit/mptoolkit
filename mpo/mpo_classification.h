// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/mpo_classification.h
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

#if !defined(MPTOOLKIT_MPO_MPO_CLASSIFICATION_H)
#define MPTOOLKIT_MPO_MPO_CLASSIFICATION_H

struct OperatorClassification
{
   // indicates that the operator is zero
   bool is_null() const;

   // indicates that the operator a product state, ie a product of 1x1 MPO's
   bool is_product() const;

   // indicates that the operator is a unitary product state, ie a string operator
   bool is_unitary() const;

   // indicates that the operator is proportional to a unitary product state,
   // up to some factor
   bool is_prop_unitary() const;

   // indicates that the operator is proportional to the identity operator
   bool is_prop_identity() const;

   // returns true if the operator is the identity multiplied by a complex phase factor of magnitude 1
   bool is_complex_identity() const
   {
      return this->is_prop_identity() && norm_frob(norm_frob(this->factor())-1.0) < 1E-12;
   }

   // indicates that the operator is equal to the identity
   bool is_identity() const;

   // returns true only if the operator fits into no other classification
   bool is_unclassified() const;

   // for operators that are proportional to the identity, returns the factor
   std::complex<double> factor() const;

   // private use only
   std::complex<double> Factor_;
   bool Product_;
   bool Unitary_;
   bool Identity_;
   bool PropUnitary_;
   bool PropIdentity_;
   bool Null_;

   OperatorClassification();
};

std::ostream& operator<<(std::ostream& out, OperatorClassification const& Class);

OperatorClassification classify(CompressedMPO const& Op, double UnityEpsilon);

inline
OperatorClassification classify(CompressedMPO const& Op)
{
   return classify(Op, DefaultClassifyUnityEpsilon);
}

OperatorClassification classify(OperatorComponent c, double UnityEpsilon);

inline
OperatorClassification classify(OperatorComponent c)
{
   return classify(c, DefaultClassifyUnityEpsilon);
}

#endif
