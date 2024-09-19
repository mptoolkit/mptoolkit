// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/linearoperator.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
//
// LinearOperator: wrapper around a MPOpCompressed for a linear matrix product operator.
//

#if !defined(LINEAROPERATOR_H_JDCHJKEHY589758YUER89H489)
#define LINEAROPERATOR_H_JDCHJKEHY589758YUER89H489

#include "mpopcompressed.h"

class LinearOperator
{
   private:
      typedef MPOpCompressed data_type;

   public:
      typedef MPOpComponent             value_type; // same as data_type::value_type
      typedef value_type::basis1_type   basis_type;
      typedef data_type::const_iterator const_iterator;

      LinearOperator() {}

      LinearOperator(MPOpCompressed const& x) : Data(x) {}

      // returns the total number of sites this operator contains
      int size() const { return Data.size(); }

      // returns true if this is a zero operator
      bool empty() const { return Data.empty() || Data.front().Basis1().size() == 0; }
      bool is_null() const { return Data.empty() || Data.front().Basis1().size() == 0; }

      // returns true if the operator transforms irreducibly.  This is true
      // iff the left basis contains only a single state.
      bool is_irreducible() const;

      // precondition: is_irreducible
      QuantumNumbers::QuantumNumber TransformsAs() const;

      // returns the left-most basis.  This is guaranteed to contain each
      // quantum number at most once.
      basis_type const& Basis1() const { return Data.front().Basis1(); }

      basis_type const& Basis2() const { return Data.back().Basis2(); }

      // iterate over the MPOpComponents at each site
      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Data.front().GetSymmetryList(); }

      // direct access to the MPOpCompressed
      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

   private:
      data_type Data;
};

PStream::opstream& operator<<(PStream::opstream& out, LinearOperator const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, LinearOperator& op);

LinearOperator& operator*=(LinearOperator& x, double a);
LinearOperator& operator*=(LinearOperator& x, std::complex<double> a);

LinearOperator& operator+=(LinearOperator& x, LinearOperator const& y);
LinearOperator& operator-=(LinearOperator& x, LinearOperator const& y);

LinearOperator operator+(LinearOperator const& x, LinearOperator const& y);
LinearOperator operator-(LinearOperator const& x, LinearOperator const& y);

LinearOperator operator-(LinearOperator const& x);

LinearOperator operator*(double a, LinearOperator const& x);
LinearOperator operator*(LinearOperator const& x, double a);
LinearOperator operator*(std::complex<double> a, LinearOperator const& x);
LinearOperator operator*(LinearOperator const& x, std::complex<double> a);

LinearOperator prod(LinearOperator const& x, LinearOperator const& y, QuantumNumbers::QuantumNumber const& q);
LinearOperator prod(LinearOperator const& x, LinearOperator const& y);
LinearOperator operator*(LinearOperator const& x, LinearOperator const& y);

// dot product - takes into account the multiplicity to rescale the result
LinearOperator dot(LinearOperator const& x, LinearOperator const& y);

// project a (reducible) quantum number onto an irreducible component
LinearOperator project(LinearOperator const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.  Only useful for n small!
LinearOperator pow(LinearOperator const& x, int n);

// Conjugate
LinearOperator conj(LinearOperator const& x);

// Adjoint
LinearOperator adjoint(LinearOperator const& x);

// output to a stream
inline
std::ostream& operator<<(std::ostream& out, LinearOperator const& x)
{
   return out << x.data();
}

#endif
