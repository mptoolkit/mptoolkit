// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/wavefunctionsection.h
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
//
// WavefunctionSection: class to represent an MPS that is in (left or right) canonical form,
// but doesn't allow for boundary lambda matrices.  That is, the Basis1() and Basis2() are fixed
// so that we cannot obtain diagonal lambda matrices at the edges.
// Schematically, this looks liks
//
//  |          |                     |
// -A- lambda -A- lambda ... lambda -A-
//
// The difference from CanonicalWavefunction is that we do not have lambda matrices 
// at the boundaries.
//
// This is a base class for common functions between variants of canonical wavefunctions,
// eg left, right, symmetric.
//
// The Lambda matrix is stored at each parition, from 1 .. size()-1.  Hence the number of partitions
// is number of sites -1.
//
// A WavefunctionSection is generally going to be read-only.  Non-const
// iterators and access is defined as protected members, meant for construction
// and transformations that preserve the canonical form.
//

#if !defined(MPTOOLKIT_WAVEFUNCTION_WAVEFUNCTIONSECTION_H)
#define MPTOOLKIT_WAVEFUNCTION_WAVEFUNCTIONSECTION_H

#include "mps/state_component.h"
#include "linearalgebra/diagonalmatrix.h"
#include "pheap/pvalueiterator.h"
#include "wavefunction/canonicalwavefunction.h"

class LinearWavefunction;

// Common functions for Left and Right-orthogonalized versions
class WavefunctionSection
{
   public:
      // storage for the MPS matrices
      typedef StateComponent                                        mps_type;
      typedef pvalue_handle<mps_type>                               mps_handle_type;
      typedef std::vector<mps_handle_type>                          mps_container_type;
      typedef mps_container_type::iterator                          base_mps_iterator;
      typedef mps_container_type::const_iterator                    const_base_mps_iterator;
      typedef pvalue_handle_iterator<base_mps_iterator>             mps_iterator;
      typedef const_pvalue_handle_iterator<const_base_mps_iterator> const_mps_iterator;

      // storage for the lambda matrix
      typedef RealDiagonalOperator                                     lambda_type;
      typedef pvalue_handle<lambda_type>                               lambda_handle_type;
      typedef std::vector<lambda_handle_type>                          lambda_container_type;
      typedef lambda_container_type::iterator                          base_lambda_iterator;
      typedef lambda_container_type::const_iterator                    const_base_lambda_iterator;
      typedef pvalue_handle_iterator<base_lambda_iterator>             lambda_iterator;
      typedef const_pvalue_handle_iterator<const_base_lambda_iterator> const_lambda_iterator;

      // number of sites in the wavefunction
      int size() const { return Data.size(); }

      // returns true if the wavefunction doesn't contain any physical sites
      bool empty() const { return Data.empty(); }

      // returns the left-most basis.
      VectorBasis Basis1() const { return Basis1_; }
      // returns the right-most basis.
      VectorBasis Basis2() const { return Basis2_; }

      SymmetryList GetSymmetryList() const { return Basis1_.GetSymmetryList(); }

      const_mps_iterator begin() const { return const_mps_iterator(Data.begin()); }
      const_mps_iterator end() const { return const_mps_iterator(Data.end()); }

      const_base_mps_iterator base_begin() const { return Data.begin(); }
      const_base_mps_iterator base_end() const { return Data.end(); }

      const_lambda_iterator lambda_begin() const { return const_lambda_iterator(Lambda.begin()); }
      const_lambda_iterator lambda_end() const { return const_lambda_iterator(Lambda.end()); }

      const_base_lambda_iterator lambda_base_begin() const { return Lambda.begin(); }
      const_base_lambda_iterator lambda_base_end() const { return Lambda.end(); }

      // return the i'th MPS matrix.  Because they are stored by handle, we can't
      // return a reference, but the tensors are reference counted anyway so a copy is cheap
      mps_type operator[](int i) const 
      { DEBUG_RANGE_CHECK_OPEN(i, 0, int(Data.size())); return *Data[i].lock(); }

      // returns the lambda matrix at partition i
      lambda_type lambda(int i) const 
      { DEBUG_RANGE_CHECK_OPEN(i, 1, int(Data.size())); return *Lambda[i-1].lock(); }

      static PStream::VersionTag VersionT;

      void check_structure() const;
      void debug_check_structure() const;

   protected:
      // don't allow construction except via derived classes
      WavefunctionSection() {}
      WavefunctionSection(WavefunctionSection const& Psi);
      WavefunctionSection(CanonicalWavefunctionBase const& Psi);

      WavefunctionSection& operator=(WavefunctionSection const& Psi);

      // non-const iterators.  Note the final underscore to prevent mismatches with the const versions

      mps_iterator begin_() { return mps_iterator(Data.begin()); }
      mps_iterator end_() { return mps_iterator(Data.end()); }

      lambda_iterator lambda_begin_() { return lambda_iterator(Lambda.begin()); }
      lambda_iterator lambda_end_() { return lambda_iterator(Lambda.end()); }

      base_mps_iterator base_begin_() { return Data.begin(); }
      base_mps_iterator base_end_(){ return Data.end(); }

      base_lambda_iterator lambda_base_begin_() { return Lambda.begin(); }
      base_lambda_iterator lambda_base_end_() { return Lambda.end(); }

      void push_back(mps_handle_type const& x) { Data.push_back(x); }
      void push_back(mps_type const& x) { Data.push_back(new mps_type(x)); }

      void push_back_lambda(lambda_type const& x) { Lambda.push_back(new lambda_type(x)); }
      void push_back_lambda(lambda_handle_type const& x) { Lambda.push_back(x); }

      void pop_back() { Data.pop_back(); }
      void pop_back_lambda() { Lambda.pop_back(); }

      void set(int i, mps_type const& A) { Data[i] = new mps_type(A); }
      void set_lambda(int i, lambda_type const& L) { Lambda[i-1] = new lambda_type(L); }

      void setBasis1(VectorBasis const& B) { Basis1_ = B; }
      void setBasis2(VectorBasis const& B) { Basis2_ = B; }

      // Read from a stream, returns the version number.  NOTE: for versions < 3,
      // lambda(0) is not initialized, and must be done by the base class.
      int ReadStream(PStream::ipstream& in);
      void WriteStream(PStream::opstream& out) const;

   private:
      mps_container_type Data;
      lambda_container_type Lambda;

      VectorBasis Basis1_, Basis2_;
};

// function to extract the local basis (as a vector of BasisList) from a wavefunction
std::vector<BasisList>
ExtractLocalBasis(WavefunctionSection const& Psi);

inline
void
WavefunctionSection::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}


#endif
