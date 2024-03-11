// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/linearwavefunction.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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
// LinearWavefunction: main class to represent a linear matrix product wavefunction.
//
// An MPWavefunction is currently a typedef for a LinearWavefunction, but MPWavefunction
// has some additional properties.  Assumed to be always in 'normal' form; with all
// matrices satisfying the right-orthonormality constraint except possibly the left most A-matrix.
//
// Eventually an MPWavefunction will be a different beast, that allows Y-junctions
// (and possibly loops).
//
// Some of the operations defined below on LinearWavefunction require properly normalized
// states; these operations will be removed once the final MPWavefunction class is constructed.
//

#if !defined(LINEARWAVEFUNCTION_H_FUIYT49786Y709)
#define LINEARWAVEFUNCTION_H_FUIYT49786Y709

#include "mpstate.h"
#include "pheap/pvalueptr.h"
#include "pheap/pvalueiterator.h"
#include "interface/attributes.h"
#include "density.h"
#include "linearoperator.h"

class LinearWavefunction
{
   public:
      typedef MPStateComponent value_type;
      typedef pvalue_handle<MPStateComponent> handle_type;
      typedef std::list<handle_type> container_type;
      typedef container_type::iterator base_iterator;
      typedef container_type::const_iterator const_base_iterator;
      typedef pvalue_handle_iterator<base_iterator> iterator;
      typedef const_pvalue_handle_iterator<const_base_iterator> const_iterator;

      LinearWavefunction() {}

      explicit LinearWavefunction(SymmetryList const& sl) : SList(sl) {}

      // construction from an array of handles
      template <typename FwdIter>
      LinearWavefunction(FwdIter first, FwdIter last)
         : SList(first->GetSymmetryList()), Data(first, last) {}

      SymmetryList GetSymmetryList() const { return this->begin()->GetSymmetryList(); }

      std::size_t size() const { return Data.size(); }
      bool empty() const { return Data.empty(); }

      // Returns true if the state transforms irreducibly.  This is true
      // iff the left basis contains only a single site.
      bool is_irreducible() const;

      // precondition: is_irreducible
      QuantumNumbers::QuantumNumber TransformsAs() const;

      // returns the left-most basis.
      VectorBasis Basis1() const;
      // returns the right-most basis.
      VectorBasis Basis2() const;

      // This is probably going away eventually - it is really an MPWavefunction operation
      void normalize();

      iterator begin() { return iterator(Data.begin()); }
      iterator end() { return iterator(Data.end()); }

      const_iterator begin() const { return const_iterator(Data.begin()); }
      const_iterator end() const { return const_iterator(Data.end()); }

      void push_front(value_type const& x)
      { Data.push_front(handle_type(new value_type(x))); }

      void push_back(value_type const& x)
      { Data.push_back(handle_type(new value_type(x))); }

      // Because the items are stored by handle, we can't provide a
      // reference-returning front() or back() function.  Instead,
      // we use get/set.
      value_type get_front() const;
      value_type get_back() const;

      void set_front(value_type const& x);
      void set_back(value_type const& x);

      iterator insert(iterator pos, value_type const& x)
      {
         return iterator(Data.insert(pos.base(), handle_type(new value_type(x))));
      }

      iterator erase(iterator pos)
      {
         return iterator(Data.erase(pos.base()));
      }

      // Interface for accessing directly the handle_type

      base_iterator base_begin() { return Data.begin(); }
      base_iterator base_end() { return Data.end(); }

      const_base_iterator base_begin() const { return Data.begin(); }
      const_base_iterator base_end() const { return Data.end(); }

      void push_front(handle_type const& x)
     { Data.push_front(x); }

      void push_back(handle_type const& x)
      { Data.push_back(x); }

      base_iterator insert(base_iterator pos, handle_type const& h)
      {
         return Data.insert(pos, h);
      }

      AttributeList const& Attributes() const { return Attr; }
      AttributeList& Attributes() { return Attr; }

      AttributeList& AttributesMutable() const { return Attr; }

   private:
      SymmetryList SList;
      container_type Data;
      mutable AttributeList Attr;

   friend PStream::opstream& operator<<(PStream::opstream& out, LinearWavefunction const& psi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, LinearWavefunction& psi);
};

// Multiplies through my M on the left
LinearWavefunction inject_left_old_interface(MatrixOperator& M, LinearWavefunction const& Psi);
// Multiplies through by M on the right
LinearWavefunction inject_right_old_interface(LinearWavefunction const& Psi, MatrixOperator& M);

LinearWavefunction operator*(double a, LinearWavefunction const& x);
LinearWavefunction operator*(LinearWavefunction const& x, double a);
LinearWavefunction operator*(std::complex<double> a, LinearWavefunction const& x);
LinearWavefunction operator*(LinearWavefunction const& x, std::complex<double> a);

LinearWavefunction& operator*=(LinearWavefunction& psi, double a);
LinearWavefunction& operator*=(LinearWavefunction& psi, std::complex<double> a);

LinearWavefunction operator+(LinearWavefunction const& x, LinearWavefunction const& y);
LinearWavefunction operator-(LinearWavefunction const& x, LinearWavefunction const& y);

// project a (reducible) wavefunction onto an irreducible component
void project(LinearWavefunction& x, QuantumNumbers::QuantumNumber const& q);

LinearWavefunction conj(LinearWavefunction const& x);

// calculates the inner product <Psi1|Psi2>
std::complex<double>
overlap(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2);

// calculates the inner product <Psi1|conj(Psi2)>
std::complex<double>
overlap_conj(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2);

MPMatrix
reduced_matrix_element(LinearWavefunction const& Psi1,
                       LinearOperator const& M,
                       LinearWavefunction const& Psi2);

MPMatrix
reduced_matrix_element_conj(LinearWavefunction const& Psi1,
                            LinearOperator const& M,
                            LinearWavefunction const& Psi2);

std::complex<double>
expectation(LinearWavefunction const& Psi1,
            LinearOperator const& M,
            LinearWavefunction const& Psi2);

std::complex<double>
expectation_conj(LinearWavefunction const& Psi1,
                 LinearOperator const& M,
                 LinearWavefunction const& Psi2);

// action of an operator on a wavefunction (exact)
LinearWavefunction prod(LinearOperator const& Op, LinearWavefunction const& Psi);

LinearWavefunction prod(LinearOperator const& Op, LinearWavefunction const& Psi,
                        QuantumNumbers::QuantumNumber const& q);

// from left to right, calculates the action of the transfer operator R = A^\dagger m A
MatrixOperator transfer_from_left(MatrixOperator const& m, LinearWavefunction const& Psi);

// from right to left, calculates the action of the transfer operato R = A m A^\dagger
MatrixOperator transfer_from_right(MatrixOperator const& m, LinearWavefunction const& Psi);

HermitianProxy<LinearWavefunction>
inline
herm(LinearWavefunction const& x)
{
   return HermitianProxy<LinearWavefunction>(x);
}

// from left to right, calculates the action of the transfer operator R = A^\dagger m B
MatrixOperator
operator_prod(HermitianProxy<LinearWavefunction> const& A,
              MatrixOperator const& m,
              LinearWavefunction const& B);

// from right to left, calculates the action of the transfer operator R = A m B^\dagger
MatrixOperator
operator_prod(LinearWavefunction const& A,
              MatrixOperator const& m,
              HermitianProxy<LinearWavefunction> const& B);

double norm_2_sq(LinearWavefunction const& Psi);

inline
double norm_2(LinearWavefunction const& Psi)
{
   return std::sqrt(norm_2_sq(Psi));
}

inline
double scalar_difference_sq(LinearWavefunction const& x, LinearWavefunction const& y)
{
   return norm_2_sq(x) + norm_2_sq(y) - 2 * overlap(x,y).real();
}

inline
double scalar_difference(LinearWavefunction const& x, LinearWavefunction const& y)
{
   return std::sqrt(scalar_difference(x,y));
}

// Truncate the wavefunction to the given parameters.  If ShowStates=true, some info is written to std::cerr
// as the state is truncated.
void truncate(LinearWavefunction& Psi, StatesInfo const& SInfo, bool ShowStates = false);

// this belongs somewhere else

MatrixOperator CollapseBasis(VectorBasis const& b);

double const EigenvalueEpsilon = std::numeric_limits<double>::epsilon() * 4;

template <typename T>
typename BasicMPStateComponent<T>::OperatorType TruncateBasis1(BasicMPStateComponent<T>& A)
{
   typedef typename BasicMPStateComponent<T>::OperatorType OperatorType;
   OperatorType Trunc = ExpandBasis1(A);  // the original component is prod(Trunc, A)
   // Do the singular value decomposition via a (reduced) density matrix
   DensityMatrix<OperatorType> DM(scalar_prod(herm(Trunc), Trunc));

   //   DM.DensityMatrixReport(std::cerr);
   typename DensityMatrix<OperatorType>::const_iterator E = DM.begin();
   // take all eigenvalues that are bigger than EigenvalueEpsilon (normalized to EigenSum)
   while (E != DM.end() && E->Eigenvalue > EigenvalueEpsilon * DM.EigenSum()) ++E;
   OperatorType U = DM.ConstructTruncator(DM.begin(), E);
   // apply the truncation
   A = prod(U, A);
   return Trunc * herm(U);
}

template <typename T>
typename BasicMPStateComponent<T>::OperatorType TruncateBasis2(BasicMPStateComponent<T>& A)
{
   typedef typename BasicMPStateComponent<T>::OperatorType OperatorType;
   OperatorType Trunc = ExpandBasis2(A);  // the original component is prod(A, Trunc)
   // Do the singular value decomposition via a (reduced) density matrix
   DensityMatrix<OperatorType> DM(scalar_prod(Trunc, herm(Trunc)));

   // DM.DensityMatrixReport(std::cerr);
   typename DensityMatrix<OperatorType>::const_iterator E = DM.begin();
   // take all eigenvalues that are bigger than EigenvalueEpsilon (normalized to EigenSum)
   while (E != DM.end() && E->Eigenvalue > EigenvalueEpsilon * DM.EigenSum()) ++E;
   OperatorType U = DM.ConstructTruncator(DM.begin(), E);
   // apply the truncation
   A = prod(A, herm(U));
   return prod(U, Trunc, QuantumNumber(Trunc.GetSymmetryList()));
}

#endif
