// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/linearwavefunction.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_WAVEFUNCTION_LINEARWAVEFUNCTION_H)
#define MPTOOLKIT_WAVEFUNCTION_LINEARWAVEFUNCTION_H

#include "mps/state_component.h"
#include "pheap/pvalueptr.h"
#include "pheap/pvalueiterator.h"
#include "interface/attributes.h"
#include "mps/density.h"
//#include "linearoperator.h"

class LinearWavefunction
{
   public:
      using value_type       = StateComponent;
      using handle_type      = pvalue_handle<StateComponent>;
      using pointer_type     = pvalue_ptr<StateComponent>;
      using proxy_type       = pvalue_proxy<StateComponent>;
      using const_proxy_type = const_pvalue_proxy<StateComponent>;


      typedef std::list<handle_type> container_type;
      typedef container_type::iterator base_iterator;
      typedef container_type::const_iterator const_base_iterator;
      typedef pvalue_handle_iterator<base_iterator> iterator;
      typedef const_pvalue_handle_iterator<const_base_iterator> const_iterator;

      LinearWavefunction() noexcept {}

      explicit LinearWavefunction(SymmetryList const& sl) : SList(sl) {}

      LinearWavefunction(LinearWavefunction&& Other) noexcept = default;
      LinearWavefunction(LinearWavefunction const& Other) = default;

      LinearWavefunction& operator=(LinearWavefunction&& Other) noexcept = default;
      LinearWavefunction& operator=(LinearWavefunction const& Other) = default;

      // construction from an array of handles
      template <typename FwdIter>
      LinearWavefunction(FwdIter first, FwdIter last)
         : Data(first, last) { if (first != last) SList = first->lock()->GetSymmetryList(); }

      template <typename FwdIter>
      LinearWavefunction(SymmetryList const& sl, FwdIter first, FwdIter last)
         : SList(sl), Data(first, last) {}

      // construction from a container of StateComponent's
      template <typename FwdIter>
      static
      LinearWavefunction FromContainer(FwdIter first, FwdIter last);

      SymmetryList GetSymmetryList() const { return SList; }

      std::size_t size() const { return Data.size(); }
      bool empty() const { return Data.empty(); }

      // Returns true if the state transforms irreducibly.  This is true
      // iff the left basis contains only a single site.
      bool is_irreducible() const;

      // returns true if the state can be considered to be a pure state; which
      // implies that both the right hand basis is one-dimensional.
      // We can interpret a reducible left-hand basis also as a pure state, by summing
      // the states, so we do not require that the left-hand basis is one-dimensional.
      bool is_pure() const;

      // precondition: is_irreducible
      QuantumNumbers::QuantumNumber TransformsAs() const;

      // returns the left-most basis.
      VectorBasis Basis1() const;
      // returns the right-most basis.
      VectorBasis Basis2() const;

      iterator begin() { return iterator(Data.begin()); }
      iterator end() { return iterator(Data.end()); }

      const_iterator begin() const { return const_iterator(Data.begin()); }
      const_iterator end() const { return const_iterator(Data.end()); }

      void push_front(value_type&& x)
      {
	 if (SList.is_null()) SList = x.GetSymmetryList();
	 Data.push_front(handle_type(new value_type(std::move(x))));
      }

      void push_front(value_type const& x)
      {
	 if (SList.is_null()) SList = x.GetSymmetryList();
	 Data.push_front(handle_type(new value_type(copy(x))));
      }

      void push_back(value_type&& x)
      {
	 if (SList.is_null()) SList = x.GetSymmetryList();
	 Data.push_back(handle_type(new value_type(std::move(x))));
      }

      void push_back(value_type const& x)
      {
	 if (SList.is_null()) SList = x.GetSymmetryList();
	 Data.push_back(handle_type(new value_type(copy(x))));
      }

      void push_front(LinearWavefunction const& x)
      {
         const_iterator Ibegin = x.begin();
         const_iterator I = x.end();
         while (I != x.begin())
         {
            --I;
            this->push_front(copy(*I));
         }
         if (SList.is_null())
            SList = x.GetSymmetryList();
      }

      void push_back(LinearWavefunction const& x)
      {
         const_iterator Iend = x.end();
         for (const_iterator I = x.begin(); I != Iend; ++I)
         {
            this->push_front(*I);
         }
         if (SList.is_null())
            SList = x.GetSymmetryList();
      }

      // Because the items are stored by handle, we can't provide a
      // reference-returning front() or back() function.  Instead,
      // we use get/set.
      const_proxy_type get_front() const;
      const_proxy_type get_back() const;

      void set_front(value_type const& x);
      void set_back(value_type const& x);

      proxy_type front()
      {
	 return proxy_type(Data.front().lock());
      }

      const_proxy_type front() const
      {
	 return const_proxy_type(Data.front().lock());
      }

      proxy_type back()
      {
	 return proxy_type(Data.back().lock());
      }

      const_proxy_type back() const
      {
	 return const_proxy_type(Data.back().lock());
      }

      iterator insert(iterator pos, value_type x)
      {
         return iterator(Data.insert(pos.base(), handle_type(new value_type(std::move(x)))));
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
      {
         Data.push_front(x);
         if (SList.is_null())
            SList = x.lock()->GetSymmetryList();
      }

      void push_back(handle_type const& x)
      {
         Data.push_back(x);
         if (SList.is_null())
            SList = x.lock()->GetSymmetryList();
      }

      base_iterator insert(base_iterator pos, handle_type const& h)
      {
         return Data.insert(pos, h);
      }

      void pop_front()
      {
         Data.pop_front();
      }

      void pop_back()
      {
         Data.pop_back();
      }

   private:
      SymmetryList SList;
      container_type Data;

   friend PStream::opstream& operator<<(PStream::opstream& out, LinearWavefunction const& psi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, LinearWavefunction& psi);
};

// Given a wavefunction Psi that is in right-orthogonal form,
// Multiplies M on the left hand side of Psi, and iteratively left-orthogonalizes the wavefunction,
// returning the remainder matrix.
MatrixOperator left_orthogonalize(LinearWavefunction& Psi, int Verbose = 0);
MatrixOperator left_orthogonalize(MatrixOperator M, LinearWavefunction& Psi, int Verbose = 0);

// Multiplies M on the right hand side of Psi, and iteratively right-orthogonalizes the wavefunction,
// returning the remainder matrix
MatrixOperator right_orthogonalize(LinearWavefunction& Psi, int Verbose = 0);
MatrixOperator right_orthogonalize(LinearWavefunction& Psi, MatrixOperator M, int Verbose = 0);

LinearWavefunction operator*(double a, LinearWavefunction const& x);
LinearWavefunction operator*(LinearWavefunction const& x, double a);
LinearWavefunction operator*(std::complex<double> a, LinearWavefunction const& x);
LinearWavefunction operator*(LinearWavefunction const& x, std::complex<double> a);

LinearWavefunction& operator*=(LinearWavefunction& psi, double a);
LinearWavefunction& operator*=(LinearWavefunction& psi, std::complex<double> a);

// Addition of wavefunctions.  Precondition: x.is_pure() && y.is_pure()
LinearWavefunction operator+(LinearWavefunction const& x, LinearWavefunction const& y);
LinearWavefunction operator-(LinearWavefunction const& x, LinearWavefunction const& y);

// project a (reducible) wavefunction onto an irreducible component
void project(LinearWavefunction& x, QuantumNumbers::QuantumNumber const& q);

// N-to-1 coarsegraining of a wavefunction
LinearWavefunction coarse_grain(LinearWavefunction const& x, int N);

LinearWavefunction conj(LinearWavefunction const& x);

// calculates the inner product <Psi1|Psi2>
std::complex<double>
overlap(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2);

// calculates the inner product <Psi1|conj(Psi2)>
std::complex<double>
overlap_conj(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2);

#if 0
StateComponent
reduced_matrix_element(LinearWavefunction const& Psi1,
                       LinearOperator const& M,
                       LinearWavefunction const& Psi2);

StateComponent
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
#endif

// Calculates the operator contraction, with a matrix
// actong on the left hand side of the wavefunction.
// (the transfer matrix of the unit cell)
// Two argument version, for a different wavefunction on the left and right.
// +-Psi1*- ... Psi1*-
// |  |          |
// m  |          |
// |  |          |
// +-Psi2-- ... Psi2--
MatrixOperator
inject_left(MatrixOperator m,
            LinearWavefunction const& Psi1,
            LinearWavefunction const& Psi2);

// equivalent to inject_left(m, Psi, Psi);
MatrixOperator
inject_left(MatrixOperator m, LinearWavefunction const& Psi);

// Calculates the operator contraction, with a matrix
// actong on the left hand side of the wavefunction.
// (the transfer matrix of the unit cell)
// --Psi1-- ... Psi1--+
//    |         |     |
//    |         |     m
//    |         |     |
// --Psi2*- ... Psi2*-+
// Two argument version, for a different wavefunction on the left and right.
MatrixOperator
inject_right(MatrixOperator m,
            LinearWavefunction const& Psi1,
            LinearWavefunction const& Psi2);

// equivalent to inject_right(m, Psi, Psi)
MatrixOperator
inject_right(MatrixOperator, LinearWavefunction const& Psi);


HermitianProxy<LinearWavefunction>
inline
herm(LinearWavefunction const& x)
{
   return HermitianProxy<LinearWavefunction>(x);
}

// from left to right, calculates the action of the transfer operator R = A^\dagger m B
MatrixOperator
operator_prod(HermitianProxy<LinearWavefunction> const& A,
              MatrixOperator E,
              LinearWavefunction const& B);

// from right to left, calculates the action of the transfer operator R = A m B^\dagger
MatrixOperator
operator_prod(LinearWavefunction const& A,
              MatrixOperator F,
              HermitianProxy<LinearWavefunction> const& B);

#if 0
// from right to left, calculates the action of the transfer operator R = A m B^\dagger
// with respect to the MPO Op
StateComponent
operator_prod(LinearOperator const& Op,
              LinearWavefunction const& A,
              StateComponent const& F,
              HermitianProxy<LinearWavefunction> const& B);

StateComponent
operator_prod(LinearOperator const& Op,
              HermitianProxy<LinearWavefunction> const& A,
              StateComponentOperator const& E,
              LinearWavefunction const& B);

// These versions require that Op as one-dimensional left and right bases.
// It is a shortcut for a StateComponent with a one-dimensional local basis.
MatrixOperator
operator_prod(LinearOperator const& Op,
              LinearWavefunction const& A,
              StateComponent const& F,
              HermitianProxy<LinearWavefunction> const& B);

MatrixOperator
operator_prod(LinearOperator const& Op,
              HermitianProxy<LinearWavefunction> const& A,
              MatrixOperator const& E,
              LinearWavefunction const& B);
#endif

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
   return std::sqrt(scalar_difference_sq(x,y));
}

// Truncate the wavefunction to the given parameters.  If ShowStates=true, some info is written to std::cerr
// as the state is truncated.
void truncate(LinearWavefunction& Psi, StatesInfo const& SInfo, bool ShowStates = false);

// this belongs somewhere else

MatrixOperator CollapseBasis(VectorBasis const& b);

// function to extract the local basis (as a vector of BasisList) from a wavefunction
std::vector<BasisList>
ExtractLocalBasis(LinearWavefunction const& Psi);

// function to return a StringOperator given an array of local basis
std::vector<SimpleOperator>
make_identity_string_operator(std::vector<BasisList> const& Basis);

double const EigenvalueEpsilon = std::numeric_limits<double>::epsilon() * 4;

inline
StateComponent::operator_type
TruncateBasis1(StateComponent& A)
{
   TRACE("TruncateBasis1");
   StateComponent ACheck = copy(A);
   typedef StateComponent::operator_type operator_type;
   operator_type Trunc = ExpandBasis1(A);  // the original component is prod(Trunc, A)
   //TRACE(A)(Trunc);
   // Do the singular value decomposition via a (reduced) density matrix
   DensityMatrix<operator_type> DM(scalar_prod(herm(Trunc), Trunc));

   DensityMatrix<operator_type>::const_iterator E = DM.begin();
   TRACE(DM.EigenSum());
   // take all eigenvalues that are bigger than EigenvalueEpsilon (normalized to EigenSum)
   while (E != DM.end() && E->Eigenvalue > EigenvalueEpsilon * DM.EigenSum())
   {
      ++E;
   }
   operator_type U = DM.ConstructTruncator(DM.begin(), E);
   TRACE(U);
   inplace_conj(U);
   TRACE(U);
   // apply the truncation
   A = prod(U, A);
   TRACE(Trunc.Basis1())(Trunc.Basis2())(U.Basis1())(U.Basis2());
   TRACE(norm_frob(prod(Trunc * herm(U), A)))(norm_frob(prod(Trunc * herm(U),A) - ACheck));
   CHECK(norm_frob(prod(Trunc * herm(U), A) - ACheck) < 0.05);
   CHECK(norm_frob_sq(Trunc*herm(U)) > 0.8)(Trunc)(U);
   //TRACE(Trunc)(U);
   return Trunc * herm(U);
}

inline
StateComponent::operator_type
TruncateBasis2(StateComponent& A)
{
   typedef StateComponent::operator_type OperatorType;
   OperatorType Trunc = ExpandBasis2(A);  // the original component is prod(A, Trunc)
   // Do the singular value decomposition via a (reduced) density matrix
   DensityMatrix<OperatorType> DM(scalar_prod(Trunc, herm(Trunc)));

   // DM.DensityMatrixReport(std::cerr);
   DensityMatrix<OperatorType>::const_iterator E = DM.begin();
   // take all eigenvalues that are bigger than EigenvalueEpsilon (normalized to EigenSum)
   while (E != DM.end() && E->Eigenvalue > EigenvalueEpsilon * DM.EigenSum()) ++E;
   OperatorType U = DM.ConstructTruncator(DM.begin(), E);
   // apply the truncation
   A = prod(A, herm(U));
   return prod(U, Trunc, QuantumNumber(Trunc.GetSymmetryList()));
}

// template definitions

template <typename FwdIter>
LinearWavefunction
LinearWavefunction::FromContainer(FwdIter first, FwdIter last)
{
   LinearWavefunction Result;
   while (first != last)
   {
      Result.push_back(*first++);
   }
   return Result;
}

#endif
