// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/mpopcomponent.h
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
/* -*- C++ -*- $Id$

  Defines MPOpComponent, which represents an element of
  a matrix product operator.

  TODO:
    Fix remaining functions that still access the members of MPOpComponent directly,
    then make the data members private.

    It would be better style to remove data() completely.

    I prefer the naming 'local basis' over the original 'site basis'.  I intend to eventually
    add functions with 'site' -> 'local' everywhere and eventually remove the 'site' variants.

*/

#if !defined(MPOPCOMPONENT_H_SDHCKJHKJLRHGIURHYULUHR)
#define MPOPCOMPONENT_H_SDHCKJHKJLRHGIURHYULUHR

#include "mpstate.h"
#include "tensor/tensorproduct.h"

using namespace Tensor;

typedef IrredTensor<SimpleOperator, BasisList, BasisList> CompoundOperator;

template <typename Component>
struct MPOperatorComponent;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, MPOperatorComponent<T> const& Op);

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, MPOperatorComponent<T>& Op);

template <typename Component>
struct MPOperatorComponent
{
   typedef Component OperatorType;
   typedef typename OperatorType::basis1_type basis1_type;
   typedef typename OperatorType::basis2_type basis2_type;

   // TODO: this definition should be conditional on basis1_type == basis2_type
   typedef basis1_type BasisType;

   typedef IrredTensor<OperatorType, BasisList, BasisList> mapped_type;

   typedef std::map<QuantumNumber, mapped_type> DataType;

   typedef typename DataType::value_type value_type;
   typedef typename DataType::iterator iterator;
   typedef typename DataType::const_iterator const_iterator;

   MPOperatorComponent() {}

   MPOperatorComponent(BasisList const& SBasis, basis1_type const& V1, basis2_type const& V2);

   SymmetryList const& GetSymmetryList() const { return SBasis_.GetSymmetryList(); }

   BasisList const& SiteBasis() const { return SBasis_; }
   BasisList const& LocalBasis() const { return SBasis_; }

   basis1_type const& Basis1() const { return Basis1_; }
   basis2_type const& Basis2() const { return Basis2_; }

   mapped_type& operator[](QuantumNumber const& q);

   mapped_type operator[](QuantumNumber const& q) const;

   iterator begin() { return Data_.begin(); }
   iterator end() { return Data_.end(); }

   const_iterator begin() const { return Data_.begin(); }
   const_iterator end() const { return Data_.end(); }

   iterator find(QuantumNumber const& q) { return Data_.find(q); }
   const_iterator find(QuantumNumber const& q) const { return Data_.find(q); }

   bool is_null() const { return SBasis_.size() == 0; }

   // sets the matrix elements this(a,b)(i,j) = x(a,b)
   void set_operator(int i, int j, SimpleOperator const& x);

   MPOperatorComponent& operator*=(double x);

   MPOperatorComponent& operator*=(std::complex<double> x);

   static MPOperatorComponent ConstructFullBasis1(BasisList const& SiteBasis,
                                                  basis2_type const& Basis2);

   static MPOperatorComponent ConstructFullBasis2(basis1_type const& Basis1,
                                                  BasisList const& SiteBasis);

   DataType& data() { return Data_; }
   DataType const& data() const { return Data_; }

   void check_structure() const;
   void debug_check_structure() const;

   BasisList SBasis_;
   basis1_type Basis1_;
   basis2_type Basis2_;

   DataType Data_;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, MPOperatorComponent const& Op);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, MPOperatorComponent& Op);
};

// backwards-compatibility typedef
typedef MPOperatorComponent<SimpleOperator> MPOpComponent;

typedef MPOperatorComponent<SimpleOperator> MPSimpleOpComponent;
typedef MPOperatorComponent<MatrixOperator> MPMatrixOpComponent;

template <typename T>
inline
void MPOperatorComponent<T>::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

inline
MPOpComponent ConstructIdentity(BasisList const& B)
{
   QuantumNumber Ident(B.GetSymmetryList());
   BasisList Vac = make_vacuum_basis(B.GetSymmetryList());
   SimpleOperator IdentOperator = SimpleOperator::make_identity(Vac);

   MPOpComponent Result(B, Vac, Vac);
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      set_element(Result[Ident].data(), i,i, IdentOperator);
   }
   return Result;
}

bool IsProportionalIdentity(MPOpComponent const& Op);

std::complex<double>
IdentityScale(MPOpComponent const& Op);

template <typename T>
std::ostream& operator<<(std::ostream& out, MPOperatorComponent<T> const& Op);

template <typename T>
inline
MPOperatorComponent<T>::MPOperatorComponent(BasisList const& SBasis,
                                            basis1_type const& V1,
                                            basis2_type const& V2)
  : SBasis_(SBasis), Basis1_(V1), Basis2_(V2)
{
}

template <typename T>
inline
typename MPOperatorComponent<T>::mapped_type&
MPOperatorComponent<T>::operator[](QuantumNumber const& q)
{
   if (Data_.find(q) == Data_.end())
   {
      Data_[q] = mapped_type(this->SiteBasis(), q);
   }
   return Data_[q];
}

template <typename T>
inline
typename MPOperatorComponent<T>::mapped_type
MPOperatorComponent<T>::operator[](QuantumNumber const& q) const
{
   typename DataType::const_iterator I = Data_.find(q);
   if (I == Data_.end())
   {
      return mapped_type(this->SiteBasis(), q);
   }
   return I->second;
}

namespace LinearAlgebra
{

template <>
struct interface<MPOpComponent>
{
   typedef void type;
};

// hermitian conjugation

template <typename T>
struct Herm<MPOperatorComponent<T> >
{
   typedef HermitianProxy<MPOperatorComponent<T> > result_type;
   typedef MPOperatorComponent<T> const& argument_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x);
   }
};

// scalar_prod

// does Result' = sum_s A[s] * herm(B[s])
template <typename T>
struct ScalarProd<MPOperatorComponent<T>, HermitianProxy<MPOperatorComponent<T> > >
{
   typedef typename MPOperatorComponent<T>::OperatorType result_type;
   typedef MPOperatorComponent<T> const& first_argument_type;
   typedef HermitianProxy<MPOperatorComponent<T> > const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

template <typename T>
struct ScalarProd<HermitianProxy<MPOperatorComponent<T> >, MPOperatorComponent<T> >
{
   typedef typename MPOperatorComponent<T>::OperatorType result_type;
   typedef HermitianProxy<MPOperatorComponent<T> > const& first_argument_type;
   typedef MPOperatorComponent<T> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

} // namespace LinearAlgebra


// Constructs a MPOpComponent that represents the sum of A and B.
// The resulting state has Result'[s] = A[s] \oplus B[s]
MPOpComponent
tensor_sum(MPOpComponent const& A, MPOpComponent const& B,
           SumBasis<BasisList> const& B1, SumBasis<BasisList> const& B2);

// Constructs a MPOpComponent that represents the sum of A and B,
// at the left boundary of the matrix product state.
// Precondition: A.Basis1() == B.Basis1()
// The resulting state has Result'[s] = (A[s], B[s])  (row-wise concatenation)
MPOpComponent tensor_row_sum(MPOpComponent const& A,
                             MPOpComponent const& B,
                             SumBasis<BasisList> const& B2);

// Constructs a MPOpComponent that represents the sum of A and B,
// at the right boundary of the matrix product state.
// Precondition: A.Basis2() == B.Basis2()
// The resulting state has Result'[s] = ( A[s] )
//                                      ( B[s] )  (column-wise concatenation)
MPOpComponent tensor_col_sum(MPOpComponent const& A,
                             MPOpComponent const& B,
                             SumBasis<BasisList> const& B1);

MPOpComponent prod(MPOpComponent const& A, SimpleOperator const& Op);
MPOpComponent prod(SimpleOperator const& Op, MPOpComponent const& A);

MPOpComponent prod(MPOpComponent const& A, HermitianProxy<SimpleOperator> const& Op);
MPOpComponent prod(HermitianProxy<SimpleOperator> const& Op, MPOpComponent const& A);

// constructs the tensor product in the local basis, the matrices
// C[(s's),(t't)] = A[(s's)] * B[(t't)]
MPOpComponent local_tensor_prod(MPOpComponent const& A, MPOpComponent const& B);

inline
MPOpComponent operator*(MPOpComponent const& A, SimpleOperator const& Op)
{
   return prod(A, Op);
}

inline
MPOpComponent operator*(SimpleOperator const& Op, MPOpComponent const& A)
{
   return prod(Op, A);
}

inline
MPOpComponent operator*(MPOpComponent const& A, HermitianProxy<SimpleOperator> const& Op)
{
   return prod(A, Op);
}

inline
MPOpComponent operator*(HermitianProxy<SimpleOperator> const& Op, MPOpComponent const& A)
{
   return prod(Op, A);
}

// The 'E' matrix has the same structure as a MPStateComponent
typedef MPStateComponent MPMatrix;

MPMatrix operator_prod(MPOpComponent const& M,
                       MPStateComponent const& A,
                       MPMatrix const& E,
                       HermitianProxy<MPStateComponent> const& B);


MPMatrix operator_prod(HermitianProxy<MPOpComponent> const& M,
                       HermitianProxy<MPStateComponent> const& A,
                       MPMatrix const& E,
                       MPStateComponent const& B);

MPMatrix local_operator_prod(MPOpComponent const& M,
                             MPStateComponent const& E,
                             MPMatrix const& A,
                             HermitianProxy<MPStateComponent> const& F);

#if 0
// not implemented yet
MPMatrix local_operator_prod(HermitianProxy<MPOpComponent> const& M,
                             HermitianProxy<MPStateComponent> const& E,
                             MPMatrix const& A,
                             MPStateComponent const& F);
#endif

template <typename Component1, typename Component2>
struct MProd;

template <typename C1, typename C2>
typename MProd<C1, C2>::result_type
mp_prod(C1 const& x, C2 const& y,
        ProductBasis<typename C1::BasisType, typename C2::BasisType> const& b1,
        ProductBasis<typename C1::BasisType, typename C2::BasisType> const& b2)
{
   return MProd<C1,C2>()(x,y,b1,b2);
}

// product of a matrix product operator and a matrix product wavefunction.
// This is about the only place that the coupling coefficients enter at such a high
// level: the reason is the local basis label of a matrix product state
// transforms covariantly, but in a matrix*vector multiply the vector
// must transform contravariantly.  The conj_phase() fixes this.
template <>
struct MProd<MPOpComponent, MPStateComponent>
{
   typedef MPStateComponent result_type;
   typedef MPOpComponent const& first_argument_type;
   typedef MPStateComponent const& second_argument_type;

   result_type operator()(first_argument_type M, second_argument_type A,
                          ProductBasis<BasisList, VectorBasis> const& B1,
                          ProductBasis<BasisList, VectorBasis> const& B2);
};

// product of two matrix product operators
template <>
struct MProd<MPOpComponent, MPOpComponent>
{
   typedef MPOpComponent result_type;
   typedef MPOpComponent const& first_argument_type;
   typedef MPOpComponent const& second_argument_type;

   result_type operator()(first_argument_type M, second_argument_type N,
                          ProductBasis<BasisList, BasisList> const& B1,
                          ProductBasis<BasisList, BasisList> const& B2);
};

SimpleOperator ExpandBasis1(MPOpComponent& A, Normalization n = Intensive);

SimpleOperator ExpandBasis2(MPOpComponent& A, Normalization n = Intensive);

// returns the matrix that is sum_s M^{ss}
SimpleOperator
local_trace(MPOpComponent const& A);

// Basis truncation via removal of parallel components.
SimpleOperator TruncateBasis1(MPOpComponent& A);
SimpleOperator TruncateBasis2(MPOpComponent& A);

MPMatrixOpComponent
outer_prod(MPStateComponent const& a, HermitianProxy<MPStateComponent> const& b);

MatrixOperator
element_prod(MPMatrixOpComponent const& a, HermitianProxy<MPMatrixOpComponent> const& b);

template <typename Component>
MPOperatorComponent<Component>
operator+(MPOperatorComponent<Component> const& a, MPOperatorComponent<Component> const& b)
{
   typedef typename MPOperatorComponent<Component>::const_iterator const_iterator;
   MPOperatorComponent<Component> Result(a);
   for (const_iterator I = b.begin(); I != b.end(); ++I)
   {
      Result[I->first] += I->second;
   }
   return Result;
}

MPOpComponent
SumFix1(MPOpComponent const& A, MPOpComponent const& B,
        SimpleOperator& AMap, SimpleOperator& BMap);

// conjugate of the components
MPOpComponent conj(MPOpComponent const& x);

// Given M^{s's}_{a'a}, calculates M^{s,s'}_{adjoint(a'),adjoint(a)}
// That is, adjoint of the local operator indices, but take the
// adjoint matrix basis without transposing matrix indices.
MPOpComponent local_adjoint(MPOpComponent const& x);

// returns the list of local operators corresponding to a given row
// of the matrix basis.  This function only works if the MPOpComponent is irreducible,
// in the sense that the local operator at each entry (i,j) in the matrix basis
// is irreducible.
LinearAlgebra::Vector<SimpleOperator>
project_row(MPOpComponent const& x, int r);

// returns the list of local operators corresponding to a given column
// of the matrix basis.   This function only works if the MPOpComponent is irreducible,
// in the sense that the local operator at each entry (i,j) in the matrix basis
// is irreducible.
LinearAlgebra::Vector<SimpleOperator>
project_column(MPOpComponent const& x, int c);

#if 0
typedef CompoundOperator MPOpMatrix;

// returns the n'th component of
// x in the form of an MPOpMatrix.  This has the matrix basis as the outer index.
MPOpMatrix as_matrix(MPOpComponent::ComponentType const& x);
#endif

#include "mpopcomponent.cc"

#endif
