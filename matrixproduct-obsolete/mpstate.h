// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/mpstate.h
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
// I prefer the naming 'local basis' over the original 'site basis'.  I intend to eventually
// add functions with 'site' -> 'local' everywhere and eventually remove the 'site' variants.

#if !defined(MPSTATE_H_SDHCKJHKJLRHGIURHYULUHR)
#define MPSTATE_H_SDHCKJHKJLRHGIURHYULUHR

#include "tensor/tensor.h"
#include "tensor/tensorsum.h"
#include "tensor/basis.h"
#include "tensor/tensorproduct.h"

using namespace Tensor;

typedef IrredTensor<std::complex<double> > SimpleOperator;

typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
                            VectorBasis, 
                            VectorBasis> MatrixOperator;

template <typename T>
struct BasicMPStateComponent;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, BasicMPStateComponent<T> const& Op);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, BasicMPStateComponent<T>& Op);

template <typename T>
struct BasicMPStateComponent
{
   typedef T value_type;
   typedef T OperatorType;

   typedef LinearAlgebra::Vector<value_type> DataType;

   typedef typename DataType::const_iterator const_iterator;
   typedef typename DataType::iterator       iterator;

   typedef VectorBasis BasisType;

   BasicMPStateComponent() {}

   BasicMPStateComponent(BasisList const& SBasis_, VectorBasis const& V1, 
		    VectorBasis const& V2);

   QuantumNumbers::SymmetryList GetSymmetryList() const { return SBasis.GetSymmetryList(); }

   // returns true if this is a zero matrix
   bool is_null() const { return SBasis.size() == 0; }

   BasisList SiteBasis() const { return SBasis; }
   BasisList LocalBasis() const { return SBasis; }

   BasisType Basis1() const { return VBasis1; }
   BasisType Basis2() const { return VBasis2; }

   std::size_t size() const { return Data.size(); }

   // Transforms *this so that it satisfies the left-handed orthogonality constraint;
   // mp_prod_left(*this', *this') is the identity operator
   // and the state is recovered by (*this) = (*this') * result'
   //   OperatorType OrthogonalizeLeft();

   // Transforms *this so that it satisfies the right-handed orthogonality constraint;
   // mp_prod_right(*this', *this') is the identity operator
   // and the state is recovered by (*this) = result' * (*this')
   //   OperatorType OrthogonalizeRight();

   const_iterator begin() const { return Data.begin(); }
   const_iterator end() const { return Data.end(); }

   iterator begin() { return Data.begin(); }
   iterator end() { return Data.end(); }

   value_type const& operator[](int s) const { return Data[s]; }
   value_type& operator[](int s) { return Data[s]; }

   // front() and back() functions are useful for Hamiltonian matrix elements
   value_type& front() { return Data[0]; }
   value_type const& front() const { return Data[0]; }

   value_type& back() { return Data[Data.size()-1]; }
   value_type const& back() const { return Data[Data.size()-1]; }

   void delta_shift(QuantumNumber const& q, QuantumNumbers::Projection const& Delta);

   static BasicMPStateComponent<T> ConstructFullBasis1(BasisList const& SiteBasis, 
                                                       VectorBasis const& Basis2);
   static BasicMPStateComponent<T> ConstructFullBasis2(VectorBasis const& Basis1, 
                                                       BasisList const& SiteBasis);
   
   BasisList SBasis;
   VectorBasis VBasis1, VBasis2;
   DataType Data;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, BasicMPStateComponent const& Op);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, BasicMPStateComponent& Op);

};

template <typename T>
void BasicMPStateComponent<T>::delta_shift(QuantumNumbers::QuantumNumber const& q,
					   QuantumNumbers::Projection const& Delta)
{
   DEBUG_TRACE("before")(scalar_prod(*this, herm(*this)))(scalar_prod(herm(*this), *this));
   VBasis1 = DeltaShift(VBasis1, Delta);
   VBasis2 = DeltaShift(VBasis2, Delta);
   for (typename BasicMPStateComponent<T>::iterator I = this->begin(); I != this->end(); ++I)
   {
      DEBUG_TRACE(*I)(I->Basis1())(I->Basis2());
      *I = Tensor::delta_shift(*I, q, Delta, VBasis1, VBasis2);
      DEBUG_TRACE("after component")(scalar_prod(*I, herm(*I)))(scalar_prod(herm(*I), *I))(I->TransformsAs());
      DEBUG_TRACE(*I)(I->Basis1())(I->Basis2());
   }
   DEBUG_TRACE("after")(scalar_prod(*this, herm(*this)))(scalar_prod(herm(*this), *this));
}

typedef BasicMPStateComponent<MatrixOperator> MPStateComponent;

inline
MPStateComponent& operator*=(MPStateComponent& x, double y)
{
   for (unsigned i = 0; i < x.size(); ++i)
      x[i] *= y;
   return x;
}

inline
MPStateComponent& operator*=(MPStateComponent& x, std::complex<double> y)
{
   for (unsigned i = 0; i < x.size(); ++i)
      x[i] *= y;
   return x;
}

std::ostream& operator<<(std::ostream& out, MPStateComponent const& Psi);

MPStateComponent make_vacuum_state(QuantumNumbers::SymmetryList const& S);

MPStateComponent make_vacuum_state(QuantumNumbers::QuantumNumber const& Q);

inline
MatrixOperator make_vacuum_matrix(QuantumNumbers::SymmetryList const& SList)
{
   VectorBasis B(SList);
   B.push_back(QuantumNumbers::QuantumNumber(SList), 1);
   return MatrixOperator::make_identity(B);
}

inline
SimpleOperator make_vacuum_simple(QuantumNumbers::SymmetryList const& SList)
{
   BasisList B(SList);
   B.push_back(QuantumNumbers::QuantumNumber(SList));
   return SimpleOperator::make_identity(B);
}

MPStateComponent
operator+(MPStateComponent const& x, MPStateComponent const& y);

MPStateComponent
operator-(MPStateComponent const& x, MPStateComponent const& y);

MPStateComponent&
operator+=(MPStateComponent& x, MPStateComponent const& y);

MPStateComponent&
operator-=(MPStateComponent& x, MPStateComponent const& y);

namespace LinearAlgebra
{

template <>
struct interface<MPStateComponent>
{
   typedef void type;
};

// hermitian conjugation

template <>
struct Herm<MPStateComponent>
{
   typedef HermitianProxy<MPStateComponent> result_type;
   typedef MPStateComponent const& argument_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x);
   }
};

// complex conjugation

inline
MPStateComponent conj(MPStateComponent const& x)
{
   MPStateComponent Result(x);
   for (MPStateComponent::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = conj(*I);
   }
   return Result;
}

// scalar_prod

// does Result' = sum_s A[s] * herm(B[s])
// generalized ScalarProd.  This means that,
// for example, scalar_direct_prod works automatically.

template <typename Func>
struct ScalarProd<MPStateComponent, HermitianProxy<MPStateComponent>, Func>
{
   typedef typename Func::result_type result_type;
   typedef MPStateComponent const& first_argument_type;
   typedef HermitianProxy<MPStateComponent> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
   result_type operator()(first_argument_type x, second_argument_type y, Func f) const;
};

template <typename Func>
struct ScalarProd<HermitianProxy<MPStateComponent>, MPStateComponent, Func>
{
   typedef typename Func::result_type result_type;
   typedef HermitianProxy<MPStateComponent> const& first_argument_type;
   typedef MPStateComponent const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
   result_type operator()(first_argument_type x, second_argument_type y, Func f) const;
};

// the common case is explicitly specialized
template <>
struct ScalarProd<MPStateComponent, HermitianProxy<MPStateComponent> >
{
   typedef MatrixOperator result_type;
   typedef MPStateComponent const& first_argument_type;
   typedef HermitianProxy<MPStateComponent> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

template <>
struct ScalarProd<HermitianProxy<MPStateComponent>, MPStateComponent>
{
   typedef MatrixOperator result_type;
   typedef HermitianProxy<MPStateComponent> const& first_argument_type;
   typedef MPStateComponent const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

} // namespace LinearAlgebra


// some quick hacks
inline
std::complex<double>
inner_prod(MPStateComponent const& x, MPStateComponent const& y)
{
   return trace(scalar_prod(x, herm(y)));
}

inline
double
norm_frob_sq(MPStateComponent const& x)
{
   return trace(scalar_prod(x, herm(x))).real();
}

inline
double
norm_frob(MPStateComponent const& x)
{
   return std::sqrt(norm_frob_sq(x));
}

inline
MPStateComponent operator*(MPStateComponent const& x, double y)
{
   MPStateComponent Res(x);
   Res *= y;
   return Res;
}

inline
MPStateComponent operator*(MPStateComponent const& x, std::complex<double> y)
{
   MPStateComponent Res(x);
   Res *= y;
   return Res;
}

inline
MPStateComponent operator*(double y, MPStateComponent const& x)
{
   MPStateComponent Res(x);
   Res *= y;
   return Res;
}

inline
MPStateComponent operator*(std::complex<double> y, MPStateComponent const& x)
{
   MPStateComponent Res(x);
   Res *= y;
   return Res;
}

// does Result' = sum_{s,t} M(t,s) * A^t * herm(B^s)
MatrixOperator operator_prod(SimpleOperator const& M, 
                             MPStateComponent const& A, 
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B);

// does Result' = sum_{s,t} M(t,s) * herm(A^t) * B^s
MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M, 
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& A, 
                             MPStateComponent const& B);

// does Result' = sum_{s,t} M(t,s) * A^t * E * herm(B^s)
MatrixOperator operator_prod(SimpleOperator const& M, 
                             MPStateComponent const& A, 
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B,
                             QuantumNumbers::QuantumNumber const& q);

#if 0
// including a delta shift.
MatrixOperator operator_prod_delta(SimpleOperator const& M, 
                                   MPStateComponent const& A, 
                                   MatrixOperator const& E,
                                   LinearAlgebra::HermitianProxy<MPStateComponent> const& B,
                                   QuantumNumbers::QuantumNumber const& q,
                                   QuantumNumbers::Projection const& p);
#endif

MatrixOperator operator_prod(SimpleOperator const& M, 
                             MPStateComponent const& A, 
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B);

// does Result' = sum_{s,t} M(t,s) * herm(A^t) * E * B^s
MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M, 
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& A, 
                             MatrixOperator const& E,
                             MPStateComponent const& B,
                             QuantumNumbers::QuantumNumber const& q);

#if 0
// does Result' = sum_{s,t} M(t,s) * herm(A^t) * E * B^s,
// including a delta shift.
MatrixOperator operator_prod_delta(LinearAlgebra::HermitianProxy<SimpleOperator> const& M, 
                                   LinearAlgebra::HermitianProxy<MPStateComponent> const& A, 
                                   MatrixOperator const& E,
                                   MPStateComponent const& B,
                                   QuantumNumbers::QuantumNumber const& q,
                                   QuantumNumbers::Projection const& p);
#endif

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M, 
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& A, 
                             MatrixOperator const& E,
                             MPStateComponent const& B);

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<MPStateComponent> const& A, 
                             MatrixOperator const& E,
                             MPStateComponent const& B);

MatrixOperator operator_prod(MPStateComponent const& A, 
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B);

// returns the operator M(s,t) = trace(herm(A^s) B^t)
// This is necessarily a scalar operator (the trace of a non-scalar operator is zero).
SimpleOperator trace_prod(LinearAlgebra::HermitianProxy<MPStateComponent> const& A,
                          MPStateComponent const& B);

// returns the operator M(s,t) = trace(A^s herm(B^t))
// This is necessarily a scalar operator (the trace of a non-scalar operator is zero).
SimpleOperator trace_prod(MPStateComponent const& A,
                          LinearAlgebra::HermitianProxy<MPStateComponent> const& B);

// product of a MPStateComponent and a MatrixOperator.
// The BasicOperator must transform as a scalar - eg typically it will be unitary.
MPStateComponent prod(MPStateComponent const& A, MatrixOperator const& Op);

MPStateComponent prod(MPStateComponent const& A, HermitianProxy<MatrixOperator> const& Op);

// product of a BasicOperator and a MPStateComponent.
// The BasicOperator must transform as a scalar - eg typically it will be unitary.
MPStateComponent prod(MatrixOperator const& Op, MPStateComponent const& A);

MPStateComponent prod(HermitianProxy<MatrixOperator> const& Op, MPStateComponent const& A);

// returns the operator a'^s' = sum_s A^s * x(s,s').  x must transform as a scalar.
MPStateComponent local_prod(MPStateComponent const& A, SimpleOperator const& x);

// returns the operator a'^s' = sum_s x(s',s) * A^s.  x must transform as a scalar.
MPStateComponent local_prod(SimpleOperator const& x, MPStateComponent const& A);

// constructs the tensor product in the local basis, the matrices
// C[k1,k2] = A[k1] * B[k2]
MPStateComponent local_tensor_prod(MPStateComponent const& A, MPStateComponent const& B);

MPStateComponent triple_prod(MatrixOperator const& Op1, 
                             MPStateComponent const& A, 
                             LinearAlgebra::HermitianProxy<MatrixOperator> const&Op2);

MPStateComponent triple_prod(LinearAlgebra::HermitianProxy<MatrixOperator> const& Op1, 
                             MPStateComponent const& A, 
                             MatrixOperator const&Op2);

// Constructs a MPStateComponent that represents the sum of A and B.
// The resulting state has Result'[s] = A[s] \oplus B[s]
MPStateComponent tensor_sum(MPStateComponent const& A, MPStateComponent const& B, 
                            SumBasis<VectorBasis> const& B1, 
                            SumBasis<VectorBasis> const& B2);

// Constructs a MPStateComponent that represents the sum of A and B,
// at the left boundary of the matrix product state.
// Precondition: A.Basis1() == B.Basis1()
// The resulting state has Result'[s] = (A[s], B[s])  (row-wise concatenation)
MPStateComponent tensor_row_sum(MPStateComponent const& A, 
                                MPStateComponent const& B, 
                                SumBasis<VectorBasis> const& B2);

// Constructs a MPStateComponent that represents the sum of A and B,
// at the right boundary of the matrix product state.
// Precondition: A.Basis2() == B.Basis2()
// The resulting state has Result'[s] = ( A[s] )
//                                      ( B[s] )  (column-wise concatenation)
MPStateComponent tensor_col_sum(MPStateComponent const& A, 
                                MPStateComponent const& B, 
                                SumBasis<VectorBasis> const& B1);

// Returns the diagonal components of the operator F given by
// F(x) = operator_prod(A, x, herm(B))
MatrixOperator extract_diagonal(MPStateComponent const& A, 
                                LinearAlgebra::HermitianProxy<MPStateComponent> const& B);

// shift the basis of an MPStateComponent by some quantum number; the local basis
// is shifted by QL, the incoming matrix basis is shifted by QM, the outgoing matrix basis
// is shifted by QM+QL.  This only works if QL and QM are degree 1 reps.
MPStateComponent ShiftLocalBasis(MPStateComponent const& Op, QuantumNumber QL, QuantumNumber QM);

MPStateComponent delta_shift(MPStateComponent const& Op, QuantumNumber const& q);

// scales the quantum numbers of all bases of an MPStateComponent by some factor.  
// This only makes sense for U(1) quantum numbers.  Name is the name of the quantum number
// to scale.
MPStateComponent ScaleBasisU1(MPStateComponent const& Op, 
			      std::string const& Name, 
			      double Factor);

// returns a MPStateComponent that is identical except for the modified SymmetryList,
// which must have the same symmetry types in the same order, but is allowed to be
// renamed.
MPStateComponent RenameSymmetry(MPStateComponent const& Op, SymmetryList const& NewSL);

// re-orders the basis
MPStateComponent ReorderLocalBasis(MPStateComponent const& Op, std::list<int> const& NewOrder);

// returns a MPStateComponent with a symmetry list that is compatible with the old one
// in the sense of all old symmetries are in the new list (possibly in a different order).
// new symmetries are allowed, and get quantum number 0.
MPStateComponent CoerceSymmetryList(MPStateComponent const& Op, SymmetryList const& NewSL);

// expand the basis

enum Normalization { Intensive, Extensive };

MatrixOperator ExpandBasis1(MPStateComponent& A, Normalization n = Intensive);
MatrixOperator ExpandBasis2(MPStateComponent& A, Normalization n = Intensive);

// Constructs an 'expanded' basis given the input right hand basis
MPStateComponent ConstructFromRightBasis(BasisList const& LocalBasis,
					 VectorBasis const& Right);

// Constructs an 'expanded' basis given the input left hand basis
MPStateComponent ConstructFromLeftBasis(BasisList const& LocalBasis,
					VectorBasis const& LeftBasis);

// MPStateProject
// A functor to project onto the given component of a MPStateComponent
struct MPStateProject
{
   typedef MPStateComponent::OperatorType result_type;
   typedef MPStateComponent argument_type;
   MPStateProject(int s_) : s(s_) {}

   result_type operator()(MPStateComponent const& x) const { return x[s]; }

   int s;
};

template <typename FwdIter>
MPStateComponent tensor_col_accumulate(FwdIter first, FwdIter last, 
                                       SumBasis<VectorBasis> const& B1)
{
   //   PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());
   //   PRECONDITION_EQUAL(A.Basis2(), B.Basis2());
   typedef boost::transform_iterator<MPStateProject, FwdIter> StateIter;
   MPStateComponent Result(first->SiteBasis(), B1, first->Basis2());
   for (int s = 0; s < first->SiteBasis().size(); ++s)
   {
      Result[s] = tensor_col_accumulate(StateIter(first, MPStateProject(s)), 
					StateIter(last, MPStateProject(s)), B1);
   }
   return Result;
}

template <typename FwdIter>
MPStateComponent tensor_row_accumulate(FwdIter first, FwdIter last,
                                       SumBasis<VectorBasis> const& B2)
{
   //   PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());
   //   PRECONDITION_EQUAL(A.Basis1(), B.Basis1());
   typedef boost::transform_iterator<MPStateProject, FwdIter> StateIter;
   MPStateComponent Result(first->SiteBasis(), first->Basis1(), B2);
   for (int s = 0; s < first->size(); ++s)
   {
      Result[s] = tensor_row_accumulate(StateIter(first, MPStateProject(s)), 
					StateIter(last, MPStateProject(s)), B2);
   }
   return Result;
}

template <typename FwdIter>
MPStateComponent tensor_accumulate(FwdIter first, FwdIter last,
                                   SumBasis<VectorBasis> const& B1, 
                                   SumBasis<VectorBasis> const& B2)
{
   //   PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());
   typedef boost::transform_iterator<MPStateProject, FwdIter> StateIter;
   MPStateComponent Result(first->SiteBasis(), B1, B2);
   for (int s = 0; s < first->SiteBasis().size(); ++s)
   {
      Result[s] = tensor_accumulate(StateIter(first, MPStateProject(s)), 
				    StateIter(last, MPStateProject(s)), B1, B2);
   }
   return Result;
}

MatrixOperator MakeRandomMatrixOperator(VectorBasis const& B1, VectorBasis const& B2,
					QuantumNumber q);

inline
MatrixOperator MakeRandomMatrixOperator(VectorBasis const& B1, VectorBasis const& B2)
{
   return MakeRandomMatrixOperator(B1, B2, QuantumNumber(B1.GetSymmetryList()));
}

#include "mpstate.cc"

//
// typedef MPWavefunction
//

//typedef MatrixProduct<MPStateComponent> MPWavefunction;

#endif
