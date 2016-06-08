// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/operator_component.h
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

// New version for an MPO component.
// This version looks like a matrix of reducible operators, which is a total reversal
// of the order of the indices from the original version.
// It can also support operators that act between different local basis sets,
// via the LocalBasis1() and LocalBasis2() functions.

#if !defined(MPO_COMPONENT_H_JKVH8RY894367R789YH9P8H)
#define MPO_COMPONENT_H_JKVH8RY894367R789YH9P8H

#include "tensor/tensor.h"
#include "tensor/reducible.h"
#include "tensor/tensorsum.h"
#include "mps/state_component.h"

using Tensor::BasisList;
using Tensor::VectorBasis;
using Tensor::ReducibleTensor;
using Tensor::IrredTensor;
using Tensor::SumBasis;
using QuantumNumbers::SymmetryList;
using LinearAlgebra::size_type;

// default epsilon for detecting whether an eigenvalue is equal to 1 for
// operator classifications
double const DefaultClassifyUnityEpsilon = 1E-14;

// a simple reducible operator typedef
typedef ReducibleTensor<std::complex<double>, BasisList, BasisList> SimpleRedOperator;

typedef IrredTensor<std::complex<double> > SimpleOperator;

typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
                            VectorBasis, 
                            VectorBasis> MatrixOperator;

class OperatorComponent
{
   public:
      typedef SimpleRedOperator value_type;
      typedef value_type operator_type;

      typedef Tensor::BasisList basis1_type;
      typedef Tensor::BasisList basis2_type;

      typedef LinearAlgebra::SparseMatrix<value_type> data_type;

      typedef LinearAlgebra::iterator<data_type>::type       iterator;
      typedef LinearAlgebra::const_iterator<data_type>::type const_iterator;

      typedef LinearAlgebra::inner_iterator<data_type>::type       inner_iterator;
      typedef LinearAlgebra::const_inner_iterator<data_type>::type const_inner_iterator;

      OperatorComponent() {}
   
      // Construction of an empty OperatorComponent.  The local basis comes first.
      OperatorComponent(BasisList const& LocalB, 
			BasisList const& B1, BasisList const& B2);
      OperatorComponent(BasisList const& LocalB1, BasisList const& LocalB2, 
			BasisList const& B1, BasisList const& B2);

      SymmetryList const& GetSymmetryList() const 
      { return LocalBasis1_.GetSymmetryList(); }

      BasisList const& Basis1() const { return Basis1_; }
      BasisList const& Basis2() const { return Basis2_; }

      BasisList const& LocalBasis1() const { return LocalBasis1_; }
      BasisList const& LocalBasis2() const { return LocalBasis2_; }

      template <typename U>
      OperatorComponent operator*=(U const& x)
      { Data_ *= x; return *this; }

      OperatorComponent& operator+=(OperatorComponent const& x);
      OperatorComponent& operator-=(OperatorComponent const& x);

      // element access.  Returns a zero operator if the component does not exist.
      value_type operator()(int i, int j) const;

      // element access, inserts the element if it does not exist
      value_type& operator()(int i, int j);

      // depreciated - exists for backwards compatibility, don't use.
      void set_operator(int i, int j, value_type const& x)
      { this->operator()(i,j) = x; }

      // returns true if this matrix is in lower-triangular form.  This is only
      // useful if the matrix is square.
      bool is_lower_triangular() const;

      // returns the top-left entry, equivalent to operator()(0,0)
      value_type top_left() const;

      // returns the bottom-right entry, equivalent to operator()(size1()-1,size2()-2)
      value_type bottom_right() const;

      size_type size1() const { return Basis1_.size(); }
      size_type size2() const { return Basis2_.size(); }

      bool is_null() const { return is_zero(Data_); }

      const_inner_iterator iterate_at(int i, int j) const 
      { return LinearAlgebra::iterate_at(Data_, i, j); }

      data_type& data() { return Data_; }
      data_type const& data() const { return Data_; }

      // Makes a 1x1 MPO with the local identity operator
      static OperatorComponent make_identity(BasisList const& LocalBasis);

      // Makes a diagonl MPO with the local identity operator repeated in each position of the AuxBasis
      static OperatorComponent make_identity(BasisList const& LocalBasis, BasisList const& AuxBasis);

      void check_structure() const;
      void debug_check_structure() const;

   friend PStream::opstream& 
          operator<<(PStream::opstream& out, OperatorComponent const& Op);
   friend PStream::ipstream& 
          operator>>(PStream::ipstream& in, OperatorComponent& Op);

   private:
      BasisList LocalBasis1_, LocalBasis2_;
      BasisList Basis1_, Basis2_;
      data_type Data_;
};

std::ostream&
operator<<(std::ostream& out, OperatorComponent const& op);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(OperatorComponent const& Op, std::ostream& out, double UnityEpsilon);

inline
void print_structure(OperatorComponent const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

inline
OperatorComponent 
OperatorComponent::make_identity(BasisList const& LocalBasis)
{
   BasisList bl = Tensor::make_vacuum_basis(LocalBasis.GetSymmetryList());
   OperatorComponent Result(LocalBasis, LocalBasis, bl, bl);
   Result(0,0) = SimpleRedOperator::make_identity(LocalBasis);
   return Result;
}

inline
OperatorComponent 
OperatorComponent::make_identity(BasisList const& LocalBasis, BasisList const& AuxBasis)
{
   OperatorComponent Result(LocalBasis, LocalBasis, AuxBasis, AuxBasis);
   SimpleRedOperator I = SimpleRedOperator::make_identity(LocalBasis);
   for (unsigned i = 0; i < AuxBasis.size(); ++i)
   {
      Result(i,i) = I;
   }
   return Result;
}

// Constructs the swap gate, that maps from the tensor product basis of B1 * B2 into
// the basis B2 * B1, such that state |i,j> maps into |j,i>
SimpleOperator
swap_gate(BasisList const& B1, BasisList const& B2, 
	  ProductBasis<BasisList, BasisList> const& Basis_21,
	  ProductBasis<BasisList, BasisList> const& Basis_12);



inline
SimpleOperator
swap_gate(BasisList const& B1, BasisList const& B2)
{
   return swap_gate(B1, B2, ProductBasis<BasisList, BasisList>(B2,B1), ProductBasis<BasisList, BasisList>(B1,B2));
}

// Fermionic swap gate that inserts (-1)^{N_1 * N_2} into the swap gate.
// Parity1 and Parity2 are the diagonal components of the P operator, ie Parity[i] == 1
// iff B1(i) is bosonic, and Parity1[i] == -1 iff B1(i) is fermionic.
SimpleOperator
swap_gate_fermion(BasisList const& B1, LinearAlgebra::Vector<double> const& Parity1, 
		  BasisList const& B2, LinearAlgebra::Vector<double> const& Parity2,
		  ProductBasis<BasisList, BasisList> const& Basis_21,
		  ProductBasis<BasisList, BasisList> const& Basis_12);

inline
SimpleOperator
swap_gate(BasisList const& B1, LinearAlgebra::Vector<double> const& Parity1, 
	  BasisList const& B2, LinearAlgebra::Vector<double> const& Parity2)
{
   return swap_gate(B1, B2, ProductBasis<BasisList, BasisList>(B2,B1), ProductBasis<BasisList, BasisList>(B1,B2));
}

// Constructs an MPO that represents a shift operator
//     |
//      \
// ---   --- IncomingBasis
//    \
//     |
//    ThisBasis

OperatorComponent
translate_left(BasisList const& LeftBasis, BasisList const& ThisBasis);

// Constructs an MPO that represents a translation to the left
//               |
//              /
// LeftBasis ---   --- 
//                /
//               |
//           ThisBasis
//
// The result is an MPO with local basis (LeftBasis, ThisBasis), and
// auxiliary basis (LeftBasis, ThisBasis)

OperatorComponent
translate_right(BasisList const& LeftBasis, BasisList const& ThisBasis);

namespace LinearAlgebra
{

template <>
struct interface<OperatorComponent>
{
   typedef void type;
};

template <>
struct Iterate<OperatorComponent&>
{
   typedef OperatorComponent& argument_type;
   typedef Iterate<OperatorComponent::data_type&>::result_type result_type;

   result_type operator()(argument_type x) const
   {
      return iterate(x.data());
   }
};

template <>
struct Iterate<OperatorComponent>
{
   typedef OperatorComponent const& argument_type;
   typedef Iterate<OperatorComponent::data_type>::result_type result_type;

   result_type operator()(argument_type x) const
   {
      return iterate(x.data());
   }
};

// hermitian conjugation

template <>
struct Herm<OperatorComponent>
{
   typedef HermitianProxy<OperatorComponent> result_type;
   typedef OperatorComponent const& argument_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x);
   }
};

// conjugation

template <>
struct Conj<OperatorComponent>
{
   typedef OperatorComponent const& argument_type;
   typedef OperatorComponent result_type;

   result_type operator()(argument_type x) const 
   { 
      OperatorComponent Result(x);
      for (OperatorComponent::iterator I = iterate(Result); I; ++I)
      {
	 for (OperatorComponent::inner_iterator J = iterate(I); J; ++J)
	 {
	    *J = conj(*J);
	 }
      }
      return Result;
   }
};

// The adjoint of a component takes the adjoint of the local operator,
// hence swapping the local indices.  It does a 'flip conjugation'
// in the auxiliary indices.
template <>
struct Adjoint<OperatorComponent>
{
   typedef OperatorComponent const& argument_type;
   typedef OperatorComponent result_type;

   result_type operator()(argument_type x) const 
   { 
      OperatorComponent Result(x.LocalBasis2(), x.LocalBasis1(), adjoint(x.Basis1()), adjoint(x.Basis2()));
      for (OperatorComponent::const_iterator I = iterate(x); I; ++I)
      {
	 for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
	 {
	    Result(J.index1(), J.index2()) = adjoint(*J);
	 }
      }
      Result.debug_check_structure();
      return Result;
   }
};

template <>
struct InvAdjoint<OperatorComponent>
{
   typedef OperatorComponent const& argument_type;
   typedef OperatorComponent result_type;

   result_type operator()(argument_type x) const 
   {
      OperatorComponent Result(x.LocalBasis2(), x.LocalBasis1(), adjoint(x.Basis1()), adjoint(x.Basis2()));
      for (OperatorComponent::const_iterator I = iterate(x); I; ++I)
      {
	 for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
	 {
	    Result(J.index1(), J.index2()) = inv_adjoint(*J);
	 }
      }
      Result.debug_check_structure();
      return Result;
   }
};

} // namespace LinearAlgebra

// Constructs a MPOpComponent that represents the sum of A and B.
// The resulting state has Result'[s] = A[s] \oplus B[s]
OperatorComponent
tensor_sum(OperatorComponent const& A, OperatorComponent const& B,
           SumBasis<BasisList> const& B1, SumBasis<BasisList> const& B2);
	   
// Constructs a MPOpComponent that represents the sum of A and B,
// at the left boundary of the matrix product state.
// Precondition: A.Basis1() == B.Basis1()
// The resulting state has Result'[s] = (A[s], B[s])  (row-wise concatenation)
OperatorComponent 
tensor_row_sum(OperatorComponent const& A, 
	       OperatorComponent const& B, 
	       SumBasis<BasisList> const& B2);

OperatorComponent 
tensor_row_sum(OperatorComponent const& A, 
	       OperatorComponent const& B);


// Constructs a MPOpComponent that represents the sum of A and B,
// at the right boundary of the matrix product state.
// Precondition: A.Basis2() == B.Basis2()
// The resulting state has Result'[s] = ( A[s] )
//                                      ( B[s] )  (column-wise concatenation)
OperatorComponent 
tensor_col_sum(OperatorComponent const& A, 
	       OperatorComponent const& B, 
	       SumBasis<BasisList> const& B1);

OperatorComponent 
tensor_col_sum(OperatorComponent const& A, 
	       OperatorComponent const& B);

// Multiplies the component by a SimpleOperator acting on the auxiliary space
OperatorComponent prod(OperatorComponent const& A, SimpleOperator const& Op);
OperatorComponent prod(SimpleOperator const& Op, OperatorComponent const& A);

OperatorComponent prod(OperatorComponent const& A, LinearAlgebra::HermitianProxy<SimpleOperator> const& Op);
OperatorComponent prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& Op, OperatorComponent const& A);

OperatorComponent triple_prod(SimpleOperator const& x, OperatorComponent const& Op, LinearAlgebra::HermitianProxy<SimpleOperator> const& y);

// constructs the tensor product in the local basis, the matrices
// C[(s's),(t't)] = A[(s's)] * B[(t't)]
// This is a coarse graining operation.  Ordinary product in the aux basis, tensor product in the local basis
OperatorComponent local_tensor_prod(OperatorComponent const& A, OperatorComponent const& B);

// constructs the tensor product in the auxiliary basis, and ordinary product in the local basis.
// This is the operation that gives the product of MPO's.
// This used to be called mp_prod
OperatorComponent aux_tensor_prod(OperatorComponent const& A, OperatorComponent const& B);

// constructs the tensor product in the auxiliary basis, and ordinary product in the local basis.
// This gives the action of an MPO on an MPS
// This used to be called mp_prod
StateComponent aux_tensor_prod(OperatorComponent const& A, StateComponent const& B);

// construct the tensor product in both local and aux basis
OperatorComponent global_tensor_prod(OperatorComponent const& A, OperatorComponent const& B);

// Constructs the matrix of overlaps Result'_{ij} = inner_prod(A(i,k), B(k,j))
SimpleOperator
local_inner_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B);

SimpleOperator
local_inner_prod(OperatorComponent const& A, HermitianProxy<OperatorComponent> const& B);

// constructs the tensor product matrix Result((i',j'), (i,j)) = inner(A(i',i), B(j',j))
// This is the matrix representation of the MPO transfer matrix
SimpleOperator
local_inner_tensor_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B);

// Performs the adjoint operation on the local operators.  The auxiliary basis
// transforms into the adjoint basis but preserves the basis1 / basis2 relationship
OperatorComponent
local_adjoint(OperatorComponent const& A);

// Flip-conjugates all bases
OperatorComponent flip_conj(OperatorComponent const& A);

// exchanges the local and auxiliary components of the MPO.  This corresponds to a
// transpose operation on the 4-index tensor.  
OperatorComponent exchange(OperatorComponent const& A);

// Basis truncation via removal of parallel components.
// Returns the operator P such that P * A' = A
// This function preserves the last row of the MPO intact
// (which is required to preserve upper triangular MPO's)
SimpleOperator
TruncateBasis1(OperatorComponent& A);

SimpleOperator
TruncateBasis1MkII(OperatorComponent& A, double Epsilon = 0.0);

// Basis truncation via removal of parallel components.
// Returns the operator P such that A' * P = A
// This function preserves the first column of the MPO intact
// (which is required to preserve upper triangular MPO's)
SimpleOperator
TruncateBasis2(OperatorComponent& A);

SimpleOperator
TruncateBasis2MkII(OperatorComponent& A, double Epsilon = 0.0);

OperatorComponent 
operator+(OperatorComponent const& A, OperatorComponent const& Op);

OperatorComponent 
operator-(OperatorComponent const& A, OperatorComponent const& Op);

#if 0
// dont need these, linearalgebra::multiplication<> takes care
inline
OperatorComponent 
operator*(OperatorComponent const& A, LinearAlgebra::HermitianProxy<SimpleOperator> const& Op)
{
   return prod(A, Op);
}

inline
OperatorComponent 
operator*(LinearAlgebra::HermitianProxy<SimpleOperator> const& Op, OperatorComponent const& A)
{
   return prod(Op, A);
}
#endif

inline
OperatorComponent
operator*(double x, OperatorComponent const& Op)
{
   OperatorComponent Result(Op);
   Result *= x;
   return Result;
}

inline
OperatorComponent
operator*(std::complex<double> x, OperatorComponent const& Op)
{
   OperatorComponent Result(Op);
   Result *= x;
   return Result;
}

OperatorComponent
conj(OperatorComponent const& x);

// This grouping treats A as a bra, and herm(B) as a ket

// Contraction from the right
// does Result'[a'](i',j') = M(s',s)(a',a) A[s'](i',i) E[a](i,j) herm(B[s](j',j))
// This is an `expectation value', the transfer operator applied to a state
// ** this function is deprecated -- use contract_from_right instead
#if defined(OLD_OPERATOR_PROD)
StateComponent
operator_prod(OperatorComponent const& M,
              StateComponent const& A, 
              StateComponent const& E,
              LinearAlgebra::HermitianProxy<StateComponent> const& B);
#endif

// F-matrix contraction
//
// Contracting from the right hand side.
// The MPO element is compex-conjugated because the F object will finally
// be conjugated itself (making A bra and herm(B) ket).
//
//                  i'  i
//                 --A--+
//                   |  |
//F'[a'](i'j') = a'--M*-F
//                   |  |
//                 --B*-+
//                  j'  j
//
// Result'[a'](i',j') = herm(M(s',s)(a',a)) A[s'](i',i) F[a](i,j) herm(B[s](j',j))
StateComponent
contract_from_right(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
		    StateComponent const& A, 
		    StateComponent const& F,
		    LinearAlgebra::HermitianProxy<StateComponent> const& B);

inline
StateComponent
contract_from_right(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
		    SimpleStateComponent const& A, 
		    StateComponent const& F,
		    LinearAlgebra::HermitianProxy<SimpleStateComponent> const& B)
{
   // TODO: optimize this implementation
   StateComponent AX = A;
   StateComponent BX = B.base();
   return contract_from_right(M, AX, F, herm(BX));
}

StateComponent
contract_from_right_mask(HermitianProxy<OperatorComponent> const& M,
			 StateComponent const& A, 
			 StateComponent const& F,
			 HermitianProxy<StateComponent> const& B,
			 std::vector<int> const& Mask1,
			 std::vector<int> const& Mask2);

// Contraction from the left
// Result'[a](i,j) = herm(M(s',s)(a',a)) herm(A[s'](i',i)) E[a'](i',j') B[s](j',j)
// ** this function is deprecated -- use contract_from_left instead
#if defined(OLD_OPERATOR_PROD)
StateComponent
operator_prod(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
              LinearAlgebra::HermitianProxy<StateComponent> const& A, 
              StateComponent const& F,
              StateComponent const& B);
#endif

// E-matrix contraction
//
// +--A*-- i
// |  |
// E--M--- a = E'[a](i,j)
// |  |
// +--B--- j
//
// Result'[a](i,j) = M(s',s)(a',a) herm(A[s'](i',i)) E[a'](i',j') B[s](j',j)
StateComponent
contract_from_left(OperatorComponent const& M,
		   LinearAlgebra::HermitianProxy<StateComponent> const& A, 
		   StateComponent const& E,
		   StateComponent const& B);


// Action of an operator on B
// Result[s'](i',i) = M(s',s)[a',a] E[a'](i',j') B[s](j',j) herm(F[a](i,j))
StateComponent
operator_prod_inner(OperatorComponent const& M,
                    StateComponent const& E, 
                    StateComponent const& B,
                    LinearAlgebra::HermitianProxy<StateComponent> const& F);

// Variants for a regular triangular MPO

// Precondition: The MPS is normalized, such that scalar_prod(A, herm(B)) = 1,
// F.back() = 1, M is upper triangular normal form.
StateComponent
operator_prod_regular(OperatorComponent const& M,
		      StateComponent const& A, 
		      StateComponent const& F,
		      LinearAlgebra::HermitianProxy<StateComponent> const& B);

// Precondition: The MPS is normalized, such that scalar_prod(herm(A), B) = 1,
// E.front = 1, M is upper triangular normal form.
StateComponent
operator_prod_regular(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
		      LinearAlgebra::HermitianProxy<StateComponent> const& A, 
		      StateComponent const& E,
		      StateComponent const& B);

#if 0
// does Result'[s'](i',j') = M(s',s)(a',a) E[a'](i',i)^* [s](i,j) F[a](j',j)^*
// This is the action of an operator on a state E A herm(F)
StateComponent
operator_prod_inner(OperatorComponent const& M,
                    StateComponent const& E, 
                    StateComponent const& A,
                    LinearAlgebra::HermitianProxy<StateComponent> const& F);

StateComponent
operator_prod_inner(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
                    LinearAlgebra::HermitianProxy<StateComponent> const& E, 
                    StateComponent const& A,
                    StateComponent const& F);
#endif

double
norm_frob_sq(OperatorComponent const& x);

inline
double
norm_frob(OperatorComponent const& x)
{
   return std::sqrt(norm_frob_sq(x));
}

// returns the specified row of the MPO
LinearAlgebra::MapVector<SimpleRedOperator>
project_row(OperatorComponent const& x, int r);

// returns the specified column of the MPO
LinearAlgebra::MapVector<SimpleRedOperator>
project_column(OperatorComponent const& x, int c);

// project onto the given rows of the component
OperatorComponent
project_rows(OperatorComponent const& x, std::set<int> const& Rows);

// project onto the given columns of the component
OperatorComponent
project_columns(OperatorComponent const& x, std::set<int> const& Cols);

// Utility function to decompose a SimpleOperator defined over a tensor product basis
// into a sum of products of operators over the original basis.
// This is a fine-graining operation on the operator,
// the inverse of the tensor_prod
std::vector<std::pair<SimpleRedOperator, SimpleRedOperator> >
decompose_tensor_prod(SimpleOperator const& Op, 
                      ProductBasis<BasisList, BasisList> const& B1,
                      ProductBasis<BasisList, BasisList> const& B2);

// fine-grain an OperatorComponent, the inverse of local_tensor_prod
std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(OperatorComponent const& Op, 
			    ProductBasis<BasisList, BasisList> const& B1,
			    ProductBasis<BasisList, BasisList> const& B2);

// Decompose a SimpleOperator into an MPO, with vacuum on the right basis,
// single on the left basis.
std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(SimpleOperator const& Op, 
			    ProductBasis<BasisList, BasisList> const& B1,
			    ProductBasis<BasisList, BasisList> const& B2);

// Version for when the Basis1 and Basis2 of the SimpleOperator are the same
std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(SimpleOperator const& Op, 
			    ProductBasis<BasisList, BasisList> const& B);

// Version for when the Basis1 and Basis2 of the SimpleOperator are the same,
// and the product basis is a simple cartesian product of two BasisList's.
// BasisA and BasisB will be the local_basis of the two result OperatorComponent's
std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(SimpleOperator const& Op, 
			    BasisList const& BasisA, BasisList const& BasisB);

// Take a state component A^s(i,j) and rotate it into an operator component
// M(s,0)(i,j), where (i,j) is the new local basis and the right boundary
// is closed with an identity rep.
OperatorComponent
RotateToOperatorRightBoundary(StateComponent const& x);

// Take a state component A^s(i,j) and rotate it into an operator component
// M(0,s)(\bar{i},\bar{j}), where (\bar{i},\bar{j}) is the new local basis and the left boundary
// is closed with an identity rep.  The components of M(0,s) are the flip conjugation of A^s.
OperatorComponent
RotateToOperatorLeftBoundary(StateComponent const& x);

// update_mask_basis functions.  Given a 'mask' of used basis elements for one basis,
// update a mask for the other basis.  The mask is a vector of boolean (vector<int>, to 
// avoid the vector<bool> bogosity)
void
update_mask_basis1(std::vector<int>& Mask1, SimpleOperator const& Op, std::vector<int> const& Mask2);

void
update_mask_basis1(std::vector<int>& Mask1, SimpleRedOperator const& Op, std::vector<int> const& Mask2);

void
update_mask_basis1(std::vector<int>& Mask1, OperatorComponent const& Op, std::vector<int> const& Mask2);


void
update_mask_basis2(std::vector<int> const& Mask1, SimpleOperator const& Op, std::vector<int>& Mask2);

void
update_mask_basis2(std::vector<int> const& Mask1, SimpleRedOperator const& Op, std::vector<int>& Mask2);

void
update_mask_basis2(std::vector<int> const& Mask1, OperatorComponent const& Op, std::vector<int>& Mask2);

MatrixOperator ExpandBasis1Used(StateComponent& A, OperatorComponent const& Op);
MatrixOperator ExpandBasis2Used(StateComponent& A, OperatorComponent const& Op);

// helper function to construct an operator that merges repeated
// quantum numbers in the basis.  
SimpleOperator CollapseBasis(BasisList const& b);

// helper function to construct an operator that projects onto a single quantum number component
SimpleOperator ProjectBasis(BasisList const& b, QuantumNumbers::QuantumNumber const& q);

// Utility function - if X is proportional to identity operator then return the constant of
// proportionality.  Otherwise return 0.0
std::complex<double> PropIdent(SimpleOperator const& X, double UnityEpsilon);

inline
std::complex<double> PropIdent(SimpleOperator const& X)
{
   return PropIdent(X, DefaultClassifyUnityEpsilon);
}

namespace LinearAlgebra
{

template <>
struct Multiplication<SimpleOperator, OperatorComponent>
{
   typedef SimpleOperator first_argument_type;
   typedef OperatorComponent second_argument_type;
   typedef OperatorComponent result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      return prod(x, y);
   }
};

template <>
struct Multiplication<OperatorComponent, SimpleOperator>
{
   typedef OperatorComponent first_argument_type;
   typedef SimpleOperator second_argument_type;
   typedef OperatorComponent result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      return prod(x, y);
   }
};

template <>
struct Multiplication<HermitianProxy<SimpleOperator>, OperatorComponent>
{
   typedef HermitianProxy<SimpleOperator> first_argument_type;
   typedef OperatorComponent second_argument_type;
   typedef OperatorComponent result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      return prod(x, y);
   }
};

template <>
struct Multiplication<OperatorComponent, HermitianProxy<SimpleOperator> >
{
   typedef OperatorComponent first_argument_type;
   typedef HermitianProxy<SimpleOperator> second_argument_type;
   typedef OperatorComponent result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      return prod(x, y);
   }
};

} // namespace LinearAlgebra


#include "operator_component.cc"

#endif
