// -*- C++ -*- $Id$

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

      static OperatorComponent make_identity(BasisList const& LocalBasis);

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

inline
OperatorComponent 
OperatorComponent::make_identity(BasisList const& LocalBasis)
{
   BasisList bl = Tensor::make_vacuum_basis(LocalBasis.GetSymmetryList());
   OperatorComponent Result(LocalBasis, LocalBasis, bl, bl);
   Result(0,0) = SimpleRedOperator::make_identity(LocalBasis);
   return Result;
}

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

#if 0
// do we need this ever?
// scalar_prod
// does Result' = sum_s A[s] * herm(B[s])

struct ScalarProd<OperatorComponent, HermitianProxy<OperatorComponent> >
{
   typedef SimpleOperator result_type;
   
};
#endif

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

// Constructs a MPOpComponent that represents the sum of A and B,
// at the right boundary of the matrix product state.
// Precondition: A.Basis2() == B.Basis2()
// The resulting state has Result'[s] = ( A[s] )
//                                      ( B[s] )  (column-wise concatenation)
OperatorComponent 
tensor_col_sum(OperatorComponent const& A, 
	       OperatorComponent const& B, 
	       SumBasis<BasisList> const& B1);

OperatorComponent prod(OperatorComponent const& A, SimpleOperator const& Op);
OperatorComponent prod(SimpleOperator const& Op, OperatorComponent const& A);

OperatorComponent prod(OperatorComponent const& A, LinearAlgebra::HermitianProxy<SimpleOperator> const& Op);
OperatorComponent prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& Op, OperatorComponent const& A);

OperatorComponent triple_prod(SimpleOperator const& x, OperatorComponent const& Op, LinearAlgebra::HermitianProxy<SimpleOperator> const& y);

// constructs the tensor product in the local basis, the matrices
// C[(s's),(t't)] = A[(s's)] * B[(t't)]
// This is a coarse graining operation.
OperatorComponent local_tensor_prod(OperatorComponent const& A, OperatorComponent const& B);

// constructs the tensor product in the auxiliary basis, and ordinary product in the local basis.
// This is the operation that gives the product of MPO's.
// This used to be called mp_prod
OperatorComponent aux_tensor_prod(OperatorComponent const& A, OperatorComponent const& B);

// construct the tensor product in both local and aux basis
OperatorComponent global_tensor_prod(OperatorComponent const& A, OperatorComponent const& B);

// Flip-conjugates all bases
OperatorComponent flip_conj(OperatorComponent const& A);

// exchanges the local and auxiliary components of the MPO.  This corresponds to a
// transpose operation on the 4-index tensor.  
OperatorComponent exchange(OperatorComponent const& A);

OperatorComponent 
operator+(OperatorComponent const& A, OperatorComponent const& Op);

inline
OperatorComponent 
operator*(OperatorComponent const& A, SimpleOperator const& Op)
{
   return prod(A, Op);
}

inline
OperatorComponent 
operator*(SimpleOperator const& Op, OperatorComponent const& A)
{
   return prod(Op, A);
}

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

OperatorComponent
conj(OperatorComponent const& x);

// does Result'(a)(i',j') = M(s',s)(a,b) A(s')(i',i) E(b)(i',j) B(s)(j',j)^*
// This is an `expectation value', the transfer operator applied to a state
StateComponent
operator_prod(OperatorComponent const& M,
              StateComponent const& A, 
              StateComponent const& F,
              LinearAlgebra::HermitianProxy<StateComponent> const& B);

StateComponent
operator_prod(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
              LinearAlgebra::HermitianProxy<StateComponent> const& A, 
              StateComponent const& E,
              StateComponent const& B);

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


// does Result'(s')(i',j') = M(s',s)(a,b) A(a)(i',i) E(s)(i',j) B(b)(j',j)^*
// This is the action of an operator on a state
StateComponent
operator_prod_inner(OperatorComponent const& M,
                    StateComponent const& A, 
                    StateComponent const& F,
                    LinearAlgebra::HermitianProxy<StateComponent> const& B);

StateComponent
operator_prod_inner(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
                    LinearAlgebra::HermitianProxy<StateComponent> const& A, 
                    StateComponent const& E,
                    StateComponent const& B);

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

// Given a scalar SimpleOperator defined over a tensor product basis,
// decompose it into a two-site MPO.  This is a kind of fine-graining operation
std::pair<OperatorComponent, OperatorComponent>
decompose_tensor_prod(SimpleOperator const& Op, 
                      ProductBasis<BasisList, BasisList> const& B1,
                      ProductBasis<BasisList, BasisList> const& B2);

// Version for when the Basis1 and Basis2 of the SimpleOperator are the same
std::pair<OperatorComponent, OperatorComponent>
decompose_tensor_prod(SimpleOperator const& Op, 
                      ProductBasis<BasisList, BasisList> const& B);

// Version for when the Basis1 and Basis2 of the SimpleOperator are the same,
// and the product basis is a simple cartesian product of two BasisList's.
// BasisA and BasisB will be the local_basis of the two result OperatorComponent's
std::pair<OperatorComponent, OperatorComponent>
decompose_tensor_prod(SimpleOperator const& Op, 
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

#include "operator_component.cc"

#endif
