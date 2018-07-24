// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/operator_component.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_MPO_OPERATOR_COMPONENT_H)
#define MPTOOLKIT_MPO_OPERATOR_COMPONENT_H

#include "tensor/tensor.h"
#include "tensor/reducible.h"
#include "tensor/tensorsum.h"
#include "mps/state_component.h"
#include "operator_utilities.h"

using Tensor::BasisList;
using Tensor::VectorBasis;
using Tensor::ReducibleTensor;
using Tensor::IrredTensor;
using Tensor::SumBasis;
using Tensor::HermitianProxy;
using QuantumNumbers::SymmetryList;

template <typename T>
class BasicOperatorComponent;

template <typename T>
PStream::opstream&
operator<<(PStream::opstream& out, BasicOperatorComponent<T> const& Op);

template <typename T>
PStream::ipstream&
operator>>(PStream::ipstream& in, BasicOperatorComponent<T>& Op);

template <typename T>
BasicOperatorComponent<T>
copy(BasicOperatorComponent<T> const& x);

template <typename T>
class BasicOperatorComponent
{
   public:
      using numeric_type  = T;
      using value_type    = ReducibleTensor<numeric_type, BasisList, BasisList>;
      using operator_type = value_type;

      using basis1_type = Tensor::BasisList;
      using basis2_type = Tensor::BasisList;

      using data_type      = blas::SparseMatrix<value_type>;
      using row_type       = typename data_type::row_type;
      using iterator       = typename data_type::iterator;
      using const_iterator = typename data_type::const_iterator;

      static_assert(std::is_nothrow_move_constructible<data_type>::value, "");
      static_assert(std::is_nothrow_destructible<data_type>::value, "");
      static_assert(std::is_destructible<data_type>::value, "");

      BasicOperatorComponent() noexcept = default;

      // Construction of an empty BasicOperatorComponent.  The local basis comes first.
      BasicOperatorComponent(BasisList const& LocalB,
                        BasisList const& B1, BasisList const& B2);
      BasicOperatorComponent(BasisList const& LocalB1, BasisList const& LocalB2,
                        BasisList const& B1, BasisList const& B2);

      BasicOperatorComponent(BasicOperatorComponent const& Other)
	 : LocalBasis1_(Other.LocalBasis1_),
	   LocalBasis2_(Other.LocalBasis2_),
	   Basis1_(Other.Basis1_),
	   Basis2_(Other.Basis2_),
	   Data_(copy(Other.Data_))
      {}

      BasicOperatorComponent(BasicOperatorComponent&& Other) noexcept = default;

      ~BasicOperatorComponent() noexcept = default;

      BasicOperatorComponent& operator=(BasicOperatorComponent&& Other) noexcept = default;

      BasicOperatorComponent& operator=(BasicOperatorComponent const& Other)
      {
	 LocalBasis1_ = Other.LocalBasis1_;
	 LocalBasis2_ = Other.LocalBasis2_;
	 Basis1_ = Other.Basis1_;
	 Basis2_ = Other.Basis2_;
	 Data_ = copy(Other.Data_);
	 return *this;
      }

      SymmetryList const& GetSymmetryList() const
      { return LocalBasis1_.GetSymmetryList(); }

      BasisList const& Basis1() const { return Basis1_; }
      BasisList const& Basis2() const { return Basis2_; }

      BasisList const& LocalBasis1() const { return LocalBasis1_; }
      BasisList const& LocalBasis2() const { return LocalBasis2_; }

      template <typename U>
      BasicOperatorComponent& operator*=(U const& x)
      { Data_ *= x; return *this; }

      BasicOperatorComponent& operator+=(BasicOperatorComponent const& x);
      BasicOperatorComponent& operator-=(BasicOperatorComponent const& x);

      // element access.  Returns a zero operator if the component does not exist.
      value_type operator()(int i, int j) const;

      // returns true if this matrix is in lower-triangular form.  This is only
      // useful if the matrix is square.
      bool is_lower_triangular() const;

      // returns the top-left entry, equivalent to operator()(0,0)
      value_type top_left() const;

      // returns the bottom-right entry, equivalent to operator()(size1()-1,size2()-2)
      value_type bottom_right() const;

      int size1() const { return Basis1_.size(); }
      int size2() const { return Basis2_.size(); }

      bool is_null() const { return Data_.empty(); }

      iterator begin() { return Data_.begin(); }
      iterator end() { return Data_.end(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }

      const_iterator cbegin() const { return Data_.cbegin(); }
      const_iterator cend() const { return Data_.cend(); }

      row_type& operator[](int r) { return Data_[r]; }
      row_type const& operator[](int r) const { return Data_[r]; }

      row_type& row(int r) { return Data_.row(r); }
      row_type const& row(int r) const { return Data_.row(r); }

      bool exists(int r, int c) const
      {
	 return Data_.row(r).exists(c);
      }

      void erase(int r, int c)
      {
	 Data_.row(r).erase(c);
      }

      template<typename... Args>
      void emplace(int Row, int Col, Args&&... args)
      {
	 Data_.emplace(Row, Col, std::forward<Args>(args)...);
      }

      template <typename U>
      void insert(int r, int c, U const& value)
      {
	 Data_.insert(r, c, value);
      }

      template <typename U>
      void insert(int r, int c, U&& value)
      {
	 Data_.insert(r, c, std::move(value));
      }

      template <typename U>
      void add(int r, int c, U const& value)
      {
	 Data_.add(r, c, value);
      }

      template <typename U>
      void add(int r, int c, U&& value)
      {
	 Data_.add(r, c, std::move(value));
      }

      template <typename U>
      void subtract(int r, int c, U const& value)
      {
	 Data_.subtract(r, c, value);
      }

      template <typename U>
      void subtract(int r, int c, U&& value)
      {
	 Data_.subtract(r, c, std::move(value));
      }

      QuantumNumber const& qn1(int i) const { return Basis1_[i]; }
      QuantumNumber const& qn2(int j) const { return Basis2_[j]; }

      data_type& data() { return Data_; }
      data_type const& data() const { return Data_; }

      // Makes a 1x1 MPO with the local identity operator
      static BasicOperatorComponent make_identity(BasisList const& LocalBasis);

      // Makes a diagonl MPO with the local identity operator repeated in each position of the AuxBasis
      static BasicOperatorComponent make_identity(BasisList const& LocalBasis, BasisList const& AuxBasis);

      void check_structure() const;
      void debug_check_structure() const;

      friend PStream::opstream&
      operator<< <T>(PStream::opstream& out, BasicOperatorComponent<T> const& Op);

      friend PStream::ipstream&
      operator>> <T>(PStream::ipstream& in, BasicOperatorComponent<T>& Op);

      friend BasicOperatorComponent copy<T>(BasicOperatorComponent const& x);

      static_assert(std::is_nothrow_move_constructible<BasisList>::value, "");
      static_assert(std::is_nothrow_move_constructible<data_type>::value, "");
      static_assert(std::is_nothrow_destructible<BasisList>::value, "");
      static_assert(std::is_nothrow_destructible<data_type>::value, "");

      BasicOperatorComponent(BasisList const& LocalB1, BasisList const& LocalB2,
			     BasisList const& B1, BasisList const& B2,
			     data_type const& Data)
	 : LocalBasis1_(LocalB1), LocalBasis2_(LocalB2), Basis1_(B1), Basis2_(B2),
	   Data_(copy(Data)) {}

   private:
      BasisList LocalBasis1_;
      BasisList LocalBasis2_;
      BasisList Basis1_;
      BasisList Basis2_;
      data_type Data_;
};

using OperatorComponent = BasicOperatorComponent<complex>;

static_assert(std::is_nothrow_destructible<OperatorComponent>::value, "");
static_assert(std::is_destructible<OperatorComponent>::value, "");
static_assert(std::is_nothrow_move_constructible<OperatorComponent>::value, "");


template <typename T>
BasicOperatorComponent<T>
copy(BasicOperatorComponent<T> const& x)
{
   return BasicOperatorComponent<T>(x.LocalBasis1(), x.LocalBasis2(), x.Basis1(), x.Basis2(), x.Data_);
}

template <typename T>
inline
HermitianProxy<BasicOperatorComponent<T>>
herm(BasicOperatorComponent<T> const& x)
{
   return HermitianProxy<BasicOperatorComponent<T>>(x);
}


template <typename T>
std::ostream&
operator<<(std::ostream& out, BasicOperatorComponent<T> const& op);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
template <typename T>
void print_structure(BasicOperatorComponent<T> const& Op, std::ostream& out, double UnityEpsilon);

template <typename T>
inline
void print_structure(BasicOperatorComponent<T> const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon<T>);
}

template <typename T>
inline
BasicOperatorComponent<T>
BasicOperatorComponent<T>::make_identity(BasisList const& LocalBasis)
{
   BasisList bl = Tensor::make_vacuum_basis(LocalBasis.GetSymmetryList());
   BasicOperatorComponent<T> Result(LocalBasis, LocalBasis, bl, bl);
   Result.insert(0,0, SimpleRedOperator_t<T>::make_identity(LocalBasis));
   return Result;
}

template <typename T>
inline
BasicOperatorComponent<T>
BasicOperatorComponent<T>::make_identity(BasisList const& LocalBasis, BasisList const& AuxBasis)
{
   BasicOperatorComponent<T> Result(LocalBasis, LocalBasis, AuxBasis, AuxBasis);
   auto I = SimpleRedOperator_t<T>::make_identity(LocalBasis);
   for (unsigned i = 0; i < AuxBasis.size(); ++i)
   {
      Result.insert(i,i, copy(I));
   }
   return Result;
}

// Constructs an MPO that represents a shift operator
//     |
//      \    o
// ---   --- IncomingBasis
//    \      o
//     |
//    ThisBasis

OperatorComponent
translate_left(BasisList const& LeftBasis, BasisList const& ThisBasis);

// Constructs an MPO that represents a translation to the left
//               |
//              /       o
// LeftBasis ---   ---  o
//                /     o
//               |
//           ThisBasis
//
// The result is an MPO with local basis (LeftBasis, ThisBasis), and
// auxiliary basis (LeftBasis, ThisBasis)

OperatorComponent
translate_right(BasisList const& LeftBasis, BasisList const& ThisBasis);

// hermitian conjugation

inline
OperatorComponent adjoint(OperatorComponent const& x)
{
   //   static_assert(std::is_nothrow_move_constructible<OperatorComponent>::value, "");
   OperatorComponent Result(x.LocalBasis2(), x.LocalBasis1(), adjoint(x.Basis1()), adjoint(x.Basis2()));
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
	 Result.insert(r.row(), c.col(), adjoint(c.value));
      }
   }
   Result.debug_check_structure();
   return Result;
}

inline
OperatorComponent inv_adjoint(OperatorComponent const& x)
{
   OperatorComponent Result(x.LocalBasis2(), x.LocalBasis1(), adjoint(x.Basis1()), adjoint(x.Basis2()));
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
	 Result.insert(r.row(), c.col(), inv_adjoint(c.value));
      }
   }
   Result.debug_check_structure();
   return Result;
}

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
OperatorComponent prod(OperatorComponent const& A, SimpleOperator const& Op, double Tol = 1e-14);
OperatorComponent prod(SimpleOperator const& Op, OperatorComponent const& A, double Tol = 1e-14);

OperatorComponent prod(OperatorComponent const& A, HermitianProxy<SimpleOperator> const& Op, double Tol = 1e-14);
OperatorComponent prod(HermitianProxy<SimpleOperator> const& Op, OperatorComponent const& A, double Tol = 1e-14);

inline
OperatorComponent
operator*(OperatorComponent const& A, SimpleOperator const& Op)
{
   return prod(A,Op);
}

inline
OperatorComponent
operator*(SimpleOperator const& Op, OperatorComponent const& A)
{
   return prod(Op, A);
}

inline
OperatorComponent
operator*(OperatorComponent const& A, HermitianProxy<SimpleOperator> const& Op)
{
   return prod(A, Op);
}

inline
OperatorComponent
operator*(HermitianProxy<SimpleOperator> const& Op, OperatorComponent const& A)
{
   return prod(Op, A);
}

OperatorComponent triple_prod(SimpleOperator const& x, OperatorComponent const& Op, HermitianProxy<SimpleOperator> const& y);

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

#if 0
// UNUSED
// constructs the matrix of the squared norm of the local operators.  This is NOT the same
// operation as local_inner_prod(herm(A), A), but is rather
// Result'(i,j) = norm_frob_sq(A(i,j))
RealSimpleOperator
local_norm_frob_sq(OperatorComponent const& A);
#endif

// constructs the tensor product matrix Result((i',j'), (i,j)) = inner(A(i',i), B(j',j))
// This is the matrix representation of the MPO transfer matrix
SimpleOperator
local_inner_tensor_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B);

#if 0
// UNUSED
// Performs the adjoint operation on the local operators.  The auxiliary basis
// transforms into the adjoint basis but preserves the basis1 / basis2 relationship
OperatorComponent
local_adjoint(OperatorComponent const& A);
#endif

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

template <typename T>
BasicOperatorComponent<T>
operator-(BasicOperatorComponent<T> const& A, BasicOperatorComponent<T> const& Op);

inline
OperatorComponent
operator*(double x, OperatorComponent&& Result)
{
   Result *= x;
   return std::move(Result);
}

inline
OperatorComponent
operator*(double a, OperatorComponent const& x)
{
   OperatorComponent Result(copy(x));
   Result *= a;
   return Result;
}

inline
OperatorComponent
operator*(std::complex<double> a, OperatorComponent&& Result)
{
   Result *= a;
   return std::move(Result);
}

inline
OperatorComponent
operator*(std::complex<double> a, OperatorComponent const& x)
{
   OperatorComponent Result(copy(x));
   Result *= a;
   return Result;
}

OperatorComponent
conj(OperatorComponent x);

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
              HermitianProxy<StateComponent> const& B);
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
contract_from_right(HermitianProxy<OperatorComponent> const& M,
                    StateComponent const& A,
                    StateComponent const& F,
                    HermitianProxy<StateComponent> const& B);

#if 0
inline
StateComponent
contract_from_right(HermitianProxy<OperatorComponent> const& M,
                    SimpleStateComponent const& A,
                    StateComponent const& F,
                    HermitianProxy<StateComponent> const& B)
{
   // TODO: optimize this implementation
   StateComponent AX = copy(A);
   StateComponent BX = copy(B.base());
   return contract_from_right(M, AX, F, herm(BX));
}
#endif

// version of contract_from_right where we exclude elements
// in Mask1/Mask2
StateComponent
contract_from_right_mask(HermitianProxy<OperatorComponent> const& M,
                         StateComponent const& A,
                         StateComponent const& F,
                         HermitianProxy<StateComponent> const& B,
                         std::set<int> const& Mask1,
                         std::set<int> const& Mask2);

// Contraction from the left
// Result'[a](i,j) = herm(M(s',s)(a',a)) herm(A[s'](i',i)) E[a'](i',j') B[s](j',j)
// ** this function is deprecated -- use contract_from_left instead
#if defined(OLD_OPERATOR_PROD)
StateComponent
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
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
// This function is defined in e-optim.cpp
StateComponent
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   StateComponent const& E,
                   StateComponent const& B);

// In this version, the MPO must have dimension 1, and we call using a
// MatrixOperator instead
// TODO: this isn't so efficient
inline
StateComponent
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   MatrixOperator const& E,
                   StateComponent const& B)
{
   StateComponent EE(M.Basis1(), E.Basis1(), E.Basis2());
   EE[0] = copy(E);
   return contract_from_left(M, A, EE, B);
}

// Action of an operator on B
// Result[s'](i',i) = M(s',s)[a',a] E[a'](i',j') B[s](j',j) herm(F[a](i,j))
StateComponent
operator_prod_inner(OperatorComponent const& M,
                    StateComponent const& E,
                    StateComponent const& B,
                    HermitianProxy<StateComponent> const& F);

// obtains the diagonal component of an operator
// Result[s](i',i) = M(s,s)[a',a] E[a'](i',i') herm(F[a](i,i))
StateComponent
operator_prod_inner_diagonal(OperatorComponent const& M,
			     StateComponent const& E,
			     HermitianProxy<StateComponent> const& F);

// Variants for a regular triangular MPO

// Precondition: The MPS is normalized, such that scalar_prod(A, herm(B)) = 1,
// F.back() = 1, M is upper triangular normal form.
StateComponent
operator_prod_regular(OperatorComponent const& M,
                      StateComponent const& A,
                      StateComponent const& F,
                      HermitianProxy<StateComponent> const& B);

// Precondition: The MPS is normalized, such that scalar_prod(herm(A), B) = 1,
// E.front = 1, M is upper triangular normal form.
StateComponent
operator_prod_regular(HermitianProxy<OperatorComponent> const& M,
                      HermitianProxy<StateComponent> const& A,
                      StateComponent const& E,
                      StateComponent const& B);

#if 0
// does Result'[s'](i',j') = M(s',s)(a',a) E[a'](i',i)^* [s](i,j) F[a](j',j)^*
// This is the action of an operator on a state E A herm(F)
StateComponent
operator_prod_inner(OperatorComponent const& M,
                    StateComponent const& E,
                    StateComponent const& A,
                    HermitianProxy<StateComponent> const& F);

StateComponent
operator_prod_inner(HermitianProxy<OperatorComponent> const& M,
                    HermitianProxy<StateComponent> const& E,
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

#include "operator_component.icc"

#endif
