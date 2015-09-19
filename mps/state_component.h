// -*- C++ -*- $Id$
//
// I prefer the naming 'local basis' over the original 'site basis'.  I intend to eventually
// add functions with 'site' -> 'local' everywhere and eventually remove the 'site' variants.

#if !defined(MPSTATE_H_SDHCKJHKJLRHGIURHYULUHR)
#define MPSTATE_H_SDHCKJHKJLRHGIURHYULUHR

#include "tensor/tensor.h"
#include "tensor/reducible.h"
#include "tensor/tensorsum.h"
#include "tensor/basis.h"
#include "tensor/tensorproduct.h"

using namespace Tensor;

typedef IrredTensor<std::complex<double> > SimpleOperator;

typedef ReducibleTensor<std::complex<double>, BasisList, BasisList> SimpleRedOperator;

typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
                            VectorBasis, 
                            VectorBasis> MatrixOperator;

template <typename T>
struct BasicStateComponent;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, BasicStateComponent<T> const& Op);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, BasicStateComponent<T>& Op);

template <typename T>
struct BasicStateComponent
{
   typedef T value_type;
   typedef T OperatorType;
   typedef T operator_type;

   typedef LinearAlgebra::Vector<value_type> DataType;

   typedef typename DataType::const_iterator const_iterator;
   typedef typename DataType::iterator       iterator;

   typedef VectorBasis BasisType;

   BasicStateComponent() {}

   BasicStateComponent(BasisList const& SBasis_, VectorBasis const& V1, 
		    VectorBasis const& V2);

   QuantumNumbers::SymmetryList GetSymmetryList() const { return SBasis.GetSymmetryList(); }

   // returns true if this is a zero matrix
   bool is_null() const { return SBasis.size() == 0; }

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

   static BasicStateComponent<T> ConstructFullBasis1(BasisList const& LocalBasis, 
                                                       VectorBasis const& Basis2);
   static BasicStateComponent<T> ConstructFullBasis2(VectorBasis const& Basis1, 
                                                       BasisList const& LocalBasis);

   // Verifies that Data[i].TransformsAs() == SBasis[i]
   void check_structure() const;

   void debug_check_structure() const;
   
   BasisList SBasis;
   VectorBasis VBasis1, VBasis2;
   DataType Data;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, BasicStateComponent const& Op);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, BasicStateComponent& Op);

};

template <typename T>
void BasicStateComponent<T>::delta_shift(QuantumNumbers::QuantumNumber const& q,
					   QuantumNumbers::Projection const& Delta)
{
   DEBUG_TRACE("before")(scalar_prod(*this, herm(*this)))(scalar_prod(herm(*this), *this));
   VBasis1 = delta_shift(VBasis1, Delta);
   VBasis2 = delta_shift(VBasis2, Delta);
   for (typename BasicStateComponent<T>::iterator I = this->begin(); I != this->end(); ++I)
   {
      DEBUG_TRACE(*I)(I->Basis1())(I->Basis2());
      *I = Tensor::delta_shift(*I, q, Delta, VBasis1, VBasis2);
      DEBUG_TRACE("after component")(scalar_prod(*I, herm(*I)))(scalar_prod(herm(*I), *I))(I->TransformsAs());
      DEBUG_TRACE(*I)(I->Basis1())(I->Basis2());
   }
   DEBUG_TRACE("after")(scalar_prod(*this, herm(*this)))(scalar_prod(herm(*this), *this));
}

template <typename T>
void
BasicStateComponent<T>::check_structure() const
{
   CHECK_EQUAL(Data.size(), SBasis.size())("StateComponent local basis size mismatch");
   for (unsigned i = 0; i < Data.size(); ++i)
   {
      CHECK_EQUAL(Data[i].TransformsAs(), SBasis[i])("StateComponent TransformsAs() doesn't match local basis");
   }
}

template <typename T>
void BasicStateComponent<T>::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

typedef BasicStateComponent<MatrixOperator> StateComponent;

inline
StateComponent& operator*=(StateComponent& x, double y)
{
   for (unsigned i = 0; i < x.size(); ++i)
      x[i] *= y;
   return x;
}

inline
StateComponent& operator*=(StateComponent& x, std::complex<double> y)
{
   for (unsigned i = 0; i < x.size(); ++i)
      x[i] *= y;
   return x;
}

std::ostream& operator<<(std::ostream& out, StateComponent const& Psi);

StateComponent make_vacuum_state(QuantumNumbers::SymmetryList const& S);

StateComponent make_vacuum_state(QuantumNumbers::QuantumNumber const& Q);

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

StateComponent
operator+(StateComponent const& x, StateComponent const& y);

StateComponent
operator-(StateComponent const& x, StateComponent const& y);

StateComponent&
operator+=(StateComponent& x, StateComponent const& y);

StateComponent&
operator-=(StateComponent& x, StateComponent const& y);

namespace LinearAlgebra
{

template <>
struct interface<StateComponent>
{
   typedef void type;
};

// hermitian conjugation

template <>
struct Herm<StateComponent>
{
   typedef HermitianProxy<StateComponent> result_type;
   typedef StateComponent const& argument_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x);
   }
};

// complex conjugation

inline
StateComponent conj(StateComponent const& x)
{
   StateComponent Result(x);
   for (StateComponent::iterator I = Result.begin(); I != Result.end(); ++I)
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
struct ScalarProd<StateComponent, HermitianProxy<StateComponent>, Func>
{
   typedef typename Func::result_type result_type;
   typedef StateComponent const& first_argument_type;
   typedef HermitianProxy<StateComponent> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
   result_type operator()(first_argument_type x, second_argument_type y, Func f) const;
};

template <typename Func>
struct ScalarProd<HermitianProxy<StateComponent>, StateComponent, Func>
{
   typedef typename Func::result_type result_type;
   typedef HermitianProxy<StateComponent> const& first_argument_type;
   typedef StateComponent const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
   result_type operator()(first_argument_type x, second_argument_type y, Func f) const;
};

// the common case is explicitly specialized
template <>
struct ScalarProd<StateComponent, HermitianProxy<StateComponent> >
{
   typedef MatrixOperator result_type;
   typedef StateComponent const& first_argument_type;
   typedef HermitianProxy<StateComponent> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

template <>
struct ScalarProd<HermitianProxy<StateComponent>, StateComponent>
{
   typedef MatrixOperator result_type;
   typedef HermitianProxy<StateComponent> const& first_argument_type;
   typedef StateComponent const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

} // namespace LinearAlgebra


std::complex<double>
inner_prod(StateComponent const& x, StateComponent const& y);

inline
std::complex<double>
inner_prod(StateComponent const& x, StateComponent const& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   std::complex<double> r = 0.0;
   for (unsigned i = 0; i < x.size(); ++i)
   {
      r += inner_prod(x[i], y[i]);
   }
   return r;
}

inline
double
norm_frob_sq(StateComponent const& x)
{
   double r = 0;
   for (StateComponent::const_iterator cmp = x.begin(); cmp != x.end(); ++cmp)
   {
      r += norm_frob_sq(*cmp);
   }
   return r;
}

inline
double
norm_frob(StateComponent const& x)
{
   return std::sqrt(norm_frob_sq(x));
}

inline
StateComponent operator*(StateComponent const& x, double y)
{
   StateComponent Res(x);
   Res *= y;
   return Res;
}

inline
StateComponent operator*(StateComponent const& x, std::complex<double> y)
{
   StateComponent Res(x);
   Res *= y;
   return Res;
}

inline
StateComponent operator*(double y, StateComponent const& x)
{
   StateComponent Res(x);
   Res *= y;
   return Res;
}

inline
StateComponent operator*(std::complex<double> y, StateComponent const& x)
{
   StateComponent Res(x);
   Res *= y;
   return Res;
}

// A reflection of an MPS, which is a conjugate-transpose combined with a flip-conjugation
// of the quantum numbers.  This preserves the local basis but transforms the
// auxiliary basis into the adjoint.
StateComponent reflect(StateComponent const& S);

// does Result' = sum_{s,t} M(s,t) * A^s * herm(B^t)
// Result' transforms the same way as M.
MatrixOperator operator_prod(SimpleOperator const& M, 
                             StateComponent const& A, 
                             LinearAlgebra::HermitianProxy<StateComponent> const& B);

MatrixOperator operator_prod(SimpleRedOperator const& M, 
                             StateComponent const& A, 
                             LinearAlgebra::HermitianProxy<StateComponent> const& B);

// does Result' = sum_{s,t} conj(M(s,t)) * herm(A^s) * B^t
// Result' transforms as adjoint(M).
MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M, 
                             LinearAlgebra::HermitianProxy<StateComponent> const& A, 
                             StateComponent const& B);

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleRedOperator> const& M, 
                             LinearAlgebra::HermitianProxy<StateComponent> const& A, 
                             StateComponent const& B);

// does Result' = sum_{s,t} M(t,s) * A^t * E * herm(B^s)
MatrixOperator operator_prod(SimpleOperator const& M, 
                             StateComponent const& A, 
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<StateComponent> const& B,
                             QuantumNumbers::QuantumNumber const& q);

MatrixOperator operator_prod(SimpleOperator const& M, 
                             StateComponent const& A, 
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<StateComponent> const& B);

MatrixOperator operator_prod(SimpleRedOperator const& M, 
                             StateComponent const& A, 
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<StateComponent> const& B,
                             QuantumNumbers::QuantumNumber const& q);

// does Result' = sum_{s,t} M(t,s) * herm(A^t) * E * B^s
MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M, 
                             LinearAlgebra::HermitianProxy<StateComponent> const& A, 
                             MatrixOperator const& E,
                             StateComponent const& B,
                             QuantumNumbers::QuantumNumber const& q);

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleRedOperator> const& M, 
                             LinearAlgebra::HermitianProxy<StateComponent> const& A, 
                             MatrixOperator const& E,
                             StateComponent const& B,
                             QuantumNumbers::QuantumNumber const& q);

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M, 
                             LinearAlgebra::HermitianProxy<StateComponent> const& A, 
                             MatrixOperator const& E,
                             StateComponent const& B);

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<StateComponent> const& A, 
                             MatrixOperator const& E,
                             StateComponent const& B);

MatrixOperator operator_prod(StateComponent const& A, 
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<StateComponent> const& B);

// Variant for A,B matrices arising from a triangular MPO, where 
// A.front() = identity and B.back() = identity
MatrixOperator
operator_prod_regular(StateComponent const& A, 
                      MatrixOperator const& E,
                      LinearAlgebra::HermitianProxy<StateComponent> const& B);


// returns the operator M(s,t) = trace(herm(A^s) B^t)
// This is necessarily a scalar operator (the trace of a non-scalar operator is zero).
SimpleOperator trace_prod(LinearAlgebra::HermitianProxy<StateComponent> const& A,
                          StateComponent const& B);

// returns the operator M(s,t) = trace(A^s herm(B^t))
// This is necessarily a scalar operator (the trace of a non-scalar operator is zero).
SimpleOperator trace_prod(StateComponent const& A,
                          LinearAlgebra::HermitianProxy<StateComponent> const& B);

// product of a StateComponent and a MatrixOperator.
// The BasicOperator must transform as a scalar - eg typically it will be unitary.
StateComponent prod(StateComponent const& A, MatrixOperator const& Op);

StateComponent prod(StateComponent const& A, HermitianProxy<MatrixOperator> const& Op);

// product of a BasicOperator and a StateComponent.
// The BasicOperator must transform as a scalar - eg typically it will be unitary.
StateComponent prod(MatrixOperator const& Op, StateComponent const& A);

StateComponent prod(HermitianProxy<MatrixOperator> const& Op, StateComponent const& A);

// returns the operator a'^s' = sum_s A^s * x(s,s').  x must transform as a scalar.
StateComponent local_prod(StateComponent const& A, SimpleOperator const& x);

// returns the operator a'^s' = sum_s x(s',s) * A^s.  x must transform as a scalar.
StateComponent local_prod(SimpleOperator const& x, StateComponent const& A);

// constructs the tensor product in the local basis, the matrices
// C[k1,k2] = A[k1] * B[k2]
// This is a coarse-graining operation
StateComponent local_tensor_prod(StateComponent const& A, StateComponent const& B);

// Calculates the 'reyni product' A * B
// This is the component Result'[s] = A[s] \otimes B[s]
// corresponds to the wavefunction where the amplitude of some state
// s1, s2, .... is A(s1,s2,...) * B(s1,s2,....)
StateComponent renyi_product(StateComponent const& A, StateComponent const& B);

StateComponent triple_prod(MatrixOperator const& Op1, 
                             StateComponent const& A, 
                             LinearAlgebra::HermitianProxy<MatrixOperator> const&Op2);

StateComponent triple_prod(LinearAlgebra::HermitianProxy<MatrixOperator> const& Op1, 
                             StateComponent const& A, 
                             MatrixOperator const&Op2);

// Constructs a StateComponent that represents the sum of A and B.
// The resulting state has Result'[s] = A[s] \oplus B[s]
StateComponent tensor_sum(StateComponent const& A, StateComponent const& B, 
                            SumBasis<VectorBasis> const& B1, 
                            SumBasis<VectorBasis> const& B2);

// Constructs a StateComponent that represents the sum of A and B,
// at the left boundary of the matrix product state.
// Precondition: A.Basis1() == B.Basis1()
// The resulting state has Result'[s] = (A[s], B[s])  (row-wise concatenation)
StateComponent tensor_row_sum(StateComponent const& A, 
                                StateComponent const& B, 
                                SumBasis<VectorBasis> const& B2);

// Constructs a StateComponent that represents the sum of A and B,
// at the right boundary of the matrix product state.
// Precondition: A.Basis2() == B.Basis2()
// The resulting state has Result'[s] = ( A[s] )
//                                      ( B[s] )  (column-wise concatenation)
StateComponent tensor_col_sum(StateComponent const& A, 
                                StateComponent const& B, 
                                SumBasis<VectorBasis> const& B1);

// Returns the diagonal components of the operator F given by
// F(x) = operator_prod(A, x, herm(B))
MatrixOperator extract_diagonal(StateComponent const& A, 
                                LinearAlgebra::HermitianProxy<StateComponent> const& B);

// shift the basis of an StateComponent by some quantum number; the local basis
// is shifted by QL, the incoming matrix basis is shifted by QM, the outgoing matrix basis
// is shifted by QM+QL.  This only works if QL and QM are degree 1 reps.
StateComponent ShiftLocalBasis(StateComponent const& Op, QuantumNumber QL, QuantumNumber QM);

StateComponent delta_shift(StateComponent const& Op, QuantumNumber const& q);

// scales the quantum numbers of all bases of an StateComponent by some factor.  
// This only makes sense for U(1) quantum numbers.  Name is the name of the quantum number
// to scale.
StateComponent ScaleBasisU1(StateComponent const& Op, 
			      std::string const& Name, 
			      double Factor);

// returns a StateComponent that is identical except for the modified SymmetryList,
// which must have the same symmetry types in the same order, but is allowed to be
// renamed.
StateComponent RenameSymmetry(StateComponent const& Op, SymmetryList const& NewSL);

// re-orders the basis
StateComponent ReorderLocalBasis(StateComponent const& Op, std::list<int> const& NewOrder);

// returns a StateComponent with a symmetry list that is compatible with the old one
// in the sense of all old symmetries are in the new list (possibly in a different order).
// new symmetries are allowed, and get quantum number 0.
StateComponent CoerceSymmetryList(StateComponent const& Op, SymmetryList const& NewSL);

// expand the basis

enum Normalization { Intensive, Extensive };

MatrixOperator ExpandBasis1(StateComponent& A, Normalization n = Intensive);
MatrixOperator ExpandBasis2(StateComponent& A, Normalization n = Intensive);

// Expand the basis, but incorporate only those matrix elements that are actually used in
// the A-matrix.  Zero matrix elements are excluded.
MatrixOperator ExpandBasis1Used(StateComponent& A, std::vector<int> const& Used);
MatrixOperator ExpandBasis2Used(StateComponent& A, std::vector<int> const& Used);

// Constructs an 'expanded' basis given the input right hand basis
StateComponent ConstructFromRightBasis(BasisList const& LocalBasis,
					 VectorBasis const& Right);

// Constructs an 'expanded' basis given the input left hand basis
StateComponent ConstructFromLeftBasis(BasisList const& LocalBasis,
					VectorBasis const& LeftBasis);

// StateProject
// A functor to project onto the given component of a StateComponent
struct StateProject
{
   typedef StateComponent::OperatorType result_type;
   typedef StateComponent argument_type;
   StateProject(int s_) : s(s_) {}

   result_type operator()(StateComponent const& x) const { return x[s]; }

   int s;
};

template <typename FwdIter>
StateComponent tensor_col_accumulate(FwdIter first, FwdIter last, 
                                       SumBasis<VectorBasis> const& B1)
{
   //   PRECONDITION_EQUAL(A.LocalBasis(), B.LocalBasis());
   //   PRECONDITION_EQUAL(A.Basis2(), B.Basis2());
   typedef boost::transform_iterator<StateProject, FwdIter> StateIter;
   StateComponent Result(first->LocalBasis(), B1, first->Basis2());
   for (int s = 0; s < first->LocalBasis().size(); ++s)
   {
      Result[s] = tensor_col_accumulate(StateIter(first, StateProject(s)), 
					StateIter(last, StateProject(s)), B1);
   }
   return Result;
}

template <typename FwdIter>
StateComponent tensor_row_accumulate(FwdIter first, FwdIter last,
                                       SumBasis<VectorBasis> const& B2)
{
   //   PRECONDITION_EQUAL(A.LocalBasis(), B.LocalBasis());
   //   PRECONDITION_EQUAL(A.Basis1(), B.Basis1());
   typedef boost::transform_iterator<StateProject, FwdIter> StateIter;
   StateComponent Result(first->LocalBasis(), first->Basis1(), B2);
   for (int s = 0; s < first->size(); ++s)
   {
      Result[s] = tensor_row_accumulate(StateIter(first, StateProject(s)), 
					StateIter(last, StateProject(s)), B2);
   }
   return Result;
}

template <typename FwdIter>
StateComponent tensor_accumulate(FwdIter first, FwdIter last,
                                   SumBasis<VectorBasis> const& B1, 
                                   SumBasis<VectorBasis> const& B2)
{
   //   PRECONDITION_EQUAL(A.LocalBasis(), B.LocalBasis());
   typedef boost::transform_iterator<StateProject, FwdIter> StateIter;
   StateComponent Result(first->LocalBasis(), B1, B2);
   for (int s = 0; s < first->LocalBasis().size(); ++s)
   {
      Result[s] = tensor_accumulate(StateIter(first, StateProject(s)), 
				    StateIter(last, StateProject(s)), B1, B2);
   }
   return Result;
}

// utility functions to generate random matrices.
// Typically used to initialize iterative eigensolvers etc
MatrixOperator MakeRandomMatrixOperator(VectorBasis const& B1, VectorBasis const& B2,
					QuantumNumber q);

inline
MatrixOperator MakeRandomMatrixOperator(VectorBasis const& B1, VectorBasis const& B2)
{
   return MakeRandomMatrixOperator(B1, B2, QuantumNumber(B1.GetSymmetryList()));
}

StateComponent
MakeRandomStateComponent(BasisList const& Local, VectorBasis const& B1, VectorBasis const& B2);


// this belongs somewhere else
MatrixOperator RenameSymmetry(MatrixOperator const& Op, SymmetryList const& NewSL);

#include "state_component.cc"

#endif
