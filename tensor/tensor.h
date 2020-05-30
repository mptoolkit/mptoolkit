// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensor.h
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

/*
  Tensor product library main file

  Created 2004 Ian McCulloch

*/

#if !defined(TENSOR_H_WUIYR74839Y78HV8975HY49)
#define TENSOR_H_WUIYR74839Y78HV8975HY49

#include "basis.h"
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/coefficient_trace.h"
#include "linearalgebra/coefficient_operations.h"
#include "quantumnumbers/quantumnumber.h"
#include "pstream/pstream.h"
#include <boost/mpl/assert.hpp>
#include <cmath>

// first, we declare some new operations that extend the LinearAlgebra operations
namespace LinearAlgebra
{

// A proxy class to represent the Hermitian conjugate of a tensor.
template <typename T>
struct HermitianProxy
{
   typedef T base_type;
   typedef typename T::value_type value_type;
   typedef T const& reference;

   explicit HermitianProxy(T const& x) : x_(x) {}

   reference base() const { return x_; }

   private:
      reference x_;
};

// A proxy class to represent the tensor conjugate.
template <typename T>
struct TensorConjugateProxy
{
   typedef T base_type;
   typedef typename T::value_type value_type;
   typedef T const& reference;

   explicit TensorConjugateProxy(T const& x) : x_(x) {}

   reference base() const { return x_; }

   private:
      reference x_;
};

// ScalarProd
// Calculates the dot product of two irreducible tensors.
template <typename S, typename T, typename Functor = Multiplication<typename S::value_type,
                                                                    typename T::value_type> >
struct ScalarProd {};

template <typename S, typename T>
typename ScalarProd<S, T>::result_type
scalar_prod(S const& x, T const& y)
{
   return ScalarProd<S,T>()(x,y);
}

template <typename S, typename T, typename Func>
typename ScalarProd<S, T, Func>::result_type
scalar_prod(S const& x, T const& y, Func const& f)
{
   return ScalarProd<S,T, Func>()(x,y,f);
}

// Adjoint
// By default, this is the same as Herm
// For IrredTensor, it is the 'Adjoint' operation
template <typename T>
struct Adjoint : public LinearAlgebra::Herm<T> { };

// Adjoint applied to a tensor adds a function to apply to the value_type
template <typename T, typename F = Adjoint<typename T::value_type> >
struct TensorAdjoint;

template <typename T>
inline
typename Adjoint<T>::result_type
adjoint(T const& x)
{
   return Adjoint<T>()(x);
}

template <typename T, typename Func>
inline
typename TensorAdjoint<T, Func>::result_type
adjoint(T const& x, Func const& f)
{
   return TensorAdjoint<T, Func>(f)(x);
}

// InvAdjoint
// By default, the 'inverse adjoint' is also the same as Herm
// For IrredTensor, it is the inverse of the 'Adjoint' operation
template <typename T>
struct InvAdjoint : public LinearAlgebra::Herm<T> { };

template <typename T>
inline
typename InvAdjoint<T>::result_type
inv_adjoint(T const& x)
{
   return InvAdjoint<T>()(x);
}

template <typename T, typename F = InvAdjoint<typename T::value_type> >
struct TensorInvAdjoint;

template <typename T, typename Func>
inline
typename TensorInvAdjoint<T, Func>::result_type
inv_adjoint(T const& x, Func const& f)
{
   return TensorInvAdjoint<T, Func>()(x);
}

//
// FlipConjugate
//

// for a non-tensor, the flip conjugate is just the ordinary conjugate
template <typename T>
struct FlipConjugate : public LinearAlgebra::Conj<T> { };

template <typename T, typename F = FlipConjugate<typename T::value_type> >
struct TensorFlipConjugate;

template <typename T>
inline
typename FlipConjugate<T>::result_type
flip_conj(T const& x)
{
   return FlipConjugate<T>()(x);
}

template <typename T, typename Func>
inline
typename TensorFlipConjugate<T, Func>::result_type
flip_conj(T const& x, Func const& f)
{
   return TensorFlipConjugate<T, Func>()(x);
}

// IrredProd
// Extends Multiplication to irreducible tensors.
template <typename Op1Type, typename Op2Type,
          typename Nest = Multiplication<typename Op1Type::value_type,
                                         typename Op2Type::value_type> >
struct IrredProd {};

template <typename S, typename T>
inline
typename IrredProd<S, T>::result_type
prod(S const& x, T const& y, QuantumNumbers::QuantumNumber const& q)
{
   return IrredProd<S, T>(q)(x,y);
}

template <typename S, typename T, typename Nest>
inline
typename IrredProd<S, T, Nest>::result_type
prod(S const& x, T const& y, QuantumNumbers::QuantumNumber const& q, Nest n)
{
   return IrredProd<S, T, Nest>(q)(x,y,n);
}

// dot() is a shortcut for prod(x,y,ident)
template <typename S, typename T>
inline
typename IrredProd<S, T>::result_type
dot(S const& x, T const& y)
{
   return IrredProd<S, T>(QuantumNumbers::QuantumNumber(x.GetSymmetryList()))(x,y);
}

// outer_product: outer product of tensor operators.
// We choose among the possible transform_targets the
// quantum number with the largest degree.
template <typename S, typename T>
inline
typename IrredProd<S, T>::result_type
outer(S const& x, T const& y)
{
   QuantumNumbers::QuantumNumberList L = transform_targets(x.TransformsAs(), y.TransformsAs());
   QuantumNumbers::QuantumNumber q = L[0];
   bool Unique = true;
   for (unsigned i = 1; i < L.size(); ++i)
   {
      if (degree(L[i]) > degree(q))
      {
         q = L[i];
         Unique = true;
      }
      else if (degree(L[i]) == degree(q))
      {
         Unique = false;
      }
   }
   CHECK(Unique)("outer product is not defined for these operators")
      (x.TransformsAs())(y.TransformsAs());

   // The factor here is deduced for spin-1 to spin-2 operators,
   // as the correct factor to get dot(outer(S,S), outer(S,S)) correct.  The test is that
   // S^4, when evaluated on a scalar wavefunction, should produce the same result
   // whether we use symmetries or not, and this is easy to test with no symmetries (S = Sx+Sy+Sz).
   // The non-zero components with SU(2) are
   // dot(S,S)*dot(S,S) and dot(outer(S,S), outer(S,S)).  The coupling constant below
   // was determined by matching coefficients.
   return std::sqrt(double(degree(x.TransformsAs()) + degree(y.TransformsAs())) / degree(q)) * prod(x,y,q);
}

template <typename S, typename T>
inline
typename IrredProd<S, T>::result_type
cross(S const& x, T const& y)
{
   CHECK(cross_product_exists(x.TransformsAs(), y.TransformsAs()))
      ("Cross product does not exist for these operators")
      (x.TransformsAs())(y.TransformsAs());

   return cross_product_factor(x.TransformsAs(), y.TransformsAs())
      * prod(x, y, cross_product_transforms_as(x.TransformsAs(), y.TransformsAs()));
}

template <typename S>
inline
typename IrredProd<S, S>::result_type
pow(S const& x, int n)
{
   if (n == 0)
   {
      CHECK_EQUAL(x.Basis1(), x.Basis2());
      return IrredProd<S, S>::result_type::make_identity(x.Basis1());
   }
   else if (n%2 == 0)
   {
      return pow(x*x, n/2);
   }
   else if (n == 1)
   {
      return x;
   }
   else
   {
      return x*pow(x*x, (n-1)/2);
   }

}

} // namespace LinearAlgebra

namespace Tensor
{

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::Projection;
using QuantumNumbers::SymmetryList;

using LinearAlgebra::size_type;
using LinearAlgebra::iterate;
using LinearAlgebra::iterator;
using LinearAlgebra::const_iterator;
using LinearAlgebra::inner_iterator;
using LinearAlgebra::const_inner_iterator;
using LinearAlgebra::set_element;

using LinearAlgebra::norm_frob;
using LinearAlgebra::norm_frob_sq;
using LinearAlgebra::inner_prod;

using LinearAlgebra::equal;
using LinearAlgebra::operator+;
using LinearAlgebra::operator-;
using LinearAlgebra::operator*;
using LinearAlgebra::transform;
using LinearAlgebra::trace;
using LinearAlgebra::prod;
using LinearAlgebra::inner_prod;
using LinearAlgebra::parallel_prod;
using LinearAlgebra::conj;
using LinearAlgebra::herm;
using LinearAlgebra::flip_conj;
using LinearAlgebra::scalar_prod;
using LinearAlgebra::adjoint;
using LinearAlgebra::inv_adjoint;
using LinearAlgebra::nnz;
using LinearAlgebra::exp;

// Structure traits types
struct DefaultStructure
{
   template <typename T>
   using value = LinearAlgebra::SparseMatrix<T>;
};

struct DiagonalStructure
{
   template <typename T>
   using value = LinearAlgebra::DiagonalMatrix<T>;
};

// For operations that act on the structure and potentially modify it,
// we need to map the transformed concrete type back to a structure trait.
// In principle this allows us to generalize operations to mutable structures,
// although currently none are implemented.

template <typename T>
struct StructureOf {};

template <typename T>
struct StructureOf<LinearAlgebra::SparseMatrix<T>>
{
   using type = DefaultStructure;
};

template <typename T>
struct StructureOf<LinearAlgebra::DiagonalMatrix<T>>
{
   using type = DiagonalStructure;
};

template <typename T,
          typename Basis1T = BasisList,
          typename Basis2T = Basis1T,
          typename Structure = DefaultStructure >
          //          typename Structure = LinearAlgebra::SparseMatrix<T> >
class IrredTensor;

template <typename T, typename B1, typename B2, typename S>
PStream::opstream& operator<<(PStream::opstream& out, IrredTensor<T, B1, B2, S> const& x);

template <typename T, typename B1, typename B2, typename S>
PStream::ipstream& operator>>(PStream::ipstream& in, IrredTensor<T, B1, B2, S>& x);

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
CoerceSymmetryList(IrredTensor<T, B1, B2, S> const& t, SymmetryList const& sl);
//   __attribute__((warn_unused_result));

template <typename T, typename B1, typename B2, typename S>
void
CoerceSymmetryListInPlace(IrredTensor<T, B1, B2, S>& t, SymmetryList const& sl);

template <typename T, typename Basis1T, typename Basis2T, typename Structure>
class IrredTensor
{
   public:
      typedef T                                                     value_type;
      //      typedef Structure                                             MatrixType;

      using StructureType = Structure;
      using MatrixType = typename Structure::template value<T>;

      BOOST_MPL_ASSERT((boost::is_same<value_type,
                        typename LinearAlgebra::interface<MatrixType>::value_type>));

      typedef typename iterator<MatrixType>::type iterator;
      typedef typename inner_iterator<MatrixType>::type inner_iterator;
      typedef typename const_iterator<MatrixType>::type const_iterator;
      typedef typename const_inner_iterator<MatrixType>::type const_inner_iterator;

      typedef Basis1T basis1_type;
      typedef Basis2T basis2_type;

      IrredTensor() = default;

      IrredTensor(IrredTensor const& Other) = default;
      IrredTensor(IrredTensor&& Other) = default;

      IrredTensor(basis1_type const& Basis, QuantumNumber const& Trans);

      IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2,
                  QuantumNumber const& Trans);

      IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2,
                  QuantumNumber const& Trans, MatrixType const& Data);

      // In this variant, the quantum number defaults to the identity
      IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2);

      explicit IrredTensor(basis1_type const& Basis);

      template <typename U, typename US>
      IrredTensor(IrredTensor<U, Basis1T, Basis2T, US> const& r)
         : Basis1_(r.Basis1()), Basis2_(r.Basis2()), Trans_(r.TransformsAs()),
           Data_(r.data()) {}

      IrredTensor& operator=(IrredTensor const& Other) = default;
      IrredTensor& operator=(IrredTensor&& Other) = default;

      template <typename U, typename US>
      IrredTensor& operator=(IrredTensor<U, Basis1T, Basis2T, US> const& r)
      {
         Basis1_ = r.Basis1(); Basis2_ = r.Basis2(); Trans_ = r.TransformsAs();
         Data_ = r.data();
         return *this;
      }

      bool is_null() const { return Trans_.is_null(); }

      IrredTensor& operator+=(IrredTensor const& Op);
      IrredTensor& operator-=(IrredTensor const& Op);

      template <typename U>
      IrredTensor& operator*=(U const& x) { Data_ *= x; return *this; }

      basis1_type const& Basis1() const { return Basis1_; }
      basis2_type const& Basis2() const { return Basis2_; }

      SymmetryList GetSymmetryList() const { return Basis1_.GetSymmetryList(); }

      size_type size1() const { return Basis1_.size(); }
      size_type size2() const { return Basis2_.size(); }

      typename LinearAlgebra::MatrixBracket<MatrixType&, size_type, size_type>::result_type
      operator()(size_type i, size_type j);

      typename LinearAlgebra::MatrixBracket<MatrixType, size_type, size_type>::result_type
      operator()(size_type i, size_type j) const;

      MatrixType& data() { return Data_; }
      MatrixType const& data() const { return Data_; }

      QuantumNumber const& TransformsAs() const { return Trans_; }

      QuantumNumber const& qn1(size_type i) const { return Basis1_[i]; }
      QuantumNumber const& qn2(size_type j) const { return Basis2_[j]; }

      // in-place delta_shfit
      void delta_shift(QuantumNumber const& q);

      // checks to see if there are any non-zero (structural) matrix elements
      // that are not allowed by selection rules.  If so, an assertion is
      // triggered.
      void check_structure() const;

      void debug_check_structure() const;

      static IrredTensor<T, Basis1T, Basis2T, Structure>
      make_identity(basis1_type const& b);

      // force change the symmetry list
      void CoerceSymmetryList(SymmetryList const& sl);

   private:
      basis1_type Basis1_;
      basis2_type Basis2_;
      QuantumNumber Trans_;
      MatrixType Data_;

      friend PStream::opstream& operator<< <>(PStream::opstream& out, IrredTensor const& x);
      friend PStream::ipstream& operator>> <>(PStream::ipstream& in, IrredTensor& x);

      template <typename T2, typename B12, typename B22, typename S2>
      friend IrredTensor<T2, B12, B22, S2>
      CoerceSymmetryList(IrredTensor<T2, B12, B22, S2> const& t, SymmetryList const& sl);


      //      friend IrredTensor ::Tensor::CoerceSymmetryList<>(IrredTensor const& t, SymmetryList const& sl);
};

// text I/O

template <typename T, typename B1, typename B2, typename S>
std::ostream& operator<<(std::ostream& out, IrredTensor<T, B1, B2, S> const& Op);

template <typename T, typename B1, typename B2, typename S>
std::string show_projections(IrredTensor<T, B1, B2, S> const& Op);

// construction

template <typename T, typename B, typename B2, typename S>
IrredTensor<T, B, B2, S>
IrredTensor<T, B, B2, S>::make_identity(B const& b)
{
   IrredTensor<T, B, B, S> Result(b, b, QuantumNumber(b.GetSymmetryList()));
   for (std::size_t i = 0; i < b.size(); ++i)
   {
      Result(i,i) = ::Tensor::make_identity<T>(b, i);
   }
   return Result;
}

// MakeIdentityFrom : construct an identity tensor over the same basis as a given tensor
template <typename T, typename Basis1T, typename Basis2T, typename Structure>
IrredTensor<T, Basis1T, Basis2T, Structure>
MakeIdentityFrom(IrredTensor<T, Basis1T, Basis2T, Structure> const& s)
{
   CHECK_EQUAL(s.Basis1(), s.Basis2());
   return IrredTensor<T, Basis1T, Basis2T, Structure>::make_identity(s.Basis1());
}

} // namespace Tensor

namespace LinearAlgebra
{

template <typename T, typename B1, typename B2, typename S>
struct interface<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef void type;
};

// iterators

template <typename T, typename B1, typename B2, typename S>
struct Iterate<Tensor::IrredTensor<T, B1, B2, S>&>
{
   typedef Tensor::IrredTensor<T, B1, B2, S>& argument_type;
   typedef typename Tensor::IrredTensor<T, B1, B2, S>::iterator result_type;

   result_type operator()(argument_type x) const
   {
      return iterate(x.data());
   }
};

template <typename T, typename B1, typename B2, typename S>
struct Iterate<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef typename Tensor::IrredTensor<T, B1, B2, S>::const_iterator result_type;

   result_type operator()(argument_type x) const
   {
      return iterate(x.data());
   }
};

// Assign

template <typename T1, typename B1, typename B2, typename S1, typename T2, typename S2>
struct Assign<Tensor::IrredTensor<T1, B1, B2, S1>&, Tensor::IrredTensor<T2, B1, B2, S2>>
{
   using result_type = void;
   using first_argument_type = Tensor::IrredTensor<T1, B1, B2, S1>&;
   using second_argument_type = Tensor::IrredTensor<T2, B1, B2, S2>;
   void operator()(Tensor::IrredTensor<T1, B1, B2, S1>& x,
                   Tensor::IrredTensor<T2, B1, B2, S2> const& y)
   {
      x = y;
   }
};

// zero

template <typename T, typename B1, typename B2, typename S>
struct Zero<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> result_type;
   result_type operator()() const { return result_type(); }
};

template <typename T, typename B1, typename B2, typename S>
struct IsZero<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef bool result_type;
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return x.is_null();
   }
};

template <typename T, typename B1, typename B2, typename S>
struct ZeroAll<Tensor::IrredTensor<T, B1, B2, S>&>
{
   typedef void result_type;
   typedef Tensor::IrredTensor<T, B1, B2, S>& argument_type;
   void operator()(argument_type x) const
   {
      zero_all(x.data());
   }
};

// nnz

template <typename T, typename B1, typename B2, typename S>
struct NNZ<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef std::size_t result_type;
   result_type operator()(argument_type x) const { return nnz(x.data()); }
};

// comparison

template <typename T, typename B1, typename B2, typename S1, typename S2>
struct Equal<Tensor::IrredTensor<T, B1, B2, S1>,  Tensor::IrredTensor<T, B1, B2, S2> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S1> const& first_argument_type;
   typedef Tensor::IrredTensor<T, B1, B2, S2> const& second_argument_type;
   typedef bool result_type;

   Equal(double tol = LinearAlgebra::default_tolerance()) : tol_(tol) {}

   bool operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null() && y.is_null()) return true;

      DEBUG_PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
      DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
      return x.TransformsAs() == y.TransformsAs()
         && LinearAlgebra::equal(x.data(), y.data(), tol_);
   }

   double tol_;
};

// element access

template <typename T, typename B1, typename B2, typename S, typename V>
struct SetMatrixElement<Tensor::IrredTensor<T, B1, B2, S>&, V>
{
   typedef void result_type;
   typedef Tensor::IrredTensor<T, B1, B2, S>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef V const& fourth_argument_type;
   void operator()(first_argument_type m, size_type i, size_type j, V const& x) const
   {
      DEBUG_PRECONDITION(is_transform_target(m.qn2(j), m.TransformsAs(), m.qn1(i)));
      set_element(m.data(), i, j, x);
   }
};

// arithmetic

template <typename T, typename B1, typename B2, typename S>
struct Negate<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef Tensor::IrredTensor<T, B1, B2, S> result_type;

   result_type operator()(Tensor::IrredTensor<T, B1, B2, S> const& x) const
   {
      return Tensor::IrredTensor<T, B1, B2, S>(x.Basis1(), x.Basis2(), x.TransformsAs(), -x.data());
   }
};

template <typename T1, typename T2, typename B1, typename B2, typename S>
struct Addition<Tensor::IrredTensor<T1, B1, B2, S>, Tensor::IrredTensor<T2, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T1, B1, B2, S> const& first_argument_type;
   typedef Tensor::IrredTensor<T2, B1, B2, S> const& second_argument_type;
   typedef typename result_value<Addition<T1, T2> >::type ResultValue;
   typedef Tensor::IrredTensor<ResultValue, B1, B2, S> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null()) return y;
      if (y.is_null()) return x;
      PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
      PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
      return result_type(x.Basis1(), x.Basis2(), x.TransformsAs(), x.data() + y.data());
   }
};

template <typename T1, typename T2, typename B1, typename B2, typename S>
struct Subtraction<Tensor::IrredTensor<T1, B1, B2, S>, Tensor::IrredTensor<T2, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T1, B1, B2, S> const& first_argument_type;
   typedef Tensor::IrredTensor<T2, B1, B2, S> const& second_argument_type;
   typedef typename result_value<Subtraction<T1, T2> >::type ResultValue;
   typedef Tensor::IrredTensor<ResultValue, B1, B2, S> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null()) return -y;
      if (y.is_null()) return x;
      PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
      PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
      return result_type(x.Basis1(), x.Basis2(), x.TransformsAs(), x.data() - y.data());
   }
};

// multiply by scalar

template <typename T, typename B1, typename B2, typename S, typename U>
struct Multiplication<U, Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef U first_argument_type;
   typedef Tensor::IrredTensor<T, B1, B2, S> const& second_argument_type;
   typedef typename result_value<Multiplication<U, T> >::type ResultValue;
   typedef Tensor::IrredTensor<ResultValue, B1, B2, S> result_type;

   result_type operator()(first_argument_type y, second_argument_type x) const
   {
      return result_type(x.Basis1(), x.Basis2(), x.TransformsAs(), y * x.data());
   }
};

template <typename T, typename B1, typename B2, typename S, typename U>
struct Multiplication<Tensor::IrredTensor<T, B1, B2, S>, U>
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& first_argument_type;
   typedef U second_argument_type;
   typedef typename result_value<Multiplication<T, U> >::type ResultValue;
   typedef Tensor::IrredTensor<ResultValue, B1, B2, S> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(x.Basis1(), x.Basis2(), x.TransformsAs(), x.data() * y);
   }
};



template <typename T, typename B1, typename B2, typename S, typename Func>
struct Transform<Tensor::IrredTensor<T, B1, B2, S>, Func>
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& first_argument_type;
   typedef Func second_argument_type;
   typedef typename result_value<Transform<T, Func> >::type ResultValue;
   typedef Tensor::IrredTensor<ResultValue, B1, B2, S> result_type;

   result_type operator()(first_argument_type x, second_argument_type f) const
   {
      return result_type(x.Basis1(), x.Basis2(), x.TransformsAs(), transform(x.data(), f));
   }
};

// complex conjugation

template <typename T, typename B1, typename B2, typename S>
struct Conj<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef Tensor::IrredTensor<T, B1, B2, S> result_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x.Basis1(), x.Basis2(), x.TransformsAs(), conj(x.data()));
   }
};

// hermitian conjugation

template <typename T, typename B1, typename B2, typename S>
struct Herm<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef HermitianProxy<Tensor::IrredTensor<T, B1, B2, S> > result_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x);
   }
};

// tensor conjugation.  This is a new function
template <typename T>
struct TensorConj {};

template <typename T, typename B1, typename B2, typename S>
struct TensorConj<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef TensorConjugateProxy<Tensor::IrredTensor<T, B1, B2, S> > result_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x);
   }
};

template <typename T>
typename TensorConj<T>::result_type
tensor_conj(T const& x)
{
   return TensorConj<T>()(x);
}

// trace

struct traceCoefficient
{
   typedef double result_type;
   typedef double value_type;

   traceCoefficient(Tensor::BasisList const& B) : B_(B) {}

   double operator()(int x) const { return degree(B_[x]); }

   private:
      Tensor::BasisList const& B_;
};

template <typename T, typename B1, typename B2, typename S>
struct Trace<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef typename result_value<Trace<T> >::type result_type;

   result_type operator()(argument_type M) const
   {
      if (M.is_null())
         return result_type();
      PRECONDITION_EQUAL(M.Basis1(), M.Basis2());
      PRECONDITION(QuantumNumbers::is_scalar(M.TransformsAs()))(M.TransformsAs());

      return LinearAlgebra::coefficient_trace(M.data(), traceCoefficient(M.Basis1()));
   }
};

// norm_frob

template <typename T, typename B1, typename B2, typename S>
struct NormFrobSq<Tensor::IrredTensor<T, B1, B2, S> >
{
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;
   typedef typename result_value<NormFrobSq<T> >::type result_type;

   result_type operator()(argument_type M) const
   {
      result_type Result = 0;
      typename Tensor::IrredTensor<T, B1, B2, S>::const_iterator I = iterate(M);
      while (I)
      {
         // this version is faster, but assumes S is a row-major matrix
         //         Result += QuantumNumbers::degree(M.qn1(I.index())) * norm_frob_sq(*I);

         // this version is more general
         typename Tensor::IrredTensor<T, B1, B2, S>::const_inner_iterator J = iterate(I);
         while (J)
         {
            Result += QuantumNumbers::degree(M.qn1(J.index1())) * norm_frob_sq(*J);
            ++J;
         }

         ++I;
      }
      return Result;
   }
};

// inner_prod

template <typename T1, typename B1, typename B2, typename S1, typename T2, typename S2>
struct InnerProd<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2, B1, B2, S2> >
{
   typedef Tensor::IrredTensor<T1, B1, B2, S1> first_argument_type;
   typedef Tensor::IrredTensor<T2, B1, B2, S2> second_argument_type;
   typedef typename result_value<InnerProd<T1, T2> >::type result_type;

   result_type operator()(first_argument_type const& M1, second_argument_type const& M2) const
   {
      result_type Result = 0;
      if (M1.is_null() || M2.is_null())
         return Result;

      PRECONDITION_EQUAL(M1.TransformsAs(), M2.TransformsAs());
      PRECONDITION_EQUAL(M1.Basis1(), M2.Basis1());
      PRECONDITION_EQUAL(M1.Basis2(), M2.Basis2());

      // TODO: this assumes the Tensor::IrredTensor is implemented as a row-major matrix!
      // TODO: This fails for inner_prod for RealDiagonalOperator.
      typename Tensor::IrredTensor<T1, B1, B2, S1>::const_iterator I1 = iterate(M1);
      typename Tensor::IrredTensor<T2, B1, B2, S2>::const_iterator I2 = iterate(M2);
      while (I1)
      {
         DEBUG_CHECK(I2);
         Result += double(QuantumNumbers::degree(M1.qn1(I1.index()))) * inner_prod(*I1, *I2);
         ++I1;
         ++I2;
      }
      DEBUG_CHECK(!I2);
      return Result;
   }
};

// parallel_prod

template <typename T1, typename B1, typename B2, typename S1, typename T2, typename S2,
          typename Func>
struct ParallelProd<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2, B1, B2, S2>, Func>
{
   typedef Tensor::IrredTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef Tensor::IrredTensor<T2, B1, B2, S2> const& second_argument_type;
   typedef typename result_value<Func>::type result_type;

   ParallelProd(Func f) : f_(f) {}

   result_type operator()(first_argument_type M1, second_argument_type M2) const
   {
      if (M1.is_null() || M2.is_null()) return result_type();

      PRECONDITION_EQUAL(M1.TransformsAs(), M2.TransformsAs());
      PRECONDITION_EQUAL(M1.Basis1(), M2.Basis1());
      PRECONDITION_EQUAL(M1.Basis2(), M2.Basis2());

      // TODO: this assumes the Tensor::IrredTensor is implemented as a row-major matrix!
      typename Tensor::IrredTensor<T1, B1, B2, S1>::const_iterator I1 = iterate(M1);
      typename Tensor::IrredTensor<T2, B1, B2, S2>::const_iterator I2 = iterate(M2);

      if (!I1) return result_type();
      DEBUG_CHECK(I2);
      result_type Result = double(QuantumNumbers::degree(M1.qn1(I1.index()))) *
         (parallel_prod(*I1, *I2, f_));
      ++I1;
      ++I2;

      while (I1)
      {
         DEBUG_CHECK(I2);
         Result += double(QuantumNumbers::degree(M1.qn1(I1.index()))) *
            (parallel_prod(*I1, *I2, f_));
         ++I1;
         ++I2;
      }
      DEBUG_CHECK(!I2);
      return Result;
   }

   Func f_;
};

// scalar_prod - this is a new function

template <typename T1, typename B1, typename B2,
          typename T2, typename B3, typename Functor>
struct ScalarProd<Tensor::IrredTensor<T1, B1, B2, Tensor::DefaultStructure>,
                  HermitianProxy<Tensor::IrredTensor<T2, B3, B2, Tensor::DefaultStructure> >,
                  Functor>
{
   typedef Tensor::IrredTensor<T1, B1, B2, Tensor::DefaultStructure> const& first_argument_type;
   typedef HermitianProxy<Tensor::IrredTensor<T2, B3, B2, Tensor::DefaultStructure> > const& second_argument_type;
   typedef Tensor::IrredTensor<typename result_value<Functor>::type, B1, B3, Tensor::DefaultStructure> result_type;

   ScalarProd(Functor f = Functor()) : f_(f) {}

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null() || y.base().is_null()) return result_type();

      PRECONDITION_EQUAL(x.TransformsAs(), y.base().TransformsAs());
      DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.base().Basis2());

      typedef typename const_iterator<Tensor::IrredTensor<T1, B1, B2, Tensor::DefaultStructure> >::type iter1;
      typedef typename const_iterator<Tensor::IrredTensor<T2, B3, B2, Tensor::DefaultStructure> >::type iter2;

      result_type Result(x.Basis1(), y.base().Basis1(), QuantumNumbers::QuantumNumber(x.GetSymmetryList()));
      for (iter1 I1 = iterate(x); I1; ++I1)
      {
         for (iter2 I2 = iterate(y.base()); I2; ++I2)
         {
            // accept only scalar elements (diagonal in the quantum number)
            if (Result.Basis1()[I1.index()] == Result.Basis2()[I2.index()])
               add_element(Result.data(), I1.index(), I2.index(), parallel_prod(*I1, herm(*I2), f_));
         }
      }
      return Result;
   }

   Functor f_;
};

template <typename T1, typename B1, typename B2,
          typename T2, typename B3, typename Functor>
struct ScalarProd<Tensor::IrredTensor<T1, B1, B2, Tensor::DiagonalStructure>,
                  HermitianProxy<Tensor::IrredTensor<T2, B3, B2, Tensor::DiagonalStructure> >,
                  Functor>
{
   typedef Tensor::IrredTensor<T1, B1, B2, Tensor::DiagonalStructure> const& first_argument_type;
   typedef HermitianProxy<Tensor::IrredTensor<T2, B3, B2, Tensor::DiagonalStructure> > const& second_argument_type;
   typedef Tensor::IrredTensor<typename result_value<Functor>::type, B1, B3, Tensor::DiagonalStructure>
   result_type;

   ScalarProd(Functor f = Functor()) : f_(f) {}

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null() || y.base().is_null()) return result_type();

      PRECONDITION_EQUAL(x.TransformsAs(), y.base().TransformsAs());
      DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.base().Basis2());

      typedef typename const_iterator<Tensor::IrredTensor<T1, B1, B2, Tensor::DiagonalStructure> >::type iter1;
      typedef typename const_iterator<Tensor::IrredTensor<T2, B3, B2, Tensor::DiagonalStructure> >::type iter2;

      result_type Result(x.Basis1(), y.base().Basis1(), QuantumNumbers::QuantumNumber(x.GetSymmetryList()));

      // TODO: this doesn't use the functor - we actually need to call MatrixMatrixMultiplication
      Result.data() = x.data() * herm(y.base().data()); //, f_);

      return Result;
   }

   Functor f_;
};

template <typename T1, typename B1, typename B2, typename S,
          typename T2, typename B3, typename Functor>
struct ScalarProd<HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S> >,
                  Tensor::IrredTensor<T2, B1, B3, S>,
                  Functor>
{
   typedef HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S> > const& first_argument_type;
   typedef Tensor::IrredTensor<T2, B1, B3, S> const& second_argument_type;
   typedef Tensor::IrredTensor<typename result_value<Functor>::type, B2, B3, S> result_type;

   ScalarProd(Functor f = Functor()) : f_(f) {}

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.base().is_null() || y.is_null()) return result_type();

      PRECONDITION_EQUAL(x.base().TransformsAs(), y.TransformsAs());
      DEBUG_PRECONDITION_EQUAL(x.base().Basis1(), y.Basis1());

      typedef typename const_iterator<Tensor::IrredTensor<T1, B1, B2, S> >::type iter1;
      typedef typename const_iterator<Tensor::IrredTensor<T2, B1, B3, S> >::type iter2;

      typedef typename const_inner_iterator<Tensor::IrredTensor<T1, B1, B2, S> >::type inner1;
      typedef typename const_inner_iterator<Tensor::IrredTensor<T2, B1, B3, S> >::type inner2;

      result_type Result(x.base().Basis2(), y.Basis2(), QuantumNumbers::QuantumNumber(y.GetSymmetryList()));

      iter1 xI = iterate(x.base());
      iter2 yI = iterate(y);
      for ( ; xI; ++xI, ++yI)
      {
         for (inner1 xJ = iterate(xI); xJ; ++xJ)
         {
            for (inner2 yJ = iterate(yI); yJ; ++yJ)
            {
	       double InnerDegree = degree(y.Basis1()[yJ.index1()]);
               // accept only scalar elements (diagonal in the quantum number)
               if (Result.Basis1()[xJ.index2()] == Result.Basis2()[yJ.index2()])
	       {
                  add_element(Result.data(), xJ.index2(), yJ.index2(),
                              (InnerDegree / degree(Result.Basis1()[xJ.index2()])) * f_(herm(*xJ), *yJ));
	       }
            }
         }
      }
      DEBUG_CHECK(!std::isnan(norm_frob_sq(Result)));
      return Result;
   }

   Functor f_;
};


//
// adjoint and inv_adjoint
//

template <typename T, typename B1, typename B2, typename S>
struct Adjoint<Tensor::IrredTensor<T, B1, B2, S> > : TensorAdjoint<Tensor::IrredTensor<T, B1, B2, S> > {};

template <typename T, typename B1, typename B2, typename S>
struct InvAdjoint<Tensor::IrredTensor<T, B1, B2, S> > : TensorInvAdjoint<Tensor::IrredTensor<T, B1, B2, S> > {};

struct AdjointCoeffFunctor
{
   typedef double value_type;

   AdjointCoeffFunctor(Tensor::BasisList const& Basis1, Tensor::BasisList const& Basis2,
                       QuantumNumbers::QuantumNumber const& q)
     : Basis1_(Basis1), Basis2_(Basis2), q_(q) {}

   double operator()(size_type i, size_type j) const
   {
      DEBUG_RANGE_CHECK_OPEN(i, 0U, Basis1_.size());
      DEBUG_RANGE_CHECK_OPEN(j, 0U, Basis2_.size());
      return adjoint_coefficient(Basis2_[j], q_, Basis1_[i]);
   }

   Tensor::BasisList Basis1_, Basis2_;
   QuantumNumbers::QuantumNumber q_;
};

struct InvAdjointCoeffFunctor
{
   typedef double value_type;

   InvAdjointCoeffFunctor(Tensor::BasisList const& Basis1, Tensor::BasisList const& Basis2,
                          QuantumNumbers::QuantumNumber const& q)
     : Basis1_(Basis1), Basis2_(Basis2), q_(q) {}

   double operator()(size_type i, size_type j) const
   {
      DEBUG_RANGE_CHECK_OPEN(i, 0U, Basis1_.size());
      DEBUG_RANGE_CHECK_OPEN(j, 0U, Basis2_.size());
      return inverse_adjoint_coefficient(Basis2_[j], q_, Basis1_[i]);
   }

   Tensor::BasisList Basis1_;
   Tensor::BasisList Basis2_;
   QuantumNumbers::QuantumNumber q_;
};

template <typename T, typename B1, typename B2, typename S, typename F>
struct TensorAdjoint<Tensor::IrredTensor<T, B1, B2, S>, F>
{
   TensorAdjoint(F f = F()) : f_(f) {}

   typedef Tensor::IrredTensor<typename result_value<F>::type, B2, B1> result_type;
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;

   result_type operator()(Tensor::IrredTensor<T, B1, B2, S> const& x) const
   {
      if (x.is_null()) return x;
      QuantumNumbers::QuantumNumber q = x.TransformsAs();
      result_type Result(x.Basis2(), x.Basis1(), adjoint(q));
      Result.data() = coefficient_transpose(x.data(),
                                            AdjointCoeffFunctor(x.Basis1(), x.Basis2(), q),
                                            f_);
      return Result;
   }

   F f_;
};

template <typename T, typename B1, typename B2, typename S, typename F>
struct TensorInvAdjoint<Tensor::IrredTensor<T, B1, B2, S>, F>
{
   TensorInvAdjoint(F f = F()) : f_(f) {}

   typedef Tensor::IrredTensor<typename result_value<F>::type, B2, B1> result_type;
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;

   result_type operator()(Tensor::IrredTensor<T, B1, B2, S> const& x) const
   {
      QuantumNumbers::QuantumNumber q = x.TransformsAs();
      result_type Result(x.Basis2(), x.Basis1(), adjoint(q));
      Result.data() = coefficient_transpose(x.data(),
                                            InvAdjointCoeffFunctor(x.Basis1(), x.Basis2(), q),
                                            f_);
      return Result;
   }

   F f_;
};

// FlipConjugate

template <typename T, typename B1, typename B2, typename S>
struct FlipConjugate<Tensor::IrredTensor<T, B1, B2, S> > : TensorFlipConjugate<Tensor::IrredTensor<T, B1, B2, S> > {};

template <typename T, typename B1, typename B2, typename S, typename F>
struct TensorFlipConjugate<Tensor::IrredTensor<T, B1, B2, S>, F>
{
   TensorFlipConjugate(F f = F()) : f_(f) {}

   typedef Tensor::IrredTensor<typename result_value<F>::type, B1, B2, S> result_type;
   typedef Tensor::IrredTensor<T, B1, B2, S> const& argument_type;

   result_type operator()(Tensor::IrredTensor<T, B1, B2, S> const& x) const
   {
      if (x.is_null()) return x;
      QuantumNumbers::QuantumNumber q = x.TransformsAs();
      result_type Result(adjoint(x.Basis1()), adjoint(x.Basis2()), adjoint(q));
      Result.data() = transform(x.data(), f_);
      return Result;
   }
   F f_;
};

//
// prod
//

struct ProductCoeffFunctor
{
   typedef double result_type;

   ProductCoeffFunctor(Tensor::BasisList const& Basis1,
                       Tensor::BasisList const& Basis2,
                       Tensor::BasisList const& Basis3,
                       QuantumNumbers::QuantumNumber const& q1,
                       QuantumNumbers::QuantumNumber const& q2,
                       QuantumNumbers::QuantumNumber const& q);

   double operator()(size_type i, size_type k, size_type j) const
   {
      DEBUG_RANGE_CHECK_OPEN(i, 0U, Basis1_.size());
      DEBUG_RANGE_CHECK_OPEN(k, 0U, Basis2_.size());
      DEBUG_RANGE_CHECK_OPEN(j, 0U, Basis3_.size());
      return is_transform_target(Basis3_[j], q_, Basis1_[i]) ?
         product_coefficient(q1_, q2_, q_, Basis1_[i], Basis3_[j], Basis2_[k])
         :
         0;
   }

   Tensor::BasisList Basis1_;
   Tensor::BasisList Basis2_;
   Tensor::BasisList Basis3_;
   QuantumNumbers::QuantumNumber q1_, q2_, q_;
};

inline
ProductCoeffFunctor::ProductCoeffFunctor(Tensor::BasisList const& Basis1,
                                         Tensor::BasisList const& Basis2,
                                         Tensor::BasisList const& Basis3,
                                         QuantumNumbers::QuantumNumber const& q1,
                                         QuantumNumbers::QuantumNumber const& q2,
                                         QuantumNumbers::QuantumNumber const& q)
   : Basis1_(Basis1), Basis2_(Basis2), Basis3_(Basis3), q1_(q1), q2_(q2), q_(q)
{
}

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2, typename Nest>
struct IrredProd<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2, B2, B3, S2>, Nest>
{
   typedef Tensor::IrredTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef Tensor::IrredTensor<T2, B2, B3, S2> const& second_argument_type;
   typedef Nest third_argument_type;
   typedef Tensor::IrredTensor<typename result_value<Nest>::type, B1, B3> result_type;

   IrredProd() {}

   IrredProd(QuantumNumbers::QuantumNumber const& Transform) : Trans_(Transform) {}

   result_type operator()(first_argument_type x, second_argument_type y,
                          Nest F) const
   {
      if (x.is_null() || y.is_null()) return result_type();

      DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.Basis1());
      QuantumNumbers::QuantumNumber Trans = Trans_;

      if (Trans.is_null())
      {
         QuantumNumbers::QuantumNumberList QL =
            transform_targets(x.TransformsAs(), y.TransformsAs());
         CHECK(QL.size() == 1)("Transform product is not specified and not unique")
            (x.TransformsAs())(y.TransformsAs())(QL);
         Trans = QL[0];
      }

      result_type Result(x.Basis1(), y.Basis2(), Trans);

      // early return if TransformsAs is not in the clebsch-gordan expansion of (x*y)
      if (!is_transform_target(x.TransformsAs(), y.TransformsAs(), Trans))
         return Result;

      ProductCoeffFunctor
         ProdHelper(x.Basis1(), x.Basis2(), y.Basis2(),
                    x.TransformsAs(), y.TransformsAs(), Trans);

      double Tol = std::numeric_limits<double>::epsilon() * 10;

      Result.data() = coefficient_multiply_cull(x.data(), y.data(), ProdHelper, Tol, F);

      Result.debug_check_structure();
      return Result;
   }

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return operator()(x, y, Nest());
   }

   QuantumNumbers::QuantumNumber Trans_;
};

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2>
struct Multiplication<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2, B2, B3, S2> >
   : IrredProd<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2, B2, B3, S2> > {};

} // namespace LinearAlgebra

namespace Tensor
{

//
// triple_prod
//

using LinearAlgebra::HermitianProxy;

template <typename T1, typename T2, typename T3>
struct Prod3Type
{
   typedef typename LinearAlgebra::result_value<LinearAlgebra::Multiplication<T1, T2> >::type v1;
   typedef typename LinearAlgebra::result_value<LinearAlgebra::Multiplication<v1, T3> >::type type;
};

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3>
Tensor::IrredTensor<typename Prod3Type<T1, T2, T3>::type, B2, B4>
triple_prod(HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> > const& x,
            Tensor::IrredTensor<T2, B1, B3, S2> const& E,
            Tensor::IrredTensor<T3, B3, B4, S3> const& y,
            QuantumNumber qxy,
            QuantumNumber qEp);

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3>
Tensor::IrredTensor<typename Prod3Type<T1, T2, T3>::type, B2, B4>
triple_prod(HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> > const& x,
            Tensor::IrredTensor<T2, B1, B3, S2> const& E,
            Tensor::IrredTensor<T3, B3, B4, S3> const& y);

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3>
Tensor::IrredTensor<typename Prod3Type<T1, T2, T3>::type, B1, B4>
triple_prod(Tensor::IrredTensor<T1, B1, B2, S1> const& x,
            Tensor::IrredTensor<T2, B2, B3, S2> const& E,
            HermitianProxy<Tensor::IrredTensor<T3, B4, B3, S3> > const& y,
            QuantumNumber qxy,
            QuantumNumber qEp);

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3>
Tensor::IrredTensor<typename Prod3Type<T1, T2, T3>::type, B1, B4>
triple_prod(Tensor::IrredTensor<T1, B1, B2, S1> const& x,
            Tensor::IrredTensor<T2, B2, B3, S2> const& E,
            HermitianProxy<Tensor::IrredTensor<T3, B4, B3, S3> > const& y);

// add_triple_prod

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3, typename T4, typename S4>
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> > const& x,
                Tensor::IrredTensor<T2, B1, B3, S2> const& E,
                Tensor::IrredTensor<T3, B3, B4, S3> const& y,
                QuantumNumber qxy,
                QuantumNumber qEp);

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3, typename T4, typename S4>
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> > const& x,
                Tensor::IrredTensor<T2, B1, B3, S2> const& E,
                Tensor::IrredTensor<T3, B3, B4, S3> const& y);

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3, typename T4, typename S4>
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                Tensor::IrredTensor<T1, B1, B2, S1> const& x,
                Tensor::IrredTensor<T2, B2, B3, S2> const& E,
                HermitianProxy<Tensor::IrredTensor<T3, B4, B3, S3> > const& y,
                QuantumNumber qxy,
                QuantumNumber qEp);

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2,
          typename T3, typename B4, typename S3, typename T4, typename S4>
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                Tensor::IrredTensor<T1, B1, B2, S1> const& x,
                Tensor::IrredTensor<T2, B2, B3, S2> const& E,
                HermitianProxy<Tensor::IrredTensor<T3, B4, B3, S3> > const& y);

//
// delta_shift
// Shifts the basis by some delta factor, with normalization.
//

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
delta_shift(Tensor::IrredTensor<T, B1, B2, S> const& x,
            QuantumNumbers::QuantumNumber q,
            QuantumNumbers::Projection p);

// Shift the basis by a strictly Abelian shift.
// Precondition: degree(q) == 1
template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
delta_shift(Tensor::IrredTensor<T, B1, B2, S> const& x,
            QuantumNumbers::QuantumNumber q);

// The basis here must have been obtained with the DeltaShift() function for the
// basis.  For use when the shifted basis is already constructed.
template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
delta_shift(Tensor::IrredTensor<T, B1, B2, S> const& x,
            QuantumNumbers::QuantumNumber q,
            QuantumNumbers::Projection p,
            B1 const& NewBasis1,
            B2 const& NewBasis2);

//
// prod variants for HermitianProxy.  The hermitian argument is
// required to transform as a scalar.
//

template <typename T1, typename T2>
struct IrredProd_Herm {};

template <typename T1, typename B1, typename B2, typename S1,
           typename T2, typename B3, typename S2>
struct IrredProd_Herm<HermitianProxy<IrredTensor<T1, B1, B2, S1> >,
                      IrredTensor<T2, B1, B3, S2> >
{
   typedef typename LinearAlgebra::result_value<
      LinearAlgebra::Multiplication<T1, T2> >::type result_value_type;

   typedef IrredTensor<result_value_type, B2, B3> result_type;
   typedef HermitianProxy<IrredTensor<T1, B1, B2, S1> > const& first_argument_type;
   typedef IrredTensor<T2, B1, B3, S2> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

template <typename T1, typename B1, typename B2, typename S1,
           typename T2, typename B3, typename S2>
struct IrredProd_Herm<IrredTensor<T2, B1, B2, S2>,
                      HermitianProxy<IrredTensor<T1, B3, B2, S1> > >
{
   typedef typename LinearAlgebra::result_value<
      LinearAlgebra::Multiplication<T1, T2> >::type result_value_type;

   typedef IrredTensor<result_value_type, B1, B3> result_type;
   typedef IrredTensor<T2, B1, B2, S2> const& first_argument_type;
   typedef HermitianProxy<IrredTensor<T1, B3, B2, S1> > const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};



} // namespace Tensor

namespace LinearAlgebra
{

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2>
struct Multiplication<Tensor::IrredTensor<T1, B1, B2, S1>,
                      HermitianProxy<Tensor::IrredTensor<T2, B3, B2, S2> > >
   : Tensor::IrredProd_Herm<Tensor::IrredTensor<T1, B1, B2, S1>,
                  HermitianProxy<Tensor::IrredTensor<T2, B3, B2, S2> > > {};

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename S2>
struct Multiplication<HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> >,
                      Tensor::IrredTensor<T2, B1, B3, S2> >
   : Tensor::IrredProd_Herm<HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> >,
                  Tensor::IrredTensor<T2, B1, B3, S2> > {};

} // namespace LinearAlgebra

namespace Tensor
{

//
// outer_diagonal
// I've forgotten what this does.  It looks like an attempt to construct the diagonal components
// of the operator x \otimes herm(y), but it cannot work as written, the result should always
// be scalar, for a start.
//

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
outer_diagonal(Tensor::IrredTensor<T, B1, B1, S> const& x,
               HermitianProxy<Tensor::IrredTensor<T, B2, B2, S> > const& y);

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
outer_diagonal(Tensor::IrredTensor<T, B1, B1, S> const& x,
               HermitianProxy<Tensor::IrredTensor<T, B2, B2, S> > const& y)
{
   DEBUG_PRECONDITION_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_PRECONDITION_EQUAL(y.base().Basis1(), y.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(x.TransformsAs(), y.base().TransformsAs());

   typedef typename inner_iterator<Tensor::IrredTensor<T, B1, B1, S> >::type Inner1Iter;
   typedef typename inner_iterator<Tensor::IrredTensor<T, B2, B2, S> >::type Inner2Iter;
   Tensor::IrredTensor<T, B1, B2, S> Result(x.Basis1(), y.Basis1());
   for (int i = 0; i < x.Basis1().size(); ++i)
   {
      Inner1Iter I = iterate_at(x.data(), i, i);
      if (I)
      {
         for (int j = 0; j < y.base().Basis1().size(); ++j)
         {
            Inner2Iter J = iterate_at(y.base().data(), j, j);
            if (J)
            {
               set_element(Result.data(), i, j, outer_prod(diagonal(*I), herm(diagonal(*J))));
            }
         }
      }
   }
   return Result;
}

} // namespace Tensor

#include "tensor.cc"

#endif
