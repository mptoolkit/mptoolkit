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

#if !defined(MPTOOLKIT_TENSOR_TENSOR_H)
#define MPTOOLKIT_TENSOR_TENSOR_H

#include "common/types.h"
#include "basis.h"
#include "quantumnumbers/quantumnumber.h"
#include "pstream/pstream.h"
#include <boost/mpl/assert.hpp>
#include "blas/sparsematrix.h"
#include "blas/arena.h"
#include "blas/number_traits.h"
#include "blas/diagonalmatrix.h"
#include <cmath>

using real_type = double;
using complex_type = std::complex<double>;

namespace Tensor
{

// A proxy class to represent the Hermitian conjugate of a tensor.
template <typename T>
struct HermitianProxy
{
   using base_type  = T;
   using value_type = typename T::value_type;
   using reference  = T const&;

   explicit HermitianProxy(T const& x) : x_(x) {}

   reference base() const { return x_; }

   private:
      reference x_;
};

// A proxy class to represent the tensor conjugate.
template <typename T>
struct ConjugateProxy
{
   using base_type  = T;
   using value_type = typename T::value_type;
   using reference  = T const&;

   explicit ConjugateProxy(T const& x) : x_(x) {}

   reference base() const { return x_; }

   private:
      reference x_;
};

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::SymmetryList;

template <typename T>
struct TagOf
{
   using type = typename T::tag_type;
};

template <>
struct TagOf<double>
{
   using type = blas::cpu_tag;
};

template <>
struct TagOf<std::complex<double>>
{
   using type = blas::cpu_tag;
};

// Structure traits types
struct DefaultStructure
{
   template <typename T>
   using value = blas::SparseMatrix<T>;

   template <typename T>
   using tag_type = typename TagOf<T>::type;
};

struct DiagonalStructure
{
   template <typename T>
   using value = blas::DiagonalMatrix<T>;

   template <typename T>
   using tag_type = typename TagOf<T>::type;
};

// For operations that act on the structure and potentially modify it,
// we need to map the transformed concrete type back to a structure trait.
// In principle this allows us to generalize operations to mutable structures,
// although currently none are implemented.

template <typename T>
struct StructureOf {};

template <typename T>
struct StructureOf<blas::SparseMatrix<T>>
{
   using type = DefaultStructure;
};

template <typename T>
struct StructureOf<blas::DiagonalMatrix<T>>
{
   using type = DiagonalStructure;
};

template <typename T,
          typename Basis1T = BasisList,
          typename Basis2T = Basis1T,
          typename Structure = DefaultStructure >
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

// forward declaration so we can make it a friend
template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
copy(IrredTensor<T, B1, B2, S> const& x);

template <typename T, typename Basis1T, typename Basis2T, typename Structure>
class IrredTensor
{
   public:
      using value_type     = T;
      using structure_type = Structure;

      using StructureType  = typename Structure::template value<T>;

      using numeric_type   = typename StructureType::value_type;

      //      using row_type       = typename StructureType::row_type;
      //using iterator       = typename StructureType::iterator;
      //using const_iterator = typename StructureType::const_iterator;

      using tag_type       = typename Structure::template tag_type<T>;

      typedef Basis1T basis1_type;
      typedef Basis2T basis2_type;

      IrredTensor() = default;

      IrredTensor(IrredTensor const& Other) = delete;

      IrredTensor(IrredTensor&& Other) noexcept = default;

      IrredTensor(basis1_type const& Basis, QuantumNumber const& Trans);

      IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2,
                  QuantumNumber const& Trans);

      // In this variant, the quantum number defaults to the identity
      IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2);


      IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2,
                  QuantumNumber const& Trans, StructureType&& Data);

      explicit IrredTensor(basis1_type const& Basis);

#if 0
      template <typename U, typename US>
      IrredTensor(IrredTensor<U, Basis1T, Basis2T, US> const& r)
         : Basis1_(r.Basis1()), Basis2_(r.Basis2()), Trans_(r.TransformsAs()),
           Data_(r.data()) {}
#endif

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

      // removes all elements, doesn't change the basis or TransformsAs
      void clear() { Data_.clear(); }

      int nnz() const { return Data_.nnz(); }

      typename StructureType::iterator begin() noexcept { return Data_.begin(); }
      typename StructureType::const_iterator begin() const noexcept { return Data_.begin(); }
      typename StructureType::const_iterator cbegin() const noexcept { return Data_.begin(); }

      typename StructureType::iterator end() noexcept { return Data_.end(); }
      typename StructureType::const_iterator end() const noexcept { return Data_.end(); }
      typename StructureType::const_iterator cend() const noexcept { return Data_.end(); }

      typename StructureType::row_type& operator[](int r) { return Data_[r]; }
      typename StructureType::row_type const& operator[](int r) const { return Data_[r]; }

      IrredTensor& operator+=(IrredTensor const& Op);
      IrredTensor& operator+=(IrredTensor&& Op);

      IrredTensor& operator-=(IrredTensor const& Op);
      IrredTensor& operator-=(IrredTensor&& Op);

      template <typename U>
      IrredTensor& operator*=(U const& x) { Data_ *= x; return *this; }

      basis1_type const& Basis1() const { return Basis1_; }
      basis2_type const& Basis2() const { return Basis2_; }

      SymmetryList GetSymmetryList() const { return Basis1_.GetSymmetryList(); }

      template<typename... Args>
      void emplace(int Row, int Col, Args&&... args)
      {
         Data_.emplace(Row, Col, std::forward<Args>(args)...);
      }

      void insert(int r, int c, T const& value)
      {
         Data_.insert(r, c, value);
      }

      template <typename U>
      void insert(int r, int c, U&& value)
      {
         Data_.insert(r, c, std::move(value));
      }

      void add(int r, int c, T const& value)
      {
         Data_.add(r, c, value);
      }

      template <typename U>
      void add(int r, int c, U&& value)
      {
         Data_.add(r, c, std::move(value));
      }

      void subtract(int r, int c, T const& value)
      {
         Data_.subtract(r, c, value);
      }

      template <typename U>
      void subtract(int r, int c, U&& value)
      {
         Data_.subtract(r, c, std::move(value));
      }

      int size1() const { return Basis1_.size(); }
      int size2() const { return Basis2_.size(); }

      value_type operator()(int r, int c) const
      {
	 auto i = Data_.row(r).find(c);
	 if (i == Data_.row(r).end())
	    return value_type();
	 // else
	 return i.value();
      }

      StructureType& data() { return Data_; }
      StructureType const& data() const { return Data_; }

      QuantumNumber const& TransformsAs() const { return Trans_; }

      QuantumNumber const& qn1(int i) const { return Basis1_[i]; }
      QuantumNumber const& qn2(int j) const { return Basis2_[j]; }

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

      static_assert(std::is_nothrow_move_constructible<basis1_type>::value, "");
      static_assert(std::is_nothrow_move_constructible<basis2_type>::value, "");
      static_assert(std::is_nothrow_move_constructible<QuantumNumber>::value, "");
      static_assert(std::is_nothrow_move_constructible<StructureType>::value, "");

  private:
      basis1_type Basis1_;
      basis2_type Basis2_;
      QuantumNumber Trans_;
      StructureType Data_;

      friend PStream::opstream& operator<< <>(PStream::opstream& out, IrredTensor const& x);
      friend PStream::ipstream& operator>> <>(PStream::ipstream& in, IrredTensor& x);

      template <typename T2, typename B12, typename B22, typename S2>
      friend IrredTensor<T2, B12, B22, S2>
      CoerceSymmetryList(IrredTensor<T2, B12, B22, S2> const& t, SymmetryList const& sl);

      friend IrredTensor copy<>(IrredTensor const& x);

      //      friend IrredTensor ::Tensor::CoerceSymmetryList<>(IrredTensor const& t, SymmetryList const& sl);
};

} // namespace Tensor

template <typename T, typename B1, typename B2, typename S>
struct ScalarTypes<Tensor::IrredTensor<T, B1, B2, S>> : ScalarTypes<T> {};

namespace Tensor
{

template <typename T, typename U, typename B1, typename B2, typename S>
void
set(IrredTensor<T, B1, B2, S>& x, IrredTensor<U, B1, B2, S> const& y)
{
   DEBUG_CHECK_EQUAL(x.Basis1(), y.Basis1());
   DEBUG_CHECK_EQUAL(x.Basis2(), y.Basis2());
   DEBUG_CHECK_EQUAL(x.TransformsAs(), y.TransformsAs());
   x.clear();
   for (auto const& ry : y)
   {
      for (auto const& cy : ry)
      {
         T Temp(cy.value.rows(), cy.value.cols());
         set(Temp, cy.value);
         x.insert(ry.row(), cy.col(), std::move(Temp));
      }
   }
   x.check_structure();
}

template <typename T, typename B1, typename B2, typename S>
HermitianProxy<IrredTensor<T, B1, B2, S>>
herm(IrredTensor<T, B1, B2, S> const& x)
{
   return HermitianProxy<IrredTensor<T, B1, B2, S>>(x);
}

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
      Result.insert(i,i, Tensor::make_identity<T>(b, i));
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

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
copy(IrredTensor<T, B1, B2, S> const& x, blas::arena const& A)
{
   return IrredTensor<T, B1, B2, S>(x.Basis1(), x.Basis2(), x.TransformsAs(), copy(x.data(), A));
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
copy(IrredTensor<T, B1, B2, S> const& x)
{
   return IrredTensor<T, B1, B2, S>(x.Basis1(), x.Basis2(), x.TransformsAs(), copy(x.data()));
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
copy(IrredTensor<T, B1, B2, S>&& x)
{
   return x;
}

// arithmetic

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
operator-(IrredTensor<T, B1, B2, S>&& x)
{
   return IrredTensor<T, B1, B2, S>(x.Basis1(), x.Basis2(), x.TransformsAs(), -std::move(x.data()));
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
operator+(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y)
{
   if (x.is_null()) return copy(y);
   if (y.is_null()) return copy(x);
   PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
   PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
   PRECONDITION_EQUAL(x.TransformsAs(), y.TransformsAs());
   return IrredTensor<T, B1, B2, S>(x.Basis1(), x.Basis2(), x.TransformsAs(), x.data() + y.data());
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
operator-(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y)
{
   if (x.is_null()) return copy(y);
   if (y.is_null()) return copy(x);
   PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
   PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
   PRECONDITION_EQUAL(x.TransformsAs(), y.TransformsAs());
   return IrredTensor<T, B1, B2, S>(x.Basis1(), x.Basis2(), x.TransformsAs(), x.data() + y.data());
}

// multiply by scalar

template <typename T, typename B1, typename B2, typename S>
void
inplace_conj(IrredTensor<T, B1, B2, S>& x)
{
   inplace_conj(x.data());
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
conj(IrredTensor<T, B1, B2, S> x)
{
   inplace_conj(x);
   return x;
}

// trace
// we split the implementation into scalar and compound types.
template <typename T, typename B, typename S>
std::enable_if_t<!blas::is_numeric_v<T>, blas::remove_proxy_t<decltype(trace(std::declval<T>()))>>
trace(IrredTensor<T, B, B, S> const& x);

template <typename T, typename B, typename S>
std::enable_if_t<blas::is_numeric_v<T>, T>
trace(IrredTensor<T, B, B, S> const& x);

// norm_frob

template <typename T, typename B1, typename B2, typename S>
decltype(norm_frob(std::declval<T>()))
norm_frob_sq(IrredTensor<T, B1, B2, S> const& x)
{
   decltype(norm_frob(std::declval<T>())) Result{};
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
         Result += qdim(x.qn1(r.row())) * norm_frob_sq(c.value);
      }
   }
}

template <typename T, typename B1, typename B2, typename S>
auto
norm_frob(IrredTensor<T, B1, B2, S> const& x)
{
   return std::sqrt(norm_frob_sq(x));
}

// inner_prod

template <typename T, typename B1, typename B2, typename S>
typename IrredTensor<T, B1, B2, S>::numeric_type
inner_prod(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y)
{
   PRECONDITION_EQUAL(x.TransformsAs(), y.TransformsAs());
   PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
   PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
   typename IrredTensor<T, B1, B2, S>::numeric_type Result = 0.0;
   for (int r = 0; r < x.data().rows(); ++r)
   {
      Result += qdim(x.qn1(r)) * inner_prod(x.data()[r], y.data()[r]);
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
void
inner_prod(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y,
           typename IrredTensor<T, B1, B2, S>::tag_type::template async_ref<complex> Result)
{
   PRECONDITION_EQUAL(x.TransformsAs(), y.TransformsAs());
   PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
   PRECONDITION_EQUAL(x.Basis2(), y.Basis2());

   using tag_type = typename IrredTensor<T, B1, B2, S>::tag_type;

   blas::Vector<complex, tag_type> Temp(x.data().rows());
   blas::Vector<complex> QDim(x.data().rows());
   for (int r = 0; r < x.data().rows(); ++r)
   {
      inner_prod(x.data()[r], y.data()[r], Temp[r]);
      QDim[r] = qdim(x.qn1(r));
   }
   inner_prod(blas::Vector<complex, tag_type>(std::move(QDim)), Temp, Result);
}

// scalar_prod - this is a new function

template <typename T, typename B1, typename B2, typename B3>
IrredTensor<T, B2, B3, Tensor::DefaultStructure>
scalar_prod(HermitianProxy<IrredTensor<T, B1, B2, Tensor::DefaultStructure>> const& x,
            IrredTensor<T, B1, B3, Tensor::DefaultStructure> const& y)
{
   using blas::herm;
   using result_type = IrredTensor<T, B2, B3, Tensor::DefaultStructure>;

   //if (x.is_null() || y.base().is_null()) return result_type();

   result_type Result(x.base().Basis2(), y.Basis2());

   //Result(a,b) = sum_r x(r,a) y(r,b)

   for (auto const& rx : x.base())
   {
      for (auto const& cx : rx)
      {
         for (auto const& cy : y[rx.row()])
         {
            // only diagonal components
            if (Result.qn1(cx.col()) == Result.qn2(cy.col()))
            {
               Result.add(cx.col(), cy.col(), (qdim(y.qn2(rx.row())) / qdim(Result.qn1(cx.col()))) *
                          herm(cx.value) * cy.value);
            }
         }
      }
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename B3>
IrredTensor<T, B2, B3, Tensor::DefaultStructure>
scalar_prod(IrredTensor<T, B1, B3, Tensor::DefaultStructure> const& x,
            HermitianProxy<IrredTensor<T, B1, B2, Tensor::DefaultStructure>> const& y)
{
   using blas::herm;

   using result_type = IrredTensor<T, B2, B3, Tensor::DefaultStructure>;

   //if (x.is_null() || y.base().is_null()) return result_type();

   result_type Result(x.Basis1(), y.base().Basis1());

   // Result(a,b) = sum_c x(a,c) * y(b,c)
   for (auto const& rx : x.data())
   {
      for (auto const& cx : rx)
      {
         for (auto const& ry : y.base().data())
         {
            auto cy = ry.find(cx.col());
            if (cy != ry.end())
            {
               Result.add(rx.row(), ry.row(), cx.value * herm(cy.value()));
            }
         }
      }
   }
   return Result;
}

//
// adjoint and inv_adjoint
//

template <typename T>
std::enable_if_t<blas::is_numeric_v<T>, T>
adjoint(T x)
{
   using std::conj;
   return conj(x);
}

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
adjoint(Tensor::IrredTensor<T, B1, B2, S> const& x)
{
   //   if (x.is_null()) return Tensor::IrredTensor<T, B1, B2, S>();
   QuantumNumbers::QuantumNumber q = x.TransformsAs();
   Tensor::IrredTensor<T, B1, B2, S> Result(x.Basis2(), x.Basis1(), adjoint(q));
   for (auto const& rx : x.data())
   {
      for (auto const& cx : rx)
      {
         Result.data().emplace(cx.col(), rx.row(), adjoint_coefficient(x.qn2(cx.col()), q, x.qn1(rx.row())) * adjoint(cx.value));
      }
   }
   return Result;
}

template <typename T>
std::enable_if_t<blas::is_numeric_v<T>, T>
inv_adjoint(T x)
{
   using std::conj;
   return conj(x);
}

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
inv_adjoint(Tensor::IrredTensor<T, B1, B2, S> const& x)
{
   //   if (x.is_null()) return x;
   QuantumNumbers::QuantumNumber q = x.TransformsAs();
   Tensor::IrredTensor<T, B1, B2, S> Result(x.Basis2(), x.Basis1(), adjoint(q));
   for (auto const& rx : x.data())
   {
      for (auto const& cx : rx)
      {
         Result.data().emplace(cx.col(), rx.row(), inverse_adjoint_coefficient(x.qn2(cx.col()), q, x.qn1(rx.row())) * inv_adjoint(cx.value));
      }
   }
   return Result;
}

// Flip conjugation - combination of complex conjugation and basis adjoint

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
flip_conj(Tensor::IrredTensor<T, B1, B2, S> const& x)
{
   //   if (x.is_null()) return Tensor::IrredTensor<T, B1, B2, S>();
   QuantumNumbers::QuantumNumber q = x.TransformsAs();
   Tensor::IrredTensor<T, B1, B2, S> Result(adjoint(x.Basis1()), adjoint(x.Basis2()), adjoint(q));
   Result.data() = conj(copy(x.data()));
   return Result;
}

// for the real case, where conj is a no-op, this does something useful
template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
flip_conj(Tensor::IrredTensor<T, B1, B2, S>&& x)
{
   //   if (x.is_null()) return Tensor::IrredTensor<T, B1, B2, S>();
   QuantumNumbers::QuantumNumber q = x.TransformsAs();
   Tensor::IrredTensor<T, B1, B2, S> Result(adjoint(x.Basis1()), adjoint(x.Basis2()), adjoint(q));
   Result.data() = conj(std::move(x.data()));
   return Result;
}

//
// prod
//

template <typename T, typename B1, typename B2, typename B3, typename S>
IrredTensor<T, B1, B3, S>
prod(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B2, B3, S> const& y, QuantumNumber const& Trans)
{
   using result_type = IrredTensor<T, B1, B3, S>;
   //   if (x.is_null() || y.is_null()) return result_type();

   DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.Basis1());

   IrredTensor<T, B1, B3, S> Result(x.Basis1(), y.Basis2(), Trans);

   // early return if TransformsAs is not in the clebsch-gordan expansion of (x*y)
   if (!is_transform_target(x.TransformsAs(), y.TransformsAs(), Trans))
      return Result;

   for (auto const& rx : x.data())
   {
      for (auto const& cx : rx)
      {
         for (auto const& cy : y[cx.col()])
         {
            real_type r = product_coefficient(x.TransformsAs(), y.TransformsAs(), Trans,
                                              x.qn1(rx.row()), y.qn2(cy.col()), x.qn2(cx.col()));
            Result.add(rx.row(), cy.col(), r * cx.value * cy.value);
         }
      }
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename B3, typename S, typename U>
void
add_prod(U Factor, IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B2, B3, S> const& y, IrredTensor<T, B1, B3, S>& Result)
{
   DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.Basis1());
   DEBUG_PRECONDITION_EQUAL(x.Basis1(), Result.Basis1());
   DEBUG_PRECONDITION_EQUAL(y.Basis2(), Result.Basis2());

   QuantumNumber Trans = Result.TransformsAs();
   for (auto const& rx : x.data())
   {
      for (auto const& cx : rx)
      {
         for (auto const& cy : y[cx.col()])
         {
            real_type r = product_coefficient(x.TransformsAs(), y.TransformsAs(), Trans,
                                              x.qn1(rx.row()), y.qn2(cy.col()), x.qn2(cx.col()));
            Result.add(rx.row(), cy.col(), Factor * r * cx.value * cy.value);
         }
      }
   }
}

template <typename T, typename B1, typename B2, typename B3, typename S>
inline
IrredTensor<T, B1, B3, S>
operator*(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B2, B3, S> const& y)
{
   QuantumNumbers::QuantumNumberList QL =
      transform_targets(x.TransformsAs(), y.TransformsAs());
   CHECK(QL.size() == 1)("Transform product is not specified and not unique")
      (x.TransformsAs())(y.TransformsAs())(QL);
   return prod(x, y, QL[0]);
}

// multiply by scalar
template <typename T, typename B1, typename B2, typename S, typename U>
inline
std::enable_if_t<blas::is_numeric_v<U>, IrredTensor<T, B1, B2, S>>
operator*(IrredTensor<T, B1, B2, S> Result, U const& a)
{
   Result.data() *= a;
   return Result;
}

template <typename T, typename B1, typename B2, typename S, typename U>
inline
std::enable_if_t<blas::is_numeric_v<U>, IrredTensor<T, B1, B2, S>>
operator*(U const& a, IrredTensor<T, B1, B2, S> Result)
{
   Result.data() *= a;
   return Result;
}

// outer_product: outer product of tensor operators.
// We choose among the possible transform_targets the
// quantum number with the largest degree.
template <typename T, typename B1, typename B2, typename B3, typename S>
IrredTensor<T, B1, B3, S>
outer(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B2, B3, S> const& y, QuantumNumber const& Trans)
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

template <typename T, typename B1, typename B2, typename B3, typename S>
IrredTensor<T, B1, B3, S>
cross(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B2, B3, S> const& y, QuantumNumber const& Trans)
{
   CHECK(cross_product_exists(x.TransformsAs(), y.TransformsAs()))
      ("Cross product does not exist for these operators")
      (x.TransformsAs())(y.TransformsAs());

   return cross_product_factor(x.TransformsAs(), y.TransformsAs())
      * prod(x, y, cross_product_transforms_as(x.TransformsAs(), y.TransformsAs()));
}

template <typename T, typename B, typename S>
IrredTensor<T, B, B, S>
pow(S const& x, int n)
{
   if (n == 0)
   {
      CHECK_EQUAL(x.Basis1(), x.Basis2());
      return IrredTensor<T, B, B, S>::make_identity(x.Basis1());
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



//
// triple_prod
//

#if 0

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

#endif

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

} // namespace Tensor

// helper function objects, extends blas/functors.h
namespace blas
{

struct Prod
{
   Prod(QuantumNumbers::QuantumNumber const& q_) : q(q_) {}

   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return prod(x,y,q);
   }

   QuantumNumbers::QuantumNumber q;
};

} // mamespace blas

#include "tensor.icc"

#endif
