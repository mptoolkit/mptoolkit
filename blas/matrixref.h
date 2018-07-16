// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_BLAS_MATRIXREF_H)
#define MPTOOLKIT_BLAS_MATRIXREF_H

#include "vectorref.h"
#include <list>
#include <mutex>
#include <iostream>
#include "safe-conversions.h"
#include "number_traits.h"
// standard blas calls don't involve any types that have associated namespaces, so
// ADL doesn't find standard blas calls in 2-phase lookup.  So we need to include them
// before defining the wrapper functions.
#include "matrix-blas.h"
#include "matrix-lapack.h"
#include "functors.h"

namespace blas
{

// need forward declaration of VectorRef
template <typename ValueType, typename DerivedType, typename Tag>
class VectorRef;

//
// MatrixRef : generic base class for a matrix using expression templates.
// DerivedType is the ultimate derived type, using the CRTP.
// Tag is a label that identifies storage type; matrices with the same
// storage can be freely mixed in expressions.  Typical use is
// for CPU versus GPU storage.
//
// The tag is also used for blas_traits<tag>.
//


template <typename T, typename Tag = cpu_tag>
class Matrix;

template <typename T, typename Tag = cpu_tag>
class DiagonalMatrix;

template <typename ValueType, typename DerivedType, typename Tag>
class MatrixRef
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;
      using remove_proxy_t = Matrix<ValueType, Tag>;

      // default construction and move construction are defined, no copying.

      MatrixRef() = default;
      ~MatrixRef() = default;
      MatrixRef(MatrixRef&& Other) = default;
      MatrixRef(MatrixRef const&) = delete;
      MatrixRef& operator=(MatrixRef&&) = delete;
      MatrixRef& operator=(MatrixRef const&) = delete;

      derived_type& as_derived() { return *static_cast<derived_type*>(this); }
      derived_type const& as_derived() const { return *static_cast<derived_type const*>(this); }

      int rows() const { return this->as_derived().rows(); }
      int cols() const { return this->as_derived().cols(); }
      std::pair<int,int> size() const { return this->as_derived().size(); }
};

// derived class for a diagonal matrix

template <typename ValueType, typename DerivedType, typename Tag>
class DiagonalMatrixRef : public MatrixRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;
      using remove_proxy_t = DiagonalMatrix<ValueType, Tag>;
};

// Specialization of a MatrixRef for a matrix that can be directly used in BLAS-like calls,
// ie it physically exists and can be addressed (although not necessarily in main memory,
// eg it might be on some other device such as a GPU), and has a valid BLAS trans parameter
// (ie, N, T, C, R).
// The BaseType is the concrete type, which allows us to specialize blas calls for different
// devices.
// We expose the leading_dimension() function, but not the memory storage, that needs
// to be done separately for each BaseType.

template <typename ValueType, typename DerivedType, typename Tag>
class BlasMatrix : public MatrixRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      BlasMatrix() = default;
      ~BlasMatrix() = default;
      BlasMatrix(BlasMatrix&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }

      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation)
      char trans() const { return this->as_derived().trans(); }

      storage_type storage() & { return this->as_derived().storage(); }
      const_storage_type storage() const& { return this->as_derived().storage(); }
};

// forward
template <typename ValueType, typename Tag>
class matrix_range_view;

// Denotes an L-value matrix, which is a BlasMatrix with the trans() parameter equal to 'N'
template <typename ValueType, typename DerivedType, typename Tag>
class NormalMatrix : public BlasMatrix<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      NormalMatrix() = default;
      ~NormalMatrix() = default;
      NormalMatrix(NormalMatrix&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }

      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation)
      constexpr char trans() const { return 'N'; }

      matrix_range_view<ValueType, Tag>
      operator()(Range rows, Range cols);

      matrix_range_view<ValueType, Tag>
      operator()(Range rows, all_t);

      matrix_range_view<ValueType, Tag>
      operator()(all_t, Range cols);

      normal_vector_view<ValueType, Tag>
      col(int c)
      {
	 DEBUG_RANGE_CHECK_OPEN(c, 0, this->cols());
	 return normal_vector_view<ValueType, Tag>(this->rows(), this->storage()+c*this->leading_dimension());
      }

      const_normal_vector_view<ValueType, Tag>
      col(int c) const
      {
	 DEBUG_RANGE_CHECK_OPEN(c, 0, this->cols());
	 return const_normal_vector_view<ValueType, Tag>(this->rows(), this->storage()+c*this->leading_dimension());
      }

      vector_view<ValueType, Tag>
      row(int r)
      {
	 DEBUG_RANGE_CHECK_OPEN(r, 0, this->rows());
	 return vector_view<ValueType, Tag>(this->cols(), this->leading_dimension(), this->storage()+r);
      }

      const_vector_view<ValueType, Tag>
      row(int r) const
      {
	 DEBUG_RANGE_CHECK_OPEN(r, 0, this->rows());
	 return vector_view<ValueType, Tag>(this->cols(), this->leading_dimension(), this->storage()+r);
      }

      storage_type storage() { return this->as_derived().storage(); }
      const_storage_type storage() const { return this->as_derived().storage(); }
};

// proxy class for a matrix that can appear on the left-hand side of an expression as an r-value reference

template <typename ValueType, typename DerivedType, typename Tag>
class NormalMatrixProxy : public BlasMatrix<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      NormalMatrixProxy() = default;
      ~NormalMatrixProxy() = default;
      NormalMatrixProxy(NormalMatrixProxy&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }

      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation)
      constexpr char trans() const { return 'N'; }

      matrix_range_view<ValueType, Tag>
      operator()(Range rows, Range cols) &&;

      normal_vector_view<ValueType, Tag>
      col(int c) &&
      {
	 DEBUG_RANGE_CHECK_OPEN(c, 0, this->cols());
	 return normal_vector_view<ValueType, Tag>(this->rows(), this->storage()+c*this->leading_dimension());
      }

      const_normal_vector_view<ValueType, Tag>
      col(int c) const&
      {
	 DEBUG_RANGE_CHECK_OPEN(c, 0, this->cols());
	 return const_normal_vector_view<ValueType, Tag>(this->rows(), this->storage()+c*this->leading_dimension());
      }

      vector_view<ValueType, Tag>
      row(int r) &&
      {
	 DEBUG_RANGE_CHECK_OPEN(r, 0, this->rows());
	 return vector_view<ValueType, Tag>(this->cols(), this->leading_dimension(), this->storage()+r);
      }

      const_vector_view<ValueType, Tag>
      row(int r) const&
      {
	 DEBUG_RANGE_CHECK_OPEN(r, 0, this->rows());
	 return vector_view<ValueType, Tag>(this->cols(), this->leading_dimension(), this->storage()+r);
      }

      storage_type storage() && { return this->as_derived().storage(); }
      const_storage_type storage() const& { return this->as_derived().storage(); }
};

// proxy class for the transpose of a BlasMatrix

template <typename T, typename BaseType, typename Tag>
class BlasMatrixTrans : public BlasMatrix<T, BlasMatrixTrans<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      BlasMatrixTrans(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.cols(); }
      int cols() const { return Base.rows(); }
      std::pair<int,int> size() const { return {Base.cols(), Base.rows()}; }
      int leading_dimension() const { return Base.leading_dimension(); }
      char trans() const { return number_traits<value_type>::blas_trans(Base.trans()); }

      base_type const& base() const { return Base; }

      const_storage_type storage() const& { return Base.storage(); }

   private:
      base_type const& Base;
};

// proxy class for the hermitian conjugate of a BlasMatrix

template <typename T, typename BaseType, typename Tag>
class BlasMatrixHerm : public BlasMatrix<T, BlasMatrixHerm<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      BlasMatrixHerm(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.cols(); }
      int cols() const { return Base.rows(); }
      std::pair<int,int> size() const { return {Base.cols(), Base.rows()}; }
      int leading_dimension() const { return Base.leading_dimension(); }
      char trans() const { return number_traits<value_type>::blas_herm(Base.trans()); }

      base_type const& base() const { return Base; }

      const_storage_type storage() const& { return Base.storage(); }

   private:
      base_type const& Base;
};

// proxy class for the complex conjugate of a BlasMatrix

template <typename T, typename BaseType, typename Tag>
class BlasMatrixConj : public BlasMatrix<T, BlasMatrixConj<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      BlasMatrixConj(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.rows(); }
      int cols() const { return Base.cols(); }
      std::pair<int,int> size() const { return Base.size(); }
      int leading_dimension() const { return Base.leading_dimension(); }
      char trans() const { return number_traits<value_type>::blas_conj(Base.trans()); }
      T* data() { return Base.data(); }
      T const* data() const { return Base.data(); }

      base_type const& base() const { return Base; }

      const_storage_type storage() const& { return Base.storage(); }

   private:
      base_type const& Base;
};

// trans

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixTrans<T, BaseType, Tag>
trans(BlasMatrix<T, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<T, BaseType, Tag>(x.as_derived());
}

template <typename T, typename BaseType, typename Tag>
BaseType const&
trans(BlasMatrixTrans<T, BaseType, Tag> const& x)
{
   return x.base();
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixConj<T, BaseType, Tag>
trans(BlasMatrixHerm<T, BaseType, Tag> const& x)
{
   return BlasMatrixConj<T, BaseType, Tag>(x.base());
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixHerm<T, BaseType, Tag>
trans(BlasMatrixConj<T, BaseType, Tag> const& x)
{
   return BlasMatrixHerm<T, BaseType, Tag>(x.base());
}

// herm

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixHerm<T, BaseType, Tag>
herm(BlasMatrix<T, BaseType, Tag> const& x)
{
   return BlasMatrixHerm<T, BaseType, Tag>(x.as_derived());
}


template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixHerm<T, BaseType, Tag>
conj(BlasMatrixTrans<T, BaseType, Tag> const& x)
{
   return BlasMatrixHerm<T, BaseType, Tag>(x.base());
}

template <typename T, typename BaseType, typename Tag>
inline
BaseType const&
herm(BlasMatrixHerm<T, BaseType, Tag> const& x)
{
   return x.base();
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixTrans<T, BaseType, Tag>
herm(BlasMatrixConj<T, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<T, BaseType, Tag>(x.base());
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixConj<T, BaseType, Tag>
herm(BlasMatrixTrans<T, BaseType, Tag> const& x)
{
   return BlasMatrixConj<T, BaseType, Tag>(x.base());
}

// specialization for real types, herm() reduces to trans()
template <typename BaseType, typename Tag>
inline
BlasMatrixTrans<float, BaseType, Tag>
herm(BlasMatrix<float, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<float, BaseType, Tag>(x.as_derived());
}

template <typename BaseType, typename Tag>
inline
BlasMatrixTrans<double, BaseType, Tag>
herm(BlasMatrix<double, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<double, BaseType, Tag>(x.as_derived());
}

#if defined(HAVE_FLOAT128)
template <typename BaseType, typename Tag>
inline
BlasMatrixTrans<float128, BaseType, Tag>
herm(BlasMatrix<float128, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<float128, BaseType, Tag>(x.as_derived());
}
#endif

// conj

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixConj<T, BaseType, Tag>
conj(BlasMatrix<T, BaseType, Tag> const& x)
{
   return BlasMatrixConj<T, BaseType, Tag>(x.as_derived());
}

// specialization for real types, conj() reduces to a no-op
template <typename BaseType, typename Tag>
inline
BaseType const&
conj(BlasMatrix<float, BaseType, Tag> const& x)
{
   return x.as_derived();
}

template <typename BaseType, typename Tag>
inline
BaseType const&
conj(BlasMatrix<double, BaseType, Tag> const& x)
{
   return x.as_derived();
}

#if defined(HAVE_FLOAT128)
template <typename BaseType, typename Tag>
inline
BaseType const&
conj(BlasMatrix<float128, BaseType, Tag> const& x)
{
   return x.as_derived();
}
#endif

// inplace_conj

template <typename T, typename U, typename Tag>
inline
void inplace_conj(NormalMatrix<T, U, Tag>& A)
{
   matrix_conj(A.rows(), A.cols(), A.storage(), A.leading_dimension());
}

template <typename T, typename U, typename Tag>
inline
void inplace_conj(NormalMatrixProxy<T, U, Tag>&& A)
{
   matrix_conj(A.rows(), A.cols(), std::move(A).storage(), A.leading_dimension());
}

// specalizations for real types are no-ops

template <typename U, typename Tag>
inline
void inplace_conj(NormalMatrix<float, U, Tag>&)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(NormalMatrixProxy<float, U, Tag>&&)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(NormalMatrix<double, U, Tag>&)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(NormalMatrixProxy<double, U, Tag>&&)
{
}

#if defined(HAVE_FLOAT128)
template <typename U, typename Tag>
inline
void inplace_conj(NormalMatrix<float128, U, Tag>&)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(NormalMatrixProxy<float128, U, Tag>&&)
{
}
#endif

// assignment

template <typename T, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_copy(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>&& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_copy(B.as_derived(), std::move(A.as_derived()));
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>& A, DiagonalMatrixRef<V, W, Tag> const& B)
{
   clear(A.as_derived());
   vector_copy(B.as_derived().diagonal(), A.as_derived().diagonal());
   //   A.as_derived().diagonal() = B.diagonal();
   //   matrix_copy(B.as_derived(), std::move(A.as_derived()));
}

// add

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>&& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add(B.as_derived(), std::move(A.as_derived()));
}

template <typename T, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add_scaled(-number_traits<T>::identity(), B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, U, Tag>&& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add_scaled(-number_traits<T>::identity(), B.as_derived(), std::move(A.as_derived()));
}

// expression template for alpha * op(A)

template <typename T, typename BaseType, typename Tag>
class ScaledMatrix : public MatrixRef<T, ScaledMatrix<T, BaseType, Tag>, Tag>
{
   public:
      using base_type = BaseType;

      ScaledMatrix(T const& Factor_, base_type const& A_) : Factor(Factor_), A(A_) {}

      int rows() const { return A.rows(); }
      int cols() const { return A.cols(); }

      base_type const& base() const { return A; }
      T factor() const { return Factor; }

   private:
      T Factor;
      base_type const& A;
};

template <typename T, typename BaseType, typename Tag, typename X>
inline
ScaledMatrix<decltype(safe_convert<T>(std::declval<X>())), BaseType, Tag>
operator*(X const& alpha, MatrixRef<T, BaseType, Tag> const& M)
{
   return ScaledMatrix<T, BaseType, Tag>(safe_convert<T>(alpha), M.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_copy_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>&& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(B.factor(), B.base(), std::move(A.as_derived()));
}

template <typename T, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, U, Tag>& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(-B.factor(), B.base(), A.as_derived());
}

// expression template for -A

template <typename T, typename BaseType, typename Tag>
class NegatedMatrix : public MatrixRef<T, NegatedMatrix<T, BaseType, Tag>, Tag>
{
   public:
      using base_type = BaseType;

      NegatedMatrix(base_type const& A_) : A(A_) {}

      int rows() const { return A.rows(); }
      int cols() const { return A.cols(); }

      base_type const& base() const { return A; }

   private:
      base_type const& A;
};

template <typename T, typename BaseType, typename Tag>
inline
NegatedMatrix<T, BaseType, Tag>
operator-(MatrixRef<T, BaseType, Tag> const& M)
{
   return NegatedMatrix<T, BaseType, Tag>(M.as_derived());
}

template <typename T, typename BaseType, typename Tag>
inline
BaseType const&
operator-(NegatedMatrix<T, BaseType, Tag> const& M)
{
   return M.base();
}

template <typename T, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>& A, NegatedMatrix<T, V, Tag> const& B)
{
   matrix_copy_scaled(-number_traits<T>::identity(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>& A, NegatedMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(-number_traits<T>::identity(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>&& A, NegatedMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(-number_traits<T>::identity(), B.base(), std::move(A.as_derived()));
}

template <typename T, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, U, Tag>& A, NegatedMatrix<T, V, Tag> const& B)
{
   matrix_add(B.base(), A.as_derived());
}

// expression template for alpha * op(A) * op(B)

template <typename T, typename U, typename D1, typename D2, typename Tag>
struct MatrixProduct : public MatrixRef<remove_proxy_t<decltype(std::declval<T>()*std::declval<U>())>, 
					MatrixProduct<T, U, D1, D2, Tag>, Tag>
{
   using value_type = remove_proxy_t<decltype(std::declval<T>()*std::declval<U>())>;

   MatrixProduct(MatrixRef<T, D1, Tag> const& A_, MatrixRef<U, D2, Tag> const& B_)
      : Factor(number_traits<value_type>::identity()), A(A_.as_derived()), B(B_.as_derived()) {}

   MatrixProduct(value_type const& Factor_, MatrixRef<T, D1, Tag> const& A_, 
		 MatrixRef<U, D2, Tag> const& B_)
      : Factor(Factor_), A(A_.as_derived()), B(B_.as_derived()) {}

   int rows() const { return A.rows(); }
   int cols() const { return B.cols(); }
   value_type factor() const { return Factor; }

   value_type Factor;
   D1 const& A;
   D2 const& B;
};

#if 0
template <typename T, typename U, typename V, typename Tag>
typename blas_traits<Tag>::template matrix_type<T>
evaluate(MatrixProduct<T, U, V, Tag>&& x)
{
   typename blas_traits<Tag>::template matrix_type<T> Result(x.rows(), x.cols());
   assign(Result, std::move(x));
   return Result;
}
#endif

template <typename T, typename U, typename D1, typename D2, typename Tag>
inline
MatrixProduct<T, U, D1, D2, Tag>
operator*(MatrixRef<T, D1, Tag> const& A, MatrixRef<U, D2, Tag> const& B)
{
   return MatrixProduct<T, U, D1, D2, Tag>(A.as_derived(), B.as_derived());
}

template <typename T, typename U, typename D1, typename D2, typename Tag>
inline
MatrixProduct<T, U, D1, D2, Tag>
operator*(ScaledMatrix<T, D1, Tag> const& A, MatrixRef<U, D2, Tag> const& B)
{
   return MatrixProduct<T, U, D1, D2, Tag>(A.factor(), A.base(), B.as_derived());
}

// operator* involving a DiagonalMatrix and a Matrix

template <typename T, typename U, typename D1, typename D2, typename Tag>
inline
MatrixProduct<T, U, D1, D2, Tag>
operator*(MatrixRef<T, D1, Tag> const& A, DiagonalMatrixRef<U, D2, Tag> const& B)
{
   return MatrixProduct<T, U, D1, D2, Tag>(A.as_derived(), B.as_derived());
}

template <typename T, typename U, typename D1, typename D2, typename Tag>
inline
MatrixProduct<T, U, D1, D2, Tag>
operator*(DiagonalMatrixRef<T, D1, Tag> const& A, MatrixRef<U, D2, Tag> const& B)
{
   return MatrixProduct<T, U, D1, D2, Tag>(A.as_derived(), B.as_derived());
}


template <typename T, typename U, typename D1, typename D2, typename Tag>
inline
MatrixProduct<T, U, D1, D2, Tag>
operator*(ScaledMatrix<T, D1, Tag> const& A, DiagonalMatrixRef<U, D2, Tag> const& B)
{
   return MatrixProduct<T, U, D1, D2, Tag>(A.factor(), A.base(), B.as_derived());
}

#if 0
// TODO: need a DiagonalMatrixProduct class
template <typename T, typename U, typename D1, typename D2, typename Tag>
inline
MatrixProduct<T, U, D1, D2, Tag>
operator*(DiagonalMatrixRef<T, D1, Tag> const& A, DiagonalMatrixRef<U, D2, Tag> const& B)
{
   return MatrixProduct<T, U, D1, D2, Tag>(A.as_derived(), B.as_derived());
}
#endif

// matrix product with a scalar, alpha * A * B

template <typename T, typename Derived, typename U, typename V, typename D1, typename D2, typename Tag>
inline
void assign(MatrixRef<T, Derived, Tag>& C, MatrixProduct<U, V, D1, D2, Tag> const& a)
{
   gemm(a.Factor, a.A, number_traits<T>::zero(), a.B, C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename D1, typename D2, typename Tag>
inline
void add(MatrixRef<T, Derived, Tag>& C, MatrixProduct<U, V, D1, D2, Tag> const& a)
{
   gemm(a.Factor, a.A.as_derived(), number_traits<T>::identity(), a.B.as_derived(), C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename D1, typename D2, typename Tag>
inline
void subtract(MatrixRef<T, Derived, Tag>& C, MatrixProduct<U, V, D1, D2, Tag> const& a)
{
   gemm(-a.Factor, a.A, number_traits<T>::zero(), a.B, C.as_derived());
}

template <typename T, typename Derived, typename Tag>
void fill(MatrixRef<T, Derived, Tag>& C, T const& x)
{
   matrix_fill(x, C.as_derived());
}

// miscellaneous

template <typename T, typename U, typename Tag>
void
inline
trace(BlasMatrix<T, U, Tag> const& x, typename Tag::template async_proxy<T>&& r)
{
   sum(x.as_derived().diagonal(), std::move(r));
}

template <typename T, typename U, typename Tag>
void
inline
trace(BlasMatrix<T, U, Tag> const& x, typename Tag::template async_ref<T>& r)
{
   sum(x.as_derived().diagonal(), r);
}

template <typename T, typename U, typename Tag>
inline
T
trace(BlasMatrix<T, U, Tag> const& x)
{
   return sum(x.as_derived().diagonal());
}

//
// DiagonalBlasMatrix - a matrix that physically exists in memory with a stride
// Essentially a different view of a BlasVector.
//

template <typename ValueType, typename DerivedType, typename Tag>
class DiagonalBlasMatrix : public DiagonalMatrixRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      DiagonalBlasMatrix() = default;
      ~DiagonalBlasMatrix() = default;
      DiagonalBlasMatrix(DiagonalBlasMatrix&& Other) = default;

      int stride() const { return this->as_derived().stride(); }

      storage_type storage() { return this->as_derived().storage(); }
      const_storage_type storage() const { return this->as_derived().storage(); }
};

// expression template for alpha * op(A) * op(B)

template <typename T, typename U, typename V, typename Tag>
struct MatrixVectorProduct : public VectorRef<T, MatrixVectorProduct<T, U, V, Tag>, Tag>
{
   MatrixVectorProduct(MatrixRef<T, U, Tag> const& A_, VectorRef<T, V, Tag> const& B_)
      : Factor(number_traits<T>::identity()), A(A_.as_derived()), B(B_.as_derived()) {}

   MatrixVectorProduct(T const& Factor_, MatrixRef<T, U, Tag> const& A_, VectorRef<T, V, Tag> const& B_)
      : Factor(Factor_), A(A_.as_derived()), B(B_.as_derived()) {}

   int size() const { return A.rows(); }
   T factor() const { return Factor; }

   T Factor;
   U const& A;
   V const& B;
};

template <typename T, typename U, typename V, typename Tag>
MatrixVectorProduct<T, U, V, Tag>
operator*(MatrixRef<T, U, Tag> const& A, VectorRef<T, V, Tag> const& B)
{
   return MatrixVectorProduct<T, U, V, Tag>(A.as_derived(), B.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
MatrixVectorProduct<T, U, V, Tag>
operator*(ScaledMatrix<T, U, Tag> const& A, VectorRef<T, V, Tag> const& B)
{
   return MatrixVectorProduct<T, U, V, Tag>(A.factor(), A.base(), B.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
MatrixVectorProduct<T, U, V, Tag>
operator*(MatrixRef<T, U, Tag> const& A, ScaledVector<T, V, Tag> const& B)
{
   return MatrixVectorProduct<T, U, V, Tag>(B.factor(), A.as_derived(), B.base());
}

template <typename T, typename U, typename V, typename Tag>
MatrixVectorProduct<T, U, V, Tag>
operator*(ScaledMatrix<T, U, Tag> const& A, ScaledVector<T, V, Tag> const& B)
{
   return MatrixVectorProduct<T, U, V, Tag>(A.factor() * B.factor(), A.base(), B.as_derived());
}

// this overload doesn't appear to be called properly??
template <typename T, typename U, typename V, typename Tag, typename X>
MatrixVectorProduct<decltype(safe_convert<T>(std::declval<X>())), U, V, Tag>
operator*(X const& alpha, MatrixVectorProduct<T, U, V, Tag> const& x)
{
   return MatrixVectorProduct<T, U, V, Tag>(safe_convert<T>(alpha)*x.Factor, x.A, x.B);
}

// assignment involving MatrixVectorProduct

template <typename T, typename Derived, typename U, typename V, typename Tag>
void assign(VectorRef<T, Derived, Tag>& C, MatrixVectorProduct<T, U, V, Tag> const& a)
{
   gemv(a.Factor, a.A, a.B, number_traits<T>::zero(), C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
void assign(BlasVectorProxy<T, Derived, Tag>&& C, MatrixVectorProduct<T, U, V, Tag> const& a)
{
   gemv(a.Factor, a.A, a.B, number_traits<T>::zero(), std::move(C).as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
void add(VectorRef<T, Derived, Tag>& C, MatrixVectorProduct<T, U, V, Tag> const& a)
{
   gemv(a.Factor, a.A, a.B, number_traits<T>::identity(), C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
void subtract(VectorRef<T, Derived, Tag>& C, MatrixVectorProduct<T, U, V, Tag> const& a)
{
   gemv(-a.Factor, a.A, a.B, number_traits<T>::identity(), C.as_derived());
}

//
// VECTOR middle-layer BLAS wrappers, that forward from a matrix/vector ref to low-level storage
//

template <typename T, typename U, typename Tag>
inline
void clear(BlasVector<T, U, Tag>& C)
{
   vector_clear(C.size, C.storage(), C.stride());
}

template <typename T, typename U, typename Tag>
inline
void clear(BlasVectorProxy<T, U, Tag>&& C)
{
   vector_clear(C.size, std::move(C).storage(), C.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_copy_scaled(T alpha, blas::BlasVector<T, U, Tag> const& x, BlasVector<T, V, Tag>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_copy_scaled(x.size(), alpha, x.storage(), x.stride(), y.storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_copy_scaled(T alpha, blas::BlasVector<T, U, Tag> const& x, BlasVectorProxy<T, V, Tag>&& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_copy_scaled(x.size(), alpha, x.storage(), x.stride(), std::move(y).storage(), y.stride());
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void vector_copy(blas::BlasVector<T, U, Tag> const& x, BlasVector<V, W, Tag>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_copy(x.size(), x.storage(), x.stride(), y.storage(), y.stride());
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void vector_copy(blas::BlasVector<T, U, Tag> const& x, BlasVectorProxy<V, W, Tag>&& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_copy(x.size(), x.storage(), x.stride(), std::move(y).storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add_scaled(T alpha, blas::BlasVector<T, U, Tag> const& x, BlasVector<T, V, Tag>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_add_scaled(x.size(), alpha, x.storage(), x.stride(), y.storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add_scaled(T alpha, blas::BlasVector<T, U, Tag> const& x, BlasVectorProxy<T, V, Tag>&& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_add_scaled(x.size(), alpha, x.storage(), x.stride(), std::move(y).storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add(blas::BlasVector<T, U, Tag> const& x, BlasVector<T, V, Tag>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_add(x.size(), blas::number_traits<T>::identity(),
              x.storage(), x.stride(),
              y.storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add(blas::BlasVector<T, U, Tag> const& x, BlasVectorProxy<T, V, Tag>&& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_add(x.size(), blas::number_traits<T>::identity(),
              x.storage(), x.stride(),
              y.storage(), y.stride());
}

template <typename T, typename U, typename Tag>
inline
void vector_scale(T alpha, BlasVector<T, U, Tag>& y)
{
   vector_scale(y.size(), alpha, y.storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void scale(BlasVector<T, U, Tag>& C, V const& x)
{
   vector_scale(C.size(), x, C.storage(), C.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void scale(BlasVectorProxy<T, U, Tag>&& C, V const& x)
{
   vector_scale(C.size(), x, std::move(C).storage(), C.stride());
}

// TODO: these middle-layer functions should disappear, eg this is just an overload of scale()
template <typename T, typename U, typename Tag>
inline
void vector_scale(T alpha, BlasVectorProxy<T, U, Tag>&& y)
{
   vector_scale(y.size(), alpha, std::move(y).storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_parallel(blas::BlasVector<T, U, Tag> const& x, blas::BlasVector<T, U, Tag> const& y,
		     BlasVector<T, V, Tag>& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_parallel(x.size(),
		   x.storage(), x.stride(),
		   y.storage(), y.stride(),
		   z.storage(), z.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_parallel(blas::BlasVector<T, U, Tag> const& x, blas::BlasVector<T, U, Tag> const& y,
		     BlasVectorProxy<T, V, Tag>&& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_parallel(x.size(),
		   x.storage(), x.stride(),
		   y.storage(), y.stride(),
		   std::move(z).storage(), z.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add_parallel(blas::BlasVector<T, U, Tag> const& x, blas::BlasVector<T, U, Tag> const& y,
			 BlasVector<T, V, Tag>& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_add_parallel(x.size(),
		       x.storage(), x.stride(),
		       y.storage(), y.stride(),
		       z.storage(), z.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add_parallel(blas::BlasVector<T, U, Tag> const& x, blas::BlasVector<T, U, Tag> const& y,
			 BlasVectorProxy<T, V, Tag>&& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_add_parallel(x.size(),
		       x.storage(), x.stride(),
		       y.storage(), y.stride(),
		       std::move(z).storage(), z.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_parallel_scaled(T const& scale, BlasVector<T, U, Tag> const& x, 
			    BlasVector<T, U, Tag> const& y,
			    BlasVector<T, V, Tag>& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_parallel_scaled(x.size(), scale,
			  x.storage(), x.stride(),
			  y.storage(), y.stride(),
			  z.storage(), z.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_parallel_scaled(T const& scale, BlasVector<T, U, Tag> const& x, 
			    BlasVector<T, U, Tag> const& y,
			    BlasVectorProxy<T, V, Tag>&& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_parallel_scaled(x.size(), scale,
			  x.storage(), x.stride(),
			  y.storage(), y.stride(),
			  std::move(z).storage(), z.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add_parallel_scaled(T const& scale, BlasVector<T, U, Tag> const& x, 
				BlasVector<T, U, Tag> const& y,
				BlasVector<T, V, Tag>& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_add_parallel_scaled(x.size(), scale,
			      x.storage(), x.stride(),
			      y.storage(), y.stride(),
			      z.storage(), z.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_add_parallel_scaled(T const& scale, BlasVector<T, U, Tag> const& x, 
				BlasVector<T, U, Tag> const& y,
				BlasVectorProxy<T, V, Tag>&& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   DEBUG_CHECK_EQUAL(x.size(), z.size());
   vector_add_parallel_scaled(x.size(), scale,
			      x.storage(), x.stride(),
			      y.storage(), y.stride(),
			      std::move(z).storage(), z.stride());
}

template <typename T, typename U, typename Tag>
void
vector_fill(T alpha, blas::BlasVector<T, U, Tag>& y)
{
   vector_fill(alpha, y.size(), y.storage(), y.stride());
}

template <typename T, typename U, typename Tag>
void
vector_fill(T alpha, blas::BlasVectorProxy<T, U, Tag>&& y)
{
   vector_fill(alpha, y.size(), std::move(y).storage(), y.stride());
}

template <typename T, typename U, typename Tag>
void
sum(blas::BlasVector<T, U, Tag> const& x, typename Tag::template async_proxy<T>&& y)
{
   vector_sum(x.size(), x.storage(), x.stride(), std::move(y));
}

template <typename T, typename U, typename Tag>
void
sum(blas::BlasVector<T, U, Tag> const& x, typename Tag::template async_ref<T>& y)
{
   vector_sum(x.size(), x.storage(), x.stride(), y);
}

template <typename T, typename U, typename Tag>
std::remove_reference_t<typename Tag::template async_ref<T>>
sum(blas::BlasVector<T, U, Tag> const& x)
{
   //   std::remove_reference_t<typename Tag::template async_ref<T>> y{};
   auto y = Tag::template allocate_async_ref<T>();
   vector_sum(x.size(), x.storage(), x.stride(), y);
   return y;
}

template <typename T, typename U, typename Tag>
void
norm_frob_sq(blas::BlasVector<T, U, Tag> const& x, typename Tag::template async_proxy<T>&& y)
{
   vector_norm_frob_sq(x.size(), x.storage(), x.stride(), std::move(y));
}

template <typename T, typename U, typename Tag>
void
norm_frob_sq(blas::BlasVector<T, U, Tag> const& x, typename Tag::template async_ref<T>& y)
{
   vector_norm_frob_sq(x.size(), x.storage(), x.stride(), y);
}

// inplace_conj

template <typename T, typename U, typename Tag>
inline
void inplace_conj(BlasVector<T, U, Tag>& x)
{
   vector_conj(x.size(), x.storage(),x.stride());
}

template <typename T, typename U, typename Tag>
inline
void inplace_conj(BlasVectorProxy<T, U, Tag>&& x)
{
   vector_conj(x.size(), std::move(x).storage(), x.stride());
}

// inplace_conj is a no-op for real types

template <typename U, typename Tag>
inline
void inplace_conj(BlasVector<float, U, Tag>& x)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(BlasVectorProxy<float, U, Tag>&& x)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(BlasVector<double, U, Tag>& x)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(BlasVectorProxy<double, U, Tag>&& x)
{
}

#if defined(HAVE_FLOAT128)
template <typename U, typename Tag>
inline
void inplace_conj(BlasVector<float128, U, Tag>& x)
{
}

template <typename U, typename Tag>
inline
void inplace_conj(BlasVectorProxy<float128, U, Tag>&& x)
{
}
#endif

#if 0
template <typename T, typename U, typename V, typename Tag>
inline
void
inner_prod(blas::BlasVector<T, U, Tag> const& x,
	   blas::BlasVector<T, V, Tag> const& y,
	   typename Tag::template async_ref<T>& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_inner_prod(x.size(), x.storage(), x.stride(), y.storage(), y.stride(), z);
}
#endif

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
void
inner_prod_nested(blas::BlasVector<T, U, Tag> const& x,
		  blas::BlasVector<T, V, Tag> const& y,
		  typename Tag::template async_ref<T>& z,
		  Nested&& Nest)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_inner_prod_nested(x.size(), x.storage(), x.stride(), y.storage(), y.stride(), z,
			    std::forward<Nested>(Nest));
}

#if 0
template <typename T, typename U, typename V, typename Tag>
inline
void
inner_prod(blas::BlasVector<T, U, Tag> const& x,
	   blas::BlasVector<T, V, Tag> const& y,
	   typename Tag::template async_proxy<T>&& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_inner_prod(x.size(), x.storage(), x.stride(), y.storage(), y.stride(), z);
}
#endif

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
void
inner_prod_nested(blas::BlasVector<T, U, Tag> const& x,
		  blas::BlasVector<T, V, Tag> const& y,
		  typename Tag::template async_proxy<T>&& z, Nested&& Nest)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_inner_prod_nested(x.size(), x.storage(), x.stride(), y.storage(), y.stride(), z,
			    std::forward<Nested>(Nest));
}

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
auto
inner_prod_nested(blas::BlasVector<T, U, Tag> const& x,
		  blas::BlasVector<T, V, Tag> const& y,
		  Nested&& Nest)
{
   using result_value = remove_proxy_t<decltype(Nested(std::declval<T>(),std::declval<T>()))>;
   auto z = Tag::template allocate_async_ref<result_value>();
   vector_inner_prod(x,y,z, std::forward<Nested>(Nest));
   return z;
}

#if 0
template <typename T, typename U, typename V, typename Tag>
inline
void
add_inner_prod(blas::BlasVector<T, U, Tag> const& x,
	       blas::BlasVector<T, V, Tag> const& y,
	       typename Tag::template async_ref<T>& z)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_add_inner_prod_nested(x.size(), x.storage(), y.storage(), z, InnerProd());
}
#endif

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
typename Tag::template async_ref<T>
add_inner_prod_nested(blas::BlasVector<T, U, Tag> const& x,
		      blas::BlasVector<T, V, Tag> const& y, Nested&& Nest)
{
   typename Tag::template async_ref<T> z;
   vector_add_inner_prod_nested(x,y,z,std::forward<Nested>(Nest));
   return z;
}

template <typename T, typename U, typename V, typename Tag>
void assign_permutation(BlasVector<T, U, Tag>& A, BlasVector<T, V, Tag> const& B, int const* Iter)
{
   vector_permute(A.size(), B.storage(), B.stride(), A.storage(), A.stride(), Iter);
}

template <typename T, typename U, typename V, typename Tag>
void assign_permutation(BlasVectorProxy<T, U, Tag>&& A, BlasVector<T, V, Tag> const& B, int const* Iter)
{
   vector_permute(A.size(), B.storage(), B.stride(), std::move(A).storage(), A.stride(), Iter);
}

//
// MATRIX high-level BLAS wrappers, that forward from a matrix/vector ref to low-level storage
//

// assign_slice is the catch-all for assigning combinations of range, slice, and indices.  Only a few combinations
// are supported; add others as needed

template <typename T, typename U, typename V, typename Tag>
void
assign_slice(NormalMatrix<T, U, Tag>& Out, NormalMatrix<T, V, Tag> const& In,
             std::vector<int> const& RowTrans, Range ColTrans);

template <typename T, typename U, typename V, typename Tag>
inline
void
assign_slice(NormalMatrix<T, U, Tag>& Out, NormalMatrix<T, V, Tag> const& In,
             std::vector<int> const& RowTrans, Range ColTrans)
{
   CHECK_EQUAL(Out.rows(), RowTrans.size());
   CHECK_EQUAL(Out.cols(), ColTrans.size());
   for (int r = 0; r < RowTrans.size(); ++r)
   {
      Out.as_derived().row(r) = In.as_derived().row(RowTrans[r])[ColTrans];
   }
}

template <typename T, typename U, typename V, typename Tag>
inline
void
assign_slice(NormalMatrix<T, U, Tag>& Out, NormalMatrix<T, V, Tag> const& In,
             Range RowTrans, std::vector<int> const& ColTrans)
{
   CHECK_EQUAL(Out.rows(), RowTrans.size());
   CHECK_EQUAL(Out.cols(), ColTrans.size());
   for (int c = 0; c < ColTrans.size(); ++c)
   {
      Out.as_derived().col(c) = In.as_derived().col(ColTrans[c])[RowTrans];
   }
}

template <typename T, typename U, typename Tag>
inline
void clear(NormalMatrix<T, U, Tag>& C)
{
   matrix_clear(C.rows(), C.cols(), C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename Tag>
inline
void clear(NormalMatrixProxy<T, U, Tag>&& C)
{
   matrix_clear(C.rows(), C.cols(), std::move(C).storage(), C.leading_dimension());
}

template <typename T, typename U, typename Tag>
inline
void matrix_fill(T const& x, NormalMatrix<T, U, Tag>& C)
{
   matrix_fill(x, C.rows(), C.cols(), C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename Tag>
inline
void matrix_fill(T const& x, NormalMatrixProxy<T, U, Tag>&& C)
{
   matrix_fill(x, C.rows(), C.cols(), std::move(C).storage(), C.leading_dimension());
}

template <typename T, typename U, typename Tag>
inline
void scale(NormalMatrix<T, U, Tag>& C, T x)
{
   matrix_scale(C.rows(), C.cols(), x, C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename Tag>
inline
void
norm_frob_sq(blas::BlasMatrix<T, U, Tag> const& x,
	     typename Tag::template async_ref<decltype(norm_frob(std::declval<real_t<T>>()))>& z)
{
   matrix_norm_frob_sq(x.rows(), x.cols(), x.storage(), x.leading_dimension(), z);
}

template <typename T, typename U, typename Tag>
//decltype(norm_frob_sq(std::declval<T>()))
auto
   norm_frob_sq(blas::BlasMatrix<T, U, Tag> const& x)
{
   using blas::norm_frob_sq;
   typename Tag::template async_ref<decltype(norm_frob_sq(std::declval<T>()))> z;
   norm_frob_sq(x, z);
   return get_wait(z);
}

#if 0
// TODO: how to implement sqrt on the GPU?
template <typename T, typename U, typename Tag>
inline
void
norm_frob_sq(blas::BlasMatrix<T, U, Tag> const& x,
	  typename Tag::template async_ref<decltype(norm_frob(std::declval<T>()))>& z)
{
   matrix_norm_frob_sq(x.rows(), x.cols(), x.storage(), x.leading_dimension(), z);
}
#endif

template <typename T, typename U, typename Tag>
auto
norm_frob(blas::BlasMatrix<T, U, Tag> const& x)
{
   using std::sqrt;
   return sqrt(norm_frob_sq(x));
}

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
void
inner_prod_nested(blas::BlasMatrix<T, U, Tag> const& x,
		  blas::BlasMatrix<T, V, Tag> const& y,
		  typename Tag::template async_ref<T>& z, Nested&& Nest)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   matrix_inner_prod_nested(x.trans(), y.trans(), x.rows(), x.cols(),
			    x.storage(), x.leading_dimension(),
			    y.storage(), y.leading_dimension(),
			    z, std::forward<Nested>(Nest));
}

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
void
inner_prod_nested(blas::BlasMatrix<T, U, Tag> const& x,
		  blas::BlasMatrix<T, V, Tag> const& y,
		  typename Tag::template async_proxy<T>&& z, Nested&& Nest)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   matrix_inner_prod_nested(x.trans(), y.trans(), x.rows(), x.cols(),
			    x.storage(), x.leading_dimension(),
			    y.storage(), y.leading_dimension(),
			    z, std::forward<Nested>(Nest));
}

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
typename Tag::template async_ref<T>
inner_prod_nested(blas::BlasMatrix<T, U, Tag> const& x,
		  blas::BlasMatrix<T, V, Tag> const& y,
		  Nested&& Nest);
// not implemented, but used for return-type deduction

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
void
add_inner_prod_nested(blas::BlasMatrix<T, U, Tag> const& x,
		      blas::BlasMatrix<T, V, Tag> const& y,
		      typename Tag::template async_ref<T>& z,
		      Nested&& Nest)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   matrix_add_inner_prod_nested(x.trans(), y.trans(), x.rows(), x.cols(),
				x.storage(), x.leading_dimension(),
				y.storage(), y.leading_dimension(),
				z, std::forward<Nested>(Nest));
}

template <typename T, typename U, typename V, typename Tag, typename Nested>
inline
void
add_inner_prod_nested(blas::BlasMatrix<T, U, Tag> const& x,
		      blas::BlasMatrix<T, V, Tag> const& y,
		      typename Tag::template async_ref<T>&& z,
		      Nested&& Nest)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   matrix_add_inner_prod_nested(x.trans(), y.trans(), x.rows(), x.cols(),
				x.storage(), x.leading_dimension(),
				y.storage(), y.leading_dimension(),
				z, std::forward<Nested>(Nest));
}

// middle-level wrappers

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void gemv(T alpha, BlasMatrix<T, U, Tag> const& A, BlasVector<T, V, Tag> const& x, T beta,
          BlasVector<T, W, Tag>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(A.rows(), y.size());
   gemv(A.trans(), A.rows(), A.cols(), alpha, A.storage(),
	A.leading_dimension(), x.storage(), x.stride(),
        beta, y.storage(), y.stride());
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void gemv(T alpha, BlasMatrix<T, U, Tag> const& A, BlasVector<T, V, Tag> const& x, T beta,
          BlasVectorProxy<T, W, Tag>&& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(A.rows(), y.size());
   gemv(A.trans(), A.rows(), A.cols(), alpha, A.storage(),
	A.leading_dimension(), x.storage(), x.stride(),
        beta, std::move(y).storage(), y.stride());
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void gemm(T alpha, BlasMatrix<T, U, Tag> const& A,
          T beta, BlasMatrix<T, V, Tag> const& B,
          NormalMatrix<T, W, Tag>& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   gemm(A.trans(), B.trans(), A.rows(), A.cols(), B.cols(), alpha, A.storage(),
        A.leading_dimension(), B.storage(), B.leading_dimension(), beta,
        C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void gemm(T alpha, BlasMatrix<T, U, Tag> const& A,
          T beta, BlasMatrix<T, V, Tag> const& B,
          NormalMatrixProxy<T, W, Tag>&& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   gemm(A.trans(), B.trans(), A.rows(), A.cols(), B.cols(), alpha, A.storage(),
        A.leading_dimension(), B.storage(), B.leading_dimension(), beta,
        std::move(C).storage(), C.leading_dimension());
}

// for a multiplication involving a ScaledMatrix, push the factors into the
// alpha parameter

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void gemm(T alpha, ScaledMatrix<T, U, Tag> const& A,
          T beta, MatrixRef<T, V, Tag> const& B, W&& C)
{
   gemm(alpha*A.factor(), A.base(), beta, B.as_derived(), std::forward<W>(C).as_derived());
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void gemm(T alpha, MatrixRef<T, U, Tag> const& A,
          T beta, ScaledMatrix<T, V, Tag> const& B, W&& C)
{
   gemm(alpha*B.factor(), A.as_derived(), beta, B.base(), std::forward<W>(C).as_derived());
}

template <typename T, typename U, typename V, typename W, typename Tag>
inline
void gemm(T alpha, ScaledMatrix<T, U, Tag> const& A,
          T beta, ScaledMatrix<T, V, Tag> const& B, W&& C)
{
   gemm(alpha*A.factor()*B.factor(), A.base(), beta, B.base(), std::forward<W>(C).as_derived());
}

// gemm for DiagonalMatrix * Matrix

template <typename Scalar, typename T, typename U, typename V, typename W, typename X, typename Y, typename Tag>
inline
void gemm(Scalar alpha, DiagonalBlasMatrix<T, U, Tag> const& A,
          Scalar beta, NormalMatrix<V, W, Tag> const& B,
          NormalMatrix<X, Y, Tag>& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   
   // The dgmm call only supports beta=0, alpha=1 so other cases need emulation
   if (beta != number_traits<Scalar>::zero())
   {
      PANIC("unsupported");
#if 0
      Matrix<X, Tag> Temp(A*B);
      if (alpha != number_traits<Scalar>::identity())
	 Temp *= alpha;
      C += Temp;
#endif
   }
   else
   {
      dgmm(A.rows(), B.cols(), A.storage(), A.stride(), 
	   B.storage(), B.leading_dimension(),
	   C.storage(), C.leading_dimension());
      if (alpha != number_traits<Scalar>::identity())
	 C.as_derived() *= alpha;
   }
}

template <typename Scalar, typename T, typename U, typename V, typename W, typename X, typename Y, typename Tag>
inline
void gemm(Scalar alpha, DiagonalBlasMatrix<T, U, Tag> const& A,
          Scalar beta, NormalMatrix<V, W, Tag> const& B,
          NormalMatrixProxy<X, Y, Tag>&& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   
   // The dgmm call only supports beta=0, alpha=1 so other cases need emulation
   if (beta != number_traits<Scalar>::zero())
   {
      PANIC("unsupported");
#if 0
      // need a temporary
      Matrix<X, Tag> Temp(A*B);
      if (alpha != number_traits<Scalar>::identity())
	 Temp *= alpha;
      std::move(C) += Temp;
#endif
   }
   else
   {
      dgmm(A.rows(), B.cols(), A.storage(), A.stride(), 
	   B.storage(), B.leading_dimension(),
	   std::move(C).storage(), C.leading_dimension());
      if (alpha != number_traits<Scalar>::identity())
	 std::move(C).as_derived() *= alpha;
   }
}

// gemm for Matrix * DiagonalMatrix

template <typename Scalar, typename T, typename U, typename V, typename W, typename X, typename Y, typename Tag>
inline
void gemm(Scalar alpha, BlasMatrix<T, U, Tag> const& A,
          Scalar beta, DiagonalBlasMatrix<V, W, Tag> const& B,
          NormalMatrix<X, Y, Tag>& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   // The dgmm call only supports beta=0, alpha=1 so other cases need emulation
   if (beta != number_traits<Scalar>::zero())
   {
      PANIC("unsupported");
   }
   else
   {
      gdmm(A.rows(), A.cols(), A.storage(), A.leading_dimension(),
	   B.storage(), B.stride(),
	   C.storage(), C.leading_dimension());
      if (alpha != number_traits<Scalar>::identity())
	 C.as_derived() *= alpha;
   }
}

template <typename Scalar, typename T, typename U, typename V, typename W, typename X, typename Y, typename Tag>
inline
void gemm(Scalar alpha, BlasMatrix<T, U, Tag> const& A,
          Scalar beta, DiagonalBlasMatrix<V, W, Tag> const& B,
          NormalMatrixProxy<X, Y, Tag>&& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   // The dgmm call only supports beta=0, alpha=1 so other cases need emulation
   if (beta != number_traits<Scalar>::zero())
   {
      PANIC("unsupported");
   }
   else
   {
      gdmm(A.rows(), A.cols(), A.storage(), A.leading_dimension(),
	   B.storage(), B.stride(),
	   std::move(C).storage(), C.leading_dimension());
      if (alpha != number_traits<Scalar>::identity())
	 std::move(C).as_derived() *= alpha;
   }
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_copy_scaled(T alpha, BlasMatrix<T, U, Tag> const& A, NormalMatrix<T, V, Tag>& C)
{
   matrix_copy_scaled(A.trans(), A.rows(), A.cols(), alpha, A.storage(), A.leading_dimension(),
                      C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_copy_scaled(T alpha, BlasMatrix<T, U, Tag> const& A, NormalMatrixProxy<T, V, Tag>&& C)
{
   matrix_copy_scaled(A.trans(), A.rows(), A.cols(), alpha, A.storage(), A.leading_dimension(),
                      std::move(C).storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_copy(BlasMatrix<T, U, Tag> const& A, NormalMatrix<T, V, Tag>& C)
{
   matrix_copy(A.trans(), A.rows(), A.cols(), A.storage(), A.leading_dimension(),
	       C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_copy(BlasMatrix<T, U, Tag> const& A, NormalMatrixProxy<T, V, Tag>&& C)
{
   matrix_copy(A.trans(), A.rows(), A.cols(), A.storage(), A.leading_dimension(),
               std::move(C).storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_add_scaled(T alpha, BlasMatrix<T, U, Tag> const& A, NormalMatrix<T, V, Tag>& C)
{
   matrix_add_scaled(A.trans(), A.rows(), A.cols(), alpha, A.storage(), A.leading_dimension(),
                     C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_add_scaled(T alpha, BlasMatrix<T, U, Tag> const& A, NormalMatrixProxy<T, V, Tag>&& C)
{
   matrix_add_scaled(A.trans(), A.rows(), A.cols(), alpha, A.storage(), A.leading_dimension(),
                     std::move(C).storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_add(BlasMatrix<T, U, Tag> const& A, NormalMatrix<T, V, Tag>& C)
{
   matrix_add(A.trans(), A.rows(), A.cols(), A.storage(), A.leading_dimension(),
              C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename V, typename Tag>
inline
void matrix_add(BlasMatrix<T, U, Tag> const& A, NormalMatrixProxy<T, V, Tag>&& C)
{
   matrix_add(A.trans(), A.rows(), A.cols(), A.storage(), A.leading_dimension(),
              C.storage(), C.leading_dimension());
}

//
// middle-layer ARPACK wrappers, that forward from a matrix/vector ref to low-level storage
//

//
// Diagonalize
//
// Diagonalizes a matrix.  The left and right eigenvectors are returned
// as successive column vectors of Left and Right respectively.
// The input matrix is destroyed, so we have two versions were we either copy
// the input matrix, or if it is an rvalue-reference, we can destroy it.
/// On input, M is an NxN matrix, V is an N-dimensional vector,
// and Left and Right are NxN matrices.
//
// TODO: add NormalVectorProxy versions

template <typename U, typename V, typename W, typename X, typename Tag>
inline
void Diagonalize(NormalMatrix<std::complex<double>, U, Tag>&& M,
		 NormalVector<std::complex<double>, V, Tag>& v,
		 NormalMatrix<std::complex<double>, W, Tag>& Left,
		 NormalMatrix<std::complex<double>, X, Tag>& Right)
{
   CHECK_EQUAL(M.rows(), M.cols());
   CHECK_EQUAL(M.rows(), v.size());
   CHECK_EQUAL(M.rows(), Left.rows());
   CHECK_EQUAL(M.rows(), Left.cols());
   CHECK_EQUAL(M.rows(), Right.rows());
   CHECK_EQUAL(M.rows(), Right.cols());
   Diagonalize(M.rows(), std::move(M).storage(), M.leading_dimension(), v.storage(),
	       Left.storage(), Left.leading_dimension(),
	       Right.storage(), Right.leading_dimension());
}

template <typename U, typename V, typename W, typename X, typename Tag>
inline
void Diagonalize(NormalMatrix<std::complex<double>, U, Tag> const& M,
		 NormalVector<std::complex<double>, V, Tag>& v,
		 NormalMatrix<std::complex<double>, W, Tag>& Left,
		 NormalMatrix<std::complex<double>, X, Tag>& Right)
{
   Diagonalize(copy(M), v, Left, Right);
}

//
// DiagonalizeSymmetric
//
// Diagonalizes a matrix in-place.  The matrix is replaced by the transform matrix, with
// the eigenvectors as sucessive column-vectors.  For input matrix M,
// X = M' is the transform matrix, E is the diagonal matrix of eigenvalues,
// we have MX = XE
//

template <typename U, typename V, typename Tag>
inline
void DiagonalizeSymmetric(NormalMatrix<double, U, Tag>& M, NormalVector<double, V, Tag>& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   DiagonalizeSymmetric(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

// Version that takes a proxy-reference for the eigenvalues
template <typename U, typename V, typename Tag>
inline
void DiagonalizeSymmetric(NormalMatrix<double, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   DiagonalizeSymmetric(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

// TODO: we could also add versions where M is passed as a NormalMatrixProxy

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<double, U, Tag>& M, NormalVector<double, V, Tag>& v)
{
   DiagonalizeSymmetric(M.as_derived(), v.as_derived());
}

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<double, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   DiagonalizeSymmetric(M.as_derived(), static_cast<V&&>(v.as_derived()));
}

//
// complex
//

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<std::complex<double>, U, Tag>& M, NormalVector<double, V, Tag>& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   DiagonalizeHermitian(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<std::complex<double>, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   DiagonalizeHermitian(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

//
// SingularValueDecomposition
//

template <typename M, typename U, typename D, typename V, typename Tag>
void
SingularValueDecomposition(NormalMatrix<double, M, Tag>&& Mmat, NormalMatrix<double, U, Tag>& Umat,
                           NormalVector<double, D, Tag>& Dvec, NormalMatrix<double, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), std::min(Mmat.rows(), Mmat.cols()));
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.rows(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
			      Umat.storage(), Umat.leading_domension(),
                              Vmat.storage(), Vmat.leading_dimension());
}

template <typename M, typename U, typename D, typename V, typename Tag>
void
SingularValueDecomposition(NormalMatrix<std::complex<double>, M, Tag>&& Mmat,
                           NormalMatrix<std::complex<double>, U, Tag>& Umat,
                           NormalVector<double, D, Tag>& Dvec,
                           NormalMatrix<std::complex<double>, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), std::min(Mmat.rows(), Mmat.cols()));
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.rows(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
			      Umat.storage(), Umat.leading_dimension(),
                              Vmat.storage(), Vmat.leading_dimension());
}

} // namespace blas

#include "matrixref.icc"

#endif
