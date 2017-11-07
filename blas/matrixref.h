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

template <typename ValueType, typename DerivedType, typename Tag>
class MatrixRef
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;

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
};

// derived class for a diagonal matrix

template <typename ValueType, typename DerivedType, typename Tag>
class DiagonalMatrixRef : MatrixRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;
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
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

      BlasMatrix() = default;
      ~BlasMatrix() = default;
      BlasMatrix(BlasMatrix&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }

      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation)
      char trans() const { return this->as_derived().trans(); }

      storage_type storage() & { return this->as_derived().storage(); }
      const_storage_type storage() const& { return this->as_derived().storage(); }
};

// Denotes an L-value matrix, which is a BlasMatrix with the trans() parameter equal to 'N'
template <typename ValueType, typename DerivedType, typename Tag>
class NormalMatrix : public BlasMatrix<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

      NormalMatrix() = default;
      ~NormalMatrix() = default;
      NormalMatrix(NormalMatrix&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }

      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation)
      constexpr char trans() const { return 'N'; }

      storage_type storage() & { return this->as_derived().storage(); }
      const_storage_type storage() const& { return this->as_derived().storage(); }
};

// proxy class for a matrix that can appear on the left-hand side of an expression

template <typename ValueType, typename DerivedType, typename Tag>
class MatrixProxy : public BlasMatrix<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

      MatrixProxy() = default;
      ~MatrixProxy() = default;
      MatrixProxy(MatrixProxy&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }

      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation)
      char trans() const { return this->as_derived().trans(); }

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
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

      BlasMatrixTrans(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.cols(); }
      int cols() const { return Base.rows(); }
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
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

      BlasMatrixHerm(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.cols(); }
      int cols() const { return Base.rows(); }
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
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

      BlasMatrixConj(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.rows(); }
      int cols() const { return Base.cols(); }
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

// assignment

template <typename T, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_copy(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add_scaled(-number_traits<T>::identity(), B.as_derived(), A.as_derived());
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
void subtract(MatrixRef<T, U, Tag>& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(-B.factor(), B.base(), A.as_derived());
}

// expression template for alpha * op(A) * op(B)

template <typename T, typename U, typename V, typename Tag>
struct MatrixProduct : public MatrixRef<T, MatrixProduct<T, U, V, Tag>, Tag>
{
   MatrixProduct(MatrixRef<T, U, Tag> const& A_, MatrixRef<T, V, Tag> const& B_)
      : Factor(number_traits<T>::identity()), A(A_.as_derived()), B(B_.as_derived()) {}

   MatrixProduct(T const& Factor_, MatrixRef<T, U, Tag> const& A_, MatrixRef<T, V, Tag> const& B_)
      : Factor(Factor_), A(A_.as_derived()), B(B_.as_derived()) {}

   int rows() const { return A.rows(); }
   int cols() const { return B.cols(); }
   T factor() const { return Factor; }

   T Factor;
   U const& A;
   V const& B;
};

template <typename T, typename U, typename V, typename Tag>
typename blas_traits<Tag>::template matrix_type<T>
evaluate(MatrixProduct<T, U, V, Tag>&& x)
{
   typename blas_traits<Tag>::template matrix_type<T> Result(x.rows(), x.cols());
   assign(Result, std::move(x));
   return Result;
}

template <typename T, typename U, typename V, typename Tag>
inline
MatrixProduct<T, U, V, Tag>
operator*(MatrixRef<T, U, Tag> const& A, MatrixRef<T, V, Tag> const& B)
{
   return MatrixProduct<T, U, V, Tag>(A.as_derived(), B.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
MatrixProduct<T, U, V, Tag>
operator*(ScaledMatrix<T, U, Tag> const& A, MatrixRef<T, V, Tag> const& B)
{
   return MatrixProduct<T, U, V, Tag>(A.factor(), A.base(), B.as_derived());
}

// matrix product with a scalar, alpha * A * B

template <typename T, typename Derived, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, Derived, Tag>& C, MatrixProduct<T, U, V, Tag> const& a)
{
   gemm(a.Factor, a.A, number_traits<T>::zero(), a.B, C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, Derived, Tag>& C, MatrixProduct<T, U, V, Tag> const& a)
{
   gemm(a.Factor, a.A, number_traits<T>::identity(), a.B, C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, Derived, Tag>& C, MatrixProduct<T, U, V, Tag> const& a)
{
   gemm(-a.Factor, a.A, number_traits<T>::zero(), a.B, C.as_derived());
}

// miscellaneous

template <typename T, typename U, typename Tag>
void
inline
trace(BlasMatrix<T, U, Tag> const& x, typename blas_traits<Tag>::template async_proxy<T>&& r)
{
   vector_sum(x.as_derived().diagonal(), std::move(r));
}

template <typename T, typename U, typename Tag>
void
inline
trace(BlasMatrix<T, U, Tag> const& x, typename blas_traits<Tag>::template async_ref<T>& r)
{
   vector_sum(x.as_derived().diagonal(), r);
}

template <typename T, typename U, typename Tag>
inline
T
trace(BlasMatrix<T, U, Tag> const& x)
{
   return vector_sum(x.as_derived().diagonal());
}

//
// DiagonalBlasMatrix - a matrix that physically exists in memory with a stride
// Essentially a different view of a BlasVector.
//

template <typename ValueType, typename DerivedType, typename Tag>
class DiagonalBlasMatrix : DiagonalMatrixRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

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
operator*(VectorRef<T, U, Tag> const& A, ScaledVector<T, V, Tag> const& B)
{
   return MatrixVectorProduct<T, U, V, Tag>(B.factor(), A.as_derived(), B.base());
}

template <typename T, typename U, typename V, typename Tag>
MatrixVectorProduct<T, U, V, Tag>
operator*(ScaledMatrix<T, U, Tag> const& A, ScaledVector<T, V, Tag> const& B)
{
   return MatrixVectorProduct<T, U, V, Tag>(A.factor() * B.factor(), A.base(), B.as_derived());
}

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

template <typename T, typename U, typename V, typename Tag>
inline
void vector_copy(blas::BlasVector<T, U, Tag> const& x, BlasVector<T, V, Tag>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   vector_copy(x.size(), x.storage(), x.stride(), y.storage(), y.stride());
}

template <typename T, typename U, typename V, typename Tag>
inline
void vector_copy(blas::BlasVector<T, U, Tag> const& x, BlasVectorProxy<T, V, Tag>&& y)
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

template <typename T, typename U, typename Tag>
inline
void vector_scale(T alpha, BlasVectorProxy<T, U, Tag>&& y)
{
   vector_scale(y.size(), alpha, std::move(y).storage(), y.stride());
}

template <typename T, typename U, typename Tag>
void
vector_sum(blas::BlasVector<T, U, Tag> const& x, typename blas_traits<Tag>::template async_proxy<T>&& y)
{
   vector_sum(x.size(), x.storage(), x.stride(), y);
}

template <typename T, typename U, typename Tag>
void
vector_sum(blas::BlasVector<T, U, Tag> const& x, typename blas_traits<Tag>::template async_ref<T>& y)
{
   vector_sum(x.size(), x.storage(), x.stride(), y);
}

template <typename T, typename U, typename Tag>
T
vector_sum(blas::BlasVector<T, U, Tag> const& x)
{
   typename blas_traits<Tag>::template async_ref<T> y;
   vector_sum(x.size(), x.storage(), x.stride(), y);
   return y;
}

//
// MATRIX middle-layer BLAS wrappers, that forward from a matrix/vector ref to low-level storage
//

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
          MatrixProxy<T, W, Tag>&& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   gemm(A.trans(), B.trans(), A.rows(), A.cols(), B.cols(), alpha, A.storage(),
        A.leading_dimension(), B.storage(), B.leading_dimension(), beta,
        std::move(C).storage(), C.leading_dimension());
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
void matrix_copy_scaled(T alpha, BlasMatrix<T, U, Tag> const& A, MatrixProxy<T, V, Tag>&& C)
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
void matrix_copy(BlasMatrix<T, U, Tag> const& A, MatrixProxy<T, V, Tag>&& C)
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
void matrix_add_scaled(T alpha, BlasMatrix<T, U, Tag> const& A, MatrixProxy<T, V, Tag>&& C)
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
void matrix_add(BlasMatrix<T, U, Tag> const& A, MatrixProxy<T, V, Tag>&& C)
{
   matrix_add(A.trans(), A.rows(), A.cols(), A.storage(), A.leading_dimension(),
              C.storage(), C.leading_dimension());
}

//
// middle-layer ARPACK wrappers, that forward from a matrix/vector ref to low-level storage
//

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
   DiagonalizeSymmetric(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

// Version that takes a proxy-reference for the eigenvalues
template <typename U, typename V, typename Tag>
inline
void DiagonalizeSymmetric(NormalMatrix<double, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   DiagonalizeSymmetric(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

// TODO: we could also add versions where M is passed as a MatrixProxy

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
   DiagonalizeHermitian(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<std::complex<double>, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   DiagonalizeHermitian(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

} // namespace blas

#endif
