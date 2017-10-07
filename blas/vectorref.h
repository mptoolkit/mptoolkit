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

#if !defined(MPTOOLKIT_BLAS_VECTORREF_H)
#define MPTOOLKIT_BLAS_VECTORREF_H

#include "matrixref.h"
#include "safe-conversions.h"
#include "number_traits.h"

namespace blas
{

//
// VectorRef : generic base class for a vector using expression templates.
// VectorRef is a reference to a matrix of type ValueType.  BaseType is
// a concrete type, eg the two currently supported types are
// Vector<ValueType> or gpu_matrix<ValueType>.
//
// If DerivedType == BaseType, then the matrix is actually a concrete matrix
// and can be constructed, used as a temporary, etc.
//

template <typename ValueType, typename DerivedType, typename Tag>
class VectorRef
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;

      // default construction and move construction are defined, no copying.

      VectorRef() = default;
      ~VectorRef() = default;
      VectorRef(VectorRef&& Other) = default;
      VectorRef(VectorRef const&) = delete;
      VectorRef& operator=(VectorRef&&) = delete;
      VectorRef& operator=(VectorRef const&) = delete;

      derived_type&& as_derived() && { return static_cast<derived_type&&>(*static_cast<derived_type*>(this)); }

      derived_type& as_derived() & { return *static_cast<derived_type*>(this); }
      derived_type const& as_derived() const& { return *static_cast<derived_type const*>(this); }

      int size() const { return this->as_derived().size(); }
};

// assignment

template <typename T, typename U, typename V, typename Tag>
void assign(VectorRef<T, U, Tag>& A, VectorRef<T, V, Tag> const& B)
{
   vector_copy(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
void add(VectorRef<T, U, Tag>& A, VectorRef<T, V, Tag> const& B)
{
   vector_add(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
void subtract(VectorRef<T, U, Tag>& A, VectorRef<T, V, Tag> const& B)
{
   vector_add_scaled(-number_traits<T>::identity(), B.as_derived(), A.as_derived());
}

// Specialization of a VectorRef for a matrix that can be directly used in BLAS-like calls,
// ie it physically exists and can be addressed (although not necessarily in main memory,
// eg it might be on some other device such as a GPU), and has a fixed stride.
//
// The BaseType is the concrete type, which allows us to specialize blas calls for different
// devices.
// storage() returns an opaque type that contains the backing storage of the vector,
// and stride() is the array stride (the INC parameter in BLAS calls).

template <typename ValueType, typename DerivedType, typename Tag>
class BlasVector : public VectorRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

      BlasVector() = default;
      ~BlasVector() = default;
      BlasVector(BlasVector&& Other) = default;

      int stride() const { return this->as_derived().stride(); }

      storage_type storage() { return this->as_derived().storage(); }
      const_storage_type storage() const { return this->as_derived().storage(); }
};

// proxy class for the complex conjugate of a vector

template <typename T, typename BaseType, typename Tag>
class VectorConj : public VectorRef<T, VectorConj<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<base_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<base_type>::template const_storage_type<value_type>;

      VectorConj(base_type const& Base_) : Base(Base_) {}

      int size() const { return Base.size(); }

      base_type const& base() const { return Base; }

   private:
      base_type const& Base;
};

template <typename T, typename BaseType, typename Tag>
VectorConj<T, BaseType, Tag>
conj(VectorRef<T, BaseType, Tag> const& x)
{
   return VectorConj<T, BaseType, Tag>(x.as_derived());
}

// specialization for real types, conj() reduces to a no-op
template <typename BaseType, typename Tag>
BaseType const&
conj(VectorRef<float, BaseType, Tag> const& x)
{
   return x.as_derived();
}

template <typename BaseType, typename Tag>
BaseType const&
conj(VectorRef<double, BaseType, Tag> const& x)
{
   return x.as_derived();
}

#if defined(HAVE_FLOAT128)
template <typename BaseType, typename Tag>
BaseType const&
conj(VectorRef<float128, BaseType, Tag> const& x)
{
   return x.as_derived();
}
#endif

// expression template for alpha * A

template <typename T, typename BaseType, typename Tag>
class ScaledVector : public VectorRef<T, ScaledVector<T, BaseType, Tag>, Tag>
{
   public:
      using base_type = BaseType;

      ScaledVector(T const& Factor_, base_type const& A_) : Factor(Factor_), A(A_) {}

      int size() const { return A.size(); }

      base_type const& base() const { return A; }
      T factor() const { return Factor; }

   private:
      T Factor;
      base_type const& A;
};

template <typename T, typename BaseType, typename Tag, typename X>
ScaledVector<decltype(safe_convert<T>(std::declval<X>())), BaseType, Tag>
operator*(X const& alpha, VectorRef<T, BaseType, Tag> const& M)
{
   return ScaledVector<T, BaseType, Tag>(safe_convert<T>(alpha), M.as_derived());
}

// assignment involving ScaledVector

template <typename T, typename U, typename V, typename Tag>
void assign(VectorRef<T, U, Tag>& A, ScaledVector<T, V, Tag> const& B)
{
   vector_copy_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
void assign(VectorRef<T, U, Tag>&& A, ScaledVector<T, V, Tag> const& B)
{
   vector_copy_scaled(B.factor(), B.base(), static_cast<U&&>(A.as_derived()));
}

template <typename T, typename U, typename V, typename Tag>
void add(VectorRef<T, U, Tag>& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
void add(VectorRef<T, U, Tag>&& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(B.factor(), B.base(), static_cast<U&&>(A.as_derived()));
}

template <typename T, typename U, typename V, typename Tag>
void subtract(VectorRef<T, U, Tag>& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(-B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
void subtract(VectorRef<T, U, Tag>&& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(-B.factor(), B.base(), static_cast<U&&>(A.as_derived()));
}

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
void add(VectorRef<T, Derived, Tag>& C, MatrixVectorProduct<T, U, V, Tag> const& a)
{
   gemv(a.Factor, a.A, a.B, number_traits<T>::identity(), C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
void subtract(VectorRef<T, Derived, Tag>& C, MatrixVectorProduct<T, U, V, Tag> const& a)
{
   gemv(-a.Factor, a.A, a.B, number_traits<T>::identity(), C.as_derived());
}

} // namespace blas

#endif
