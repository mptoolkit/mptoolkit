// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/vectorref.h
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

#include "safe-conversions.h"
#include "number_traits.h"
#include "stride_ptr.h"

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

struct cpu_tag;

template <typename T, typename Tag = cpu_tag>
class Vector;

template <typename ValueType, typename DerivedType, typename Tag>
class VectorRef
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;
      using remove_proxy_t = Vector<ValueType, Tag>;

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

template <typename T, typename U, typename T2, typename V, typename Tag>
void assign(VectorRef<T, U, Tag>& A, VectorRef<T2, V, Tag> const& B)
{
   vector_copy(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename T2, typename V, typename Tag>
void assign(VectorRef<T, U, Tag>&& A, VectorRef<T2, V, Tag> const& B)
{
   vector_copy(B.as_derived(), std::move(A).as_derived());
}

template <typename T, typename U, typename T2, typename V, typename Tag>
void assign_copy(VectorRef<T, U, Tag>& A, VectorRef<T2, V, Tag> const& B)
{
   vector_deep_copy(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename T2, typename V, typename Tag>
void assign_copy(VectorRef<T, U, Tag>&& A, VectorRef<T2, V, Tag> const& B)
{
   vector_deep_copy(B.as_derived(), std::move(A).as_derived());
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

template <typename T, typename U, typename Tag>
void fill(VectorRef<T, U, Tag>& A, T const& x)
{
   vector_fill(x, A.as_derived());
}

// Specialization of a VectorRef for a matrix that can be directly used in BLAS-like calls,
// ie it physically exists and can be addressed (although not necessarily in main memory,
// eg it might be on some other device such as a GPU), and has a fixed stride.
//
// The BaseType is the concrete type, which allows us to specialize blas calls for different
// devices.
// storage() returns an opaque type that contains the backing storage of the vector,
// and stride() is the array stride (the INC parameter in BLAS calls).

// some forward declarations for the operator[]
template <typename T, typename Tag>
class vector_view;
template <typename T, typename Tag>
class const_vector_view;
class Range;

template <typename ValueType, typename DerivedType, typename Tag>
class BlasVector : public VectorRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      BlasVector() = default;
      ~BlasVector() = default;
      BlasVector(BlasVector&& Other) = default;

      vector_view<ValueType, Tag> operator[](Range r) &;

      const_vector_view<ValueType, Tag> operator[](Range r) const&;

      int stride() const { return this->as_derived().stride(); }

      storage_type storage() & { return this->as_derived().storage(); }
      storage_type storage() && { return this->as_derived().storage(); }
      const_storage_type storage() const & { return this->as_derived().storage(); }
};

template <typename T, typename U>
stride_ptr<T> begin(BlasVector<T,U,cpu_tag>& v)
{
   return stride_ptr<T>(v.storage(), v.stride());
}

template <typename T, typename U>
stride_ptr<T const> begin(BlasVector<T,U,cpu_tag> const& v)
{
   return stride_ptr<T const>(v.storage(), v.stride());
}

template <typename T, typename U>
stride_ptr<T const> cbegin(BlasVector<T,U,cpu_tag> const& v)
{
   return stride_ptr<T const>(v.storage(), v.stride());
}

template <typename T, typename U>
stride_ptr<T> end(BlasVector<T,U,cpu_tag>& v)
{
   return stride_ptr<T>(v.storage()+v.stride()*v.size(), v.stride());
}

template <typename T, typename U>
stride_ptr<T const> end(BlasVector<T,U,cpu_tag> const& v)
{
   return stride_ptr<T const>(v.storage()+v.stride()*v.size(), v.stride());
}

template <typename T, typename U>
stride_ptr<T const> cend(BlasVector<T,U,cpu_tag> const& v)
{
   return stride_ptr<T const>(v.storage()+v.stride()*v.size(), v.stride());
}

// Specialization of BlasVector that can act as an rvalue-reference, ie
// a temporary proxy that can appear on the left-hand side of an expression.
template <typename ValueType, typename DerivedType, typename Tag>
class BlasVectorProxy : public BlasVector<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      BlasVectorProxy() = default;
      ~BlasVectorProxy() = default;
      BlasVectorProxy(BlasVectorProxy&& Other) = default;

      int stride() const { return this->as_derived().stride(); }

      vector_view<ValueType, Tag> operator[](Range r) &&;
      const_vector_view<ValueType, Tag> operator[](Range r) const&;

      storage_type storage() && { return static_cast<derived_type&&>(this->as_derived()).storage(); }
      const_storage_type storage() const& { return this->as_derived().storage(); }
};

template <typename T, typename U, typename V, typename Tag>
void assign(BlasVectorProxy<T, U, Tag>&& A, VectorRef<T, V, Tag> const& B)
{
   vector_copy(B.as_derived(), static_cast<U&&>(A.as_derived()));
}


template <typename T, typename U>
stride_ptr<T> begin(BlasVector<T,U,cpu_tag>&& v)
{
   return stride_ptr<T>(std::move(v).storage(), v.stride());
}

template <typename T, typename U>
stride_ptr<T> end(BlasVector<T,U,cpu_tag>&& v)
{
   return stride_ptr<T>(std::move(v).storage()+v.stride()*v.size(), v.stride());
}

// Specialization of a BlasVector for a vector that exists in memory and has stride 1

template <typename ValueType, typename DerivedType, typename Tag>
class NormalVector : public BlasVector<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      NormalVector() = default;
      ~NormalVector() = default;
      NormalVector(NormalVector&& Other) = default;

      constexpr int stride() const { return 1; }

      storage_type storage() & { return this->as_derived().storage(); }
      const_storage_type storage() const & { return this->as_derived().storage(); }
};

template <typename T, typename U>
T* begin(NormalVector<T,U,cpu_tag>& v)
{
   return v.storage();
}

template <typename T, typename U>
T const* begin(NormalVector<T,U,cpu_tag> const& v)
{
   return v.storage();
}

template <typename T, typename U>
T const* cbegin(NormalVector<T,U,cpu_tag> const& v)
{
   return v.storage();
}

template <typename T, typename U>
T* end(NormalVector<T,U,cpu_tag>& v)
{
   return v.storage()+v.size();
}

template <typename T, typename U>
T const* end(NormalVector<T,U,cpu_tag> const& v)
{
   return v.storage()+v.size();
}

template <typename T, typename U>
T const* cend(NormalVector<T,U,cpu_tag> const& v)
{
   return v.storage()+v.size();
}

// Specialization of BlasVectorProxy that can act as an rvalue-reference, ie
// a temporary proxy that can appear on the left-hand side of an expression.
template <typename ValueType, typename DerivedType, typename Tag>
class NormalVectorProxy : public BlasVectorProxy<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      NormalVectorProxy() = default;
      ~NormalVectorProxy() = default;
      NormalVectorProxy(NormalVectorProxy&& Other) = default;

      constexpr int stride() const { return 1; }

      storage_type storage() && { return static_cast<derived_type&&>(this->as_derived()).storage(); }
      const_storage_type storage() const& { return this->as_derived().storage(); }
};

template <typename T, typename U>
T* begin(NormalVectorProxy<T,U,cpu_tag>&& v)
{
   return std::move(v).storage();
}

template <typename T, typename U>
T* begin(NormalVectorProxy<T,U,cpu_tag>& v)
{
   return std::move(v).storage();
}

template <typename T, typename U>
T const* begin(NormalVectorProxy<T,U,cpu_tag> const& v)
{
   return v.storage();
}

template <typename T, typename U>
T const* cbegin(NormalVectorProxy<T,U,cpu_tag> const& v)
{
   return v.storage();
}

template <typename T, typename U>
T* end(NormalVectorProxy<T,U,cpu_tag>&& v)
{
   return std::move(v).storage()+v.size();
}

template <typename T, typename U>
T* end(NormalVectorProxy<T,U,cpu_tag>& v)
{
   return std::move(v).storage()+v.size();
}

template <typename T, typename U>
T const* end(NormalVectorProxy<T,U,cpu_tag> const& v)
{
   return v.storage()+v.size();
}

template <typename T, typename U>
T const* cend(NormalVectorProxy<T,U,cpu_tag> const& v)
{
   return v.storage()+v.size();
}

// proxy class for the complex conjugate of a vector

template <typename T, typename BaseType, typename Tag>
class VectorConj : public VectorRef<T, VectorConj<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

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
void assign(BlasVectorProxy<T, U, Tag>&& A, ScaledVector<T, V, Tag> const& B)
{
   vector_copy_scaled(B.factor(), B.base(), static_cast<U&&>(A.as_derived()));
}

template <typename T, typename U, typename V, typename Tag>
void add(VectorRef<T, U, Tag>& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
void add(BlasVectorProxy<T, U, Tag>&& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(B.factor(), B.base(), static_cast<U&&>(A.as_derived()));
}

template <typename T, typename U, typename V, typename Tag>
void subtract(VectorRef<T, U, Tag>& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(-B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
void subtract(BlasVectorProxy<T, U, Tag>&& A, ScaledVector<T, V, Tag> const& B)
{
   vector_add_scaled(-B.factor(), B.base(), std::move(A).as_derived());
		     //static_cast<U&&>(A.as_derived()));
}

template <typename T, typename U, typename Tag>
void fill(BlasVectorProxy<T, U, Tag>&& A, T const& x)
{
   vector_fill(x, std::move(A).as_derived());
}

} // namespace blas

#include "vectorref.icc"

#endif
