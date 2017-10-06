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

template <typename ValueType, typename BaseType, typename DerivedType = BaseType>
class VectorRef
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using base_type      = BaseType;

      // default construction and move construction are defined, no copying.

      VectorRef() = default;
      ~VectorRef() = default;
      VectorRef(VectorRef&& Other) = default;
      VectorRef(VectorRef const&) = delete;
      VectorRef& operator=(VectorRef&&) = delete;
      VectorRef& operator=(VectorRef const&) = delete;

      derived_type&& as_derived() && { return std::move(*static_cast<derived_type*>(this)); }

      derived_type& as_derived() & { return *static_cast<derived_type*>(this); }
      derived_type const& as_derived() const& { return *static_cast<derived_type const*>(this); }

      int size() const { return this->as_derived().size(); }
};

// Specialization of a VectorRef for a matrix that can be directly used in BLAS-like calls,
// ie it physically exists and can be addressed (although not necessarily in main memory,
// eg it might be on some other device such as a GPU), and has a fixed stride.
//
// The BaseType is the concrete type, which allows us to specialize blas calls for different
// devices.
// storage() returns an opaque type that contains the backing storage of the vector,
// and stride() is the array stride (the INC parameter in BLAS calls).

template <typename ValueType, typename BaseType, typename DerivedType = BaseType>
class BlasVector : public VectorRef<ValueType, BaseType, DerivedType>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using base_type          = BaseType;
      using storage_type       = typename blas_traits<base_type>::storage_type;
      using const_storage_type = typename blas_traits<base_type>::const_storage_type;

      BlasVector() = default;
      ~BlasVector() = default;
      BlasVector(BlasVector&& Other) = default;

      int stride() const { return this->as_derived().stride(); }

      storage_type storage() { return this->as_derived().storage(); }
      const_storage_type storage() const { return this->as_derived().storage(); }
};

// proxy class for the complex conjugate of a vector

template <typename T, typename BaseType>
class VectorConj : public VectorRef<T, BaseType, VectorConj<T, BaseType>>
{
   public:
      using base_type = BaseType;
      using value_type    = T;
      using storage_type       = typename blas_traits<base_type>::storage_type;
      using const_storage_type = typename blas_traits<base_type>::const_storage_type;

      VectorConj(base_type const& Base_) : Base(Base_) {}

      int size() const { return Base.size(); }

      base_type const& base() const { return Base; }

   private:
      base_type const& Base;
};

template <typename T, typename Base, typename Derived>
VectorConj<T, Derived>
conj(VectorRef<T, Base, Derived> const& x)
{
   return VectorConj<T, Derived>(x.as_derived());
}

// specialization for real types, conj() reduces to a no-op
template <typename Base, typename Derived>
Derived const&
conj(VectorRef<float, Base, Derived> const& x)
{
   return x.as_derived();
}

template <typename Base, typename Derived>
Derived const&
conj(VectorRef<double, Base, Derived> const& x)
{
   return x.as_derived();
}

#if defined(HAVE_FLOAT128)
template <typename Base, typename Derived>
Derived const&
conj(VectorRef<float128, Base, Derived> const& x)
{
   return x.as_derived();
}
#endif

// assignment

template <typename T, typename BaseType, typename U, typename V>
void assign(VectorRef<T, BaseType, U>& A, VectorRef<T, BaseType, V> const& B)
{
   vector_copy(B.as_derived(), A.as_derived());
}

template <typename T, typename BaseType>
void add(VectorRef<T, BaseType>& A, VectorRef<T, BaseType> const& B)
{
   vector_add(B.as_derived(), A.as_derived());
}

template <typename T, typename BaseType>
void subtract(VectorRef<T, BaseType>& A, VectorRef<T, BaseType> const& B)
{
   vector_add_scaled(-number_traits<T>::identity(), B.as_derived(), A.as_derived());
}

// expression template for alpha * A

template <typename T, typename BaseType, typename SubType>
class ScaledVector : public VectorRef<T, BaseType, ScaledVector<T, BaseType, SubType>>
{
   public:
      // confusing naming here - base_type is the original object, but the
      // template parameter is named SubType.  BaseType template parameter is the concrete vector type.
      using base_type = SubType;

      ScaledVector(T const& Factor_, base_type const& A_) : Factor(Factor_), A(A_) {}

      int size() const { return A.size(); }

      base_type const& base() const { return A; }
      T factor() const { return Factor; }

   private:
      T Factor;
      base_type const& A;
};

template <typename T, typename BaseType, typename SubType, typename X>
ScaledVector<decltype(safe_convert<T>(std::declval<X>())), BaseType, SubType>
operator*(X const& alpha, VectorRef<T, BaseType, SubType> const& M)
{
   return ScaledVector<T, BaseType, SubType>(safe_convert<T>(alpha), M.as_derived());
}

template <typename T, typename BaseType, typename SubType, typename U>
void assign(VectorRef<T, BaseType, U>& A, ScaledVector<T, BaseType, SubType> const& B)
{
   vector_copy_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename BaseType, typename SubType, typename U>
void assign(VectorRef<T, BaseType, U>&& A, ScaledVector<T, BaseType, SubType> const& B)
{
   vector_copy_scaled(B.factor(), B.base(), std::move(A).as_derived());
}

template <typename T, typename BaseType, typename SubType>
void add(VectorRef<T, BaseType>& A, ScaledVector<T, BaseType, SubType> const& B)
{
   vector_add_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename BaseType, typename SubType>
void subtract(VectorRef<T, BaseType>& A, ScaledVector<T, BaseType, SubType> const& B)
{
   vector_add_scaled(-B.factor(), B.base(), A.as_derived());
}

// expression template for alpha * op(A) * op(B)

template <typename T, typename BaseType, typename U, typename V>
struct MatrixVectorProduct : public VectorRef<T, BaseType, MatrixVectorProduct<T, BaseType, U, V>>
{
   template <typename MatrixBaseType>
   MatrixVectorProduct(MatrixRef<T, MatrixBaseType, U> const& A_, VectorRef<T, BaseType, V> const& B_)
      : Factor(number_traits<T>::identity()), A(A_.as_derived()), B(B_.as_derived()) {}

   template <typename MatrixBaseType>
   MatrixVectorProduct(T const& Factor_, MatrixRef<T, MatrixBaseType, U> const& A_, VectorRef<T, BaseType, V> const& B_)
      : Factor(Factor_), A(A_.as_derived()), B(B_.as_derived()) {}

   int size() const { return A.rows(); }
   T factor() const { return Factor; }

   T Factor;
   U const& A;
   V const& B;
};

template <typename T, typename MatrixBaseType, typename BaseType, typename U, typename V>
MatrixVectorProduct<T, BaseType, U, V>
operator*(MatrixRef<T, MatrixBaseType, U> const& A, VectorRef<T, BaseType, V> const& B)
{
   return MatrixVectorProduct<T, BaseType, U, V>(A.as_derived(), B.as_derived());
}

template <typename T, typename MatrixBaseType, typename MatrixSubType, typename BaseType, typename V>
MatrixVectorProduct<T, BaseType, MatrixBaseType, V>
operator*(ScaledMatrix<T, MatrixBaseType, MatrixSubType> const& A, VectorRef<T, BaseType, V> const& B)
{
   return MatrixVectorProduct<T, BaseType, MatrixSubType, V>(A.factor(), A.base(), B.as_derived());
}

template <typename T, typename MatrixBaseType, typename BaseType, typename SubType, typename U>
MatrixVectorProduct<T, BaseType, U, BaseType>
operator*(VectorRef<T, MatrixBaseType, U> const& A, ScaledVector<T, BaseType, SubType> const& B)
{
   return MatrixVectorProduct<T, BaseType, U, SubType>(B.factor(), A.as_derived(), B.base());
}

template <typename T, typename MatrixBaseType, typename MatrixSubType, typename BaseType, typename SubType>
MatrixVectorProduct<T, BaseType, MatrixBaseType, SubType>
operator*(ScaledMatrix<T, MatrixBaseType, MatrixSubType> const& A, ScaledVector<T, BaseType, SubType> const& B)
{
   return MatrixVectorProduct<T, BaseType, MatrixSubType, SubType>(A.factor() * B.factor(),
                                                                   A.base(), B.as_derived());
}

template <typename T, typename BaseType, typename U, typename V, typename X>
MatrixVectorProduct<decltype(safe_convert<T>(std::declval<X>())), BaseType, U, V>
operator*(X const& alpha, MatrixVectorProduct<T, BaseType, U, V> const& x)
{
   return MatrixVectorProduct<T, BaseType, U, V>(safe_convert<T>(alpha)*x.Factor, x.A, x.B);
}

// assignment

template <typename T, typename BaseType, typename Derived, typename U, typename V>
void assign(VectorRef<T, BaseType, Derived>& C, MatrixVectorProduct<T, BaseType, U, V> const& a)
{
   gemv(a.Factor, a.A, a.B, number_traits<T>::zero(), C.as_derived());
}

template <typename T, typename BaseType, typename Derived, typename U, typename V>
void add(VectorRef<T, BaseType, Derived>& C, MatrixVectorProduct<T, BaseType, U, V> const& a)
{
   gemv(a.Factor, a.A, a.B, number_traits<T>::identity(), C.as_derived());
}

template <typename T, typename BaseType, typename Derived, typename U, typename V>
void subtract(VectorRef<T, BaseType, Derived>& C, MatrixVectorProduct<T, BaseType, U, V> const& a)
{
   gemv(-a.Factor, a.A, a.B, number_traits<T>::identity(), C.as_derived());
}

} // namespace blas

#endif
