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

#if !defined(MPTOOLKIT_BLAS_VECTOR_VIEW_H)
#define MPTOOLKIT_BLAS_VECTOR_VIEW_H

#include "vectorref.h"

// vector_view is a proxy class that interprets strided 'view' of a
// buffer as a vector.  vector_view can be used as an l-value.

namespace blas
{

//
// normal_vector_view
// stride-1 version of a vector_view
//

template <typename T, typename Tag>
class const_normal_vector_view;

template <typename T, typename Tag>
class const_vector_view;

template <typename T, typename Tag>
class normal_vector_view;

template <typename T, typename Tag>
class vector_view;

template <typename T, typename Tag>
class normal_vector_view : public NormalVectorProxy<T, normal_vector_view<T, Tag>, Tag>
{
   public:
      using value_type         = T;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      normal_vector_view() = delete;

      normal_vector_view(int Size_, storage_type const& Ptr_)
         : Size(Size_), Ptr(Ptr_) {}

      normal_vector_view(normal_vector_view&& Other) = default;

      normal_vector_view(normal_vector_view const&) = delete;

      normal_vector_view& operator=(normal_vector_view&&) = delete;

      ~normal_vector_view() = default;

      template <typename U>
      normal_vector_view&& operator=(blas::VectorRef<T, U, Tag> const& E) &&
      {
	 assign(static_cast<normal_vector_view&&>(*this), E.as_derived());
	 return std::move(*this);
      }

      // TODO: move assignment from a Vector

      template <typename U>
      normal_vector_view&& operator+=(blas::VectorRef<T, U, Tag> const& E) &&
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      normal_vector_view&& operator-=(blas::VectorRef<T, U, Tag> const& E) &&
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      constexpr int stride() const { return 1; }

      int size() const { return Size; }

      storage_type storage() && { return Ptr; }
      const_storage_type storage() const& { return Ptr; }

   private:
      int Size;
      storage_type Ptr;

      friend class const_normal_vector_view<T, Tag>;
      friend class const_vector_view<T, Tag>;
      friend class vector_view<T, Tag>;
};

template <typename T, typename Tag>
class const_normal_vector_view : public blas::NormalVector<T,normal_vector_view<T, Tag>, Tag>
{
   public:
      using value_type         = T;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      const_normal_vector_view() = delete;

      const_normal_vector_view(int Size_, const_storage_type const& Ptr_)
         : Size(Size_), Ptr(Ptr_) {}

      const_normal_vector_view(const_normal_vector_view&& Other) = default;

      const_normal_vector_view(normal_vector_view<T, Tag>&& Other)
         : Size(Other.Size), Ptr(std::move(Other.Ptr)) {}

      const_normal_vector_view(const_normal_vector_view const&) = delete;

      const_normal_vector_view& operator=(const_normal_vector_view&&) = delete;

      ~const_normal_vector_view() = default;

      constexpr int stride() const { return 1; }

      int size() const { return Size; }

      const_storage_type storage() const { return Ptr; }

   private:
      int Size;
      const_storage_type Ptr;

      friend class const_vector_view<T, Tag>;
};

//
// vector_view
//

template <typename T, typename Tag>
class vector_view : public BlasVectorProxy<T, vector_view<T, Tag>, Tag>
{
   public:
      using value_type         = T;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      vector_view() = delete;

      vector_view(int Size_, int Stride_, storage_type const& Ptr_)
         : Size(Size_), Stride(Stride_), Ptr(Ptr_) {}

      vector_view(vector_view&& Other) = default;

      vector_view(normal_vector_view<T, Tag>&& Other)
         : Size(Other.Size), Stride(1), Ptr(std::move(Other.Ptr)) {}

      vector_view(vector_view const&) = delete;

      vector_view&& operator=(vector_view&& E) &&
      {
	 assign(static_cast<vector_view&&>(*this), E.as_derived());
	 return std::move(*this);
      }

      ~vector_view() = default;

      template <typename U>
      vector_view&& operator=(blas::VectorRef<T, U, Tag> const& E) &&
      {
	 assign(static_cast<vector_view&&>(*this), E.as_derived());
	 return std::move(*this);
      }

      template <typename U>
      vector_view&& operator+=(blas::VectorRef<T, U, Tag> const& E) &&
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      vector_view&& operator-=(blas::VectorRef<T, U, Tag> const& E) &&
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      int stride() const { return Stride; }

      int size() const { return Size; }

      storage_type storage() && { return Ptr; }
      const_storage_type storage() const& { return Ptr; }

   private:
      int Size;
      int Stride;
      storage_type Ptr;

      friend class const_vector_view<T, Tag>;
};

template <typename T, typename Tag>
class const_vector_view : public blas::BlasVector<T,vector_view<T, Tag>, Tag>
{
   public:
      using value_type         = T;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;

      const_vector_view() = delete;

      const_vector_view(int Size_, int Stride_, const_storage_type const& Ptr_)
         : Size(Size_), Stride(Stride_), Ptr(Ptr_) {}

      const_vector_view(const_vector_view&& Other) = default;

      const_vector_view(vector_view<T, Tag>&& Other)
         : Size(Other.Size), Stride(Other.Stride), Ptr(std::move(Other.Ptr)) {}

      const_vector_view(normal_vector_view<T, Tag>&& Other)
         : Size(Other.Size), Stride(1), Ptr(std::move(Other.Ptr)) {}

      const_vector_view(const_normal_vector_view<T, Tag>&& Other)
         : Size(Other.Size), Stride(1), Ptr(std::move(Other.Ptr)) {}

      const_vector_view(const_vector_view const&) = delete;

      const_vector_view& operator=(const_vector_view&&) = delete;

      ~const_vector_view() = default;

      int stride() const { return Stride; }

      int size() const { return Size; }

      const_storage_type storage() const { return Ptr; }

   private:
      int Size;
      int Stride;
      const_storage_type Ptr;
};

#if defined(USE_PSTREAM)

template <typename T, int Format>
PStream::opstreambuf<Format>&
operator<<(PStream::opstreambuf<Format>& out, normal_vector_view<T, cpu_tag> const& x)
{
   typedef typename PStream::opstreambuf<Format>::size_type st;
   st s = x.size();
   out << s;
   for (auto const& i : x)
   {
      out << i;
   }
   return out;
}

template <typename T, int Format>
PStream::opstreambuf<Format>&
operator<<(PStream::opstreambuf<Format>& out, const_normal_vector_view<T, cpu_tag> const& x)
{
   typedef typename PStream::opstreambuf<Format>::size_type st;
   st s = x.size();
   out << s;
   for (auto const& i : x)
   {
      out << i;
   }
   return out;
}

template <typename T, int Format>
PStream::opstreambuf<Format>&
operator<<(PStream::opstreambuf<Format>& out, vector_view<T, cpu_tag> const& x)
{
   typedef typename PStream::opstreambuf<Format>::size_type st;
   st s = x.size();
   out << s;
   for (auto const& i : x)
   {
      out << i;
   }
   return out;
}

template <typename T, int Format>
PStream::opstreambuf<Format>&
operator<<(PStream::opstreambuf<Format>& out, const_vector_view<T, cpu_tag> const& x)
{
   typedef typename PStream::opstreambuf<Format>::size_type st;
   st s = x.size();
   out << s;
   for (auto const& i : x)
   {
      out << i;
   }
   return out;
}

#endif

} // namespace blas

#endif
