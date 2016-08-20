// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectormemproxy.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  vectormemproxy.h

  A proxy class for treating a memory range as a vector.  This includes
  a normalization correction for non-abelian quantum numbers.

  Created 2005-01-07 Ian McCulloch

*/

#if !defined(VECTORMEMPROXY_H_SCHJUJHVUIREGHULIY3489Y4839Y54398)
#define VECTORMEMPROXY_H_SCHJUJHVUIREGHULIY3489Y4839Y54398

#include "vectoroperationsbase.h"
#include "vecptriterator.h"
#include "noalias.h"
#include "crtp_vector.h"
#include "slice.h"
#include <complex>

namespace LinearAlgebra
{

template <typename T, typename Stride = boost::mpl::int_<1> >
class VectorMemProxy;

template <typename T>
class VectorMemProxy<T, boost::mpl::int_<1> >
   : public VectorBase<VectorMemProxy<T, boost::mpl::int_<1> > >
{
   public:
      typedef typename boost::remove_const<T>::type value_type;
      typedef T* iterator;
      typedef T const* const_iterator;
      typedef T& reference;
      typedef T const& const_reference;

      // an easy way to declare this class as being a proxy reference type
      // is to define proxy to mpl::true_, and specify what the
      // corresponding const reference type is.
      typedef boost::is_const<T>                    const_proxy;
      typedef boost::mpl::not_<boost::is_const<T> > proxy;

      typedef VectorMemProxy<T const, boost::mpl::int_<1> > const_type;

      VectorMemProxy() : Data_(NULL), Size_(0) {}

      VectorMemProxy(T* Data, size_type Size) : Data_(Data), Size_(Size) {}

      // copy-assignment
      VectorMemProxy& operator=(VectorMemProxy const& Other)
      {
         assign(*this, Other);
         return *this;
      }

      // Assignment has reference semantics
      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy<T>& >::type
      operator=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         assign(*this, x);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy<T>& >::type
      operator=(NoAliasProxy<U> const& x)
      {
         assign(*this, x.value());
         return *this;
      }

      size_type size() const { return Size_; }

      iterator begin() { return Data_; }
      iterator end() { return Data_ + Size_; }

      const_iterator begin() const { return Data_; }
      const_iterator end() const { return Data_ + Size_; }

   //      reference operator[](size_type n) { return Data_[n]; }
   //      const_reference operator[](size_type n) const { return Data_[n]; }

      T* data() { return Data_; }
      T const* data() const { return Data_; }
      difference_type stride() const { return 1; }


      operator VectorMemProxy<T const>() const { return VectorMemProxy<T const>(Data_, Size_); }

      // computed assignment must be members here, so we can use them on temporaries

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy<T>& >::type
      operator+=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         add(*this, x);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy<T>& >::type
      operator+=(NoAliasProxy<U> const& x)
      {
         add(*this, x.value());
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy<T>& >::type
      operator-=(NoAliasProxy<U> const& x)
      {
         subtract(*this, x.value());
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy<T>& >::type
      operator-=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         subtract(*this, x);
         return *this;
      }

   private:
      T* Data_;
      size_type Size_;
};

template <typename T>
class VectorMemProxy<T, tagVariable>
   : public VectorBase<VectorMemProxy<T, tagVariable> >
{
   public:

      typedef typename basic_type<T>::type value_type;
      typedef typename make_reference<T>::type reference;
      typedef typename make_const_reference<T>::type const_reference;

      // an easy way to declare this class as being a proxy reference type
      // is to define proxy to mpl::true_, and specify what the
      // corresponding const reference type is.
      typedef boost::is_const<T>                    const_proxy;
      typedef boost::mpl::not_<boost::is_const<T> > proxy;

      typedef VectorMemProxy<T const, tagVariable> const_type;

      VectorMemProxy() : Data_(NULL), Size_(0) {}

      VectorMemProxy(T* Data, size_type Size, difference_type Stride)
         : Data_(Data), Size_(Size), Stride_(Stride) {}

      // copy-assignment
      VectorMemProxy& operator=(VectorMemProxy const& Other)
      {
         assign(*this, Other);
         return *this;
      }

      // Assignment has reference semantics
      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator=(U const& x)
      {
         assign_copy(*this, x);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy& >::type
      operator=(NoAliasProxy<U> const& x)
      {
         assign(*this, x.value());
         return *this;
      }

      size_type size() const { return Size_; }

   //      reference operator[](size_type n) { return Data_[Stride_ * n]; }
   //      const_reference operator[](size_type n) const { return Data_[Stride_ * n]; }

      T* data() { return Data_; }
      T const* data() const { return Data_; }
      difference_type stride() const { return Stride_; }

#if 0
   // does this do anything?  icc warns that it is never used
      operator VectorMemProxy<T const, tagVariable>() const
   { return VectorMemProxy<T const, tagVariable>(Data_, Size_, Stride_); }
#endif

      // computed assignment must be members here, so we can use them on temporaries

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator+=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         add(*this, x);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy& >::type
      operator+=(NoAliasProxy<U> const& x)
      {
         add(*this, x.value());
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy& >::type
      operator-=(NoAliasProxy<U> const& x)
      {
         subtract(*this, x.value());
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy& >::type
      operator-=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         subtract(*this, x);
         return *this;
      }

      // operator *=
      template <typename U>
      VectorMemProxy&
      operator*=(U const& y)
      {
         multiply(*this, y);
         return *this;
      }

   private:
      T* Data_;
      difference_type Size_;
      difference_type Stride_;
};

template <typename T, int Stride_>
class VectorMemProxy<T, boost::mpl::int_<Stride_> >
: public VectorBase<VectorMemProxy<T, boost::mpl::int_<Stride_> > >
{
   public:
      typedef typename boost::remove_const<T>::type value_type;
      typedef T& reference;
      typedef T const& const_reference;

      // an easy way to declare this class as being a proxy reference type
      // is to define proxy to mpl::true_, and specify what the
      // corresponding const reference type is.
      typedef boost::is_const<T>                    const_proxy;
      typedef boost::mpl::not_<boost::is_const<T> > proxy;

      typedef VectorMemProxy<T const, boost::mpl::int_<Stride_> > const_type;

      VectorMemProxy() : Data_(NULL), Size_(0) {}

      VectorMemProxy(T* Data, size_type Size, difference_type Stride = Stride_)
         : Data_(Data), Size_(Size) { PRECONDITION_EQUAL(Stride, Stride_); }

      // copy-assignment
      VectorMemProxy& operator=(VectorMemProxy const& Other)
      {
         assign(*this, Other);
         return *this;
      }

      // Assignment has reference semantics
      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         assign(*this, x);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator=(NoAliasProxy<U> const& x)
      {
         assign(*this, x.value());
         return *this;
      }

      size_type size() const { return Size_; }

   //      reference operator[](size_type n) { return Data_[Stride_ * n]; }
   //      const_reference operator[](size_type n) const { return Data_[Stride_ * n]; }

      T* data() { return Data_; }
      T const* data() const { return Data_; }
      difference_type stride() const { return Stride_; }
      static difference_type const static_stride = Stride_;

      operator VectorMemProxy<T const,  boost::mpl::int_<Stride_> >() const
   { return VectorMemProxy<T const,  boost::mpl::int_<Stride_> >(Data_, Size_, Stride_); }

      // computed assignment must be members here, so we can use them on temporaries

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator+=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         add(*this, x);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator+=(NoAliasProxy<U> const& x)
      {
         add(*this, x.value());
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator-=(NoAliasProxy<U> const& x)
      {
         subtract(*this, x.value());
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorMemProxy&>::type
      operator-=(U const& x)
      {
         // TODO: the vector type should be something better here!
         //      std::vector<value_type> Temp(x.size());
         //      noalias(Temp) = x;
         subtract(*this, x);
         return *this;
      }

      // operator *=
      template <typename U>
      VectorMemProxy&
      operator*=(U const& y)
      {
         multiply(*this, y);
         return *this;
      }

   private:
      T* Data_;
      difference_type Size_;
};

#if 0
 // the declaration of VectorMemProxy<T>::proxy and VectorMemProxy<T>::const_type
 // means that we do not have to specify this:
template <typename T>
struct make_reference<VectorMemProxy<T> >
{
   typedef VectorMemProxy<T> type;
};

template <typename T>
struct convert_to_const_reference<VectorMemProxy<T> >
{
   typedef VectorMemProxy<T const> type;
};
#endif

//
// make_vec - a helper function to construct a VectorMemProxy
// from a pointer and size, or a pointer range.
//

template <typename T>
inline
VectorMemProxy<T> make_vec(T* Ptr, size_type Size)
{
   return VectorMemProxy<T>(Ptr, Size);
}

template <typename T>
inline
VectorMemProxy<T const> make_vec(T const* Ptr, size_type Size)
{
   return VectorMemProxy<T const>(Ptr, Size);
}

template <typename T>
inline
VectorMemProxy<T> make_vec(T* Begin, T* End)
{
   return VectorMemProxy<T>(Begin, size_type(End-Begin));
}

template <typename T>
inline
VectorMemProxy<T const> make_vec(T const* Begin, T const* End)
{
   return VectorMemProxy<T const>(Begin, size_type(End-Begin));
}

// interface

// general stride
template <typename Scalar, typename Stride>
struct interface<VectorMemProxy<Scalar, Stride> >
{
   typedef typename boost::remove_const<Scalar>::type value_type;
   typedef STRIDE_VECTOR(value_type, void) type;
};

// unit stride
template <typename Scalar>
struct interface<VectorMemProxy<Scalar> >
{
   typedef typename boost::remove_const<Scalar>::type value_type;
   typedef CONTIGUOUS_VECTOR(value_type, void) type;
};

// iterators
// we need 4 variations to cover iterators and const iterators
// over non-const and const proxies.

template <typename T, typename Stride>
struct Iterate<VectorMemProxy<T, Stride>&>
{
   typedef VectorMemProxy<T, Stride>& argument_type;
   typedef VecPtrIterator<T, Stride> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size(), 0, x.stride());
   }
};

template <typename T, typename Stride>
struct Iterate<VectorMemProxy<T, Stride> >
{
   typedef VectorMemProxy<T, Stride> const& argument_type;
   typedef VecPtrIterator<T const, Stride> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size(), 0, x.stride());
   }
};

template <typename T, typename Stride>
struct Iterate<VectorMemProxy<T const, Stride>&>
{
#if 0
   typedef VectorMemProxy<T const, Stride> const& argument_type;
   typedef VecPtrIterator<T const, Stride> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size(), 0, x.stride());
   }
#endif
};

// defaults for slices of STRIDE_VECTOR and CONTIGUOUS_VECTOR

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T, Slice, CONTIGUOUS_VECTOR(Tv, Ti)>
{
   typedef T first_argument_type;
   typedef Slice second_argument_type;
   typedef VectorMemProxy<Tv const, tagVariable> result_type;

   result_type operator()(first_argument_type const& v, Slice const& s) const
   {
      return result_type(data(v) + s.start(), s.size(), s.stride());
   }
};

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T&, Slice, CONTIGUOUS_VECTOR(Tv, Ti)>
{
   typedef T& first_argument_type;
   typedef Slice second_argument_type;
   typedef VectorMemProxy<Tv, tagVariable> result_type;

   result_type operator()(first_argument_type v, Slice const& s) const
   {
      return result_type(data(v) + s.start(), s.size(), s.stride());
   }
};

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T, Slice, STRIDE_VECTOR(Tv, Ti)>
{
   typedef T first_argument_type;
   typedef Slice second_argument_type;
   typedef VectorMemProxy<Tv const, tagVariable> result_type;

   result_type operator()(first_argument_type const& v, Slice const& s) const
   {
      return result_type(data(v) + stride(v) * s.start(), s.size(),
                         stride(v) * s.stride());
   }
};

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T&, Slice, STRIDE_VECTOR(Tv, Ti)>
{
   typedef T& first_argument_type;
   typedef Slice second_argument_type;
   typedef VectorMemProxy<Tv, tagVariable> result_type;

   result_type operator()(first_argument_type v, Slice const& s) const
   {
      return result_type(data(v) +  stride(v) * s.start(), s.size(),
                         stride(v) * s.stride());
   }
};



// defaults for ranges of STRIDE_VECTOR and CONTIGUOUS_VECTOR

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T, Range, CONTIGUOUS_VECTOR(Tv, Ti)>
{
   typedef T first_argument_type;
   typedef Range second_argument_type;
   typedef VectorMemProxy<Tv const> result_type;

   result_type operator()(first_argument_type const& v, Range const& r) const
   {
      return result_type(data(v) + r.first(), r.size());
   }
};

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T&, Range, CONTIGUOUS_VECTOR(Tv, Ti)>
{
   typedef T& first_argument_type;
   typedef Range second_argument_type;
   typedef VectorMemProxy<Tv> result_type;

   result_type operator()(first_argument_type v, Range const& r) const
   {
      return result_type(data(v) + r.first(), r.size());
   }
};

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T, Range, STRIDE_VECTOR(Tv, Ti)>
{
   typedef T first_argument_type;
   typedef Range second_argument_type;
   typedef VectorMemProxy<Tv const, tagVariable> result_type;

   result_type operator()(first_argument_type const& v, Range const& r) const
   {
      return result_type(data(v) + stride(v) * r.first(), r.size(), stride(v));
   }
};

template <typename T, typename Tv, typename Ti>
struct VectorBracketInterface<T&, Range, STRIDE_VECTOR(Tv, Ti)>
{
   typedef T& first_argument_type;
   typedef Range second_argument_type;
   typedef VectorMemProxy<Tv, tagVariable> result_type;

   result_type operator()(first_argument_type v, Range const& r) const
   {
      return result_type(data(v) + stride(v) * r.first(), r.size(), stride(v));
   }
};

// A glorious hack to express real() and imag() of STRIDE_VECTOR(std::complex<T>) as
// STRIDE_VECTOR(T)

// Real for CONTIGUOUS_VECTOR

template <typename T, typename Tv, typename Ti>
struct RealInterface<T, CONTIGUOUS_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv const, boost::mpl::int_<2> > result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
      { return result_type(reinterpret_cast<Tv const*>(data(x)), size(x)); }
};

template <typename T, typename Tv, typename Ti>
struct RealInterface<T&, CONTIGUOUS_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv, boost::mpl::int_<2> > result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const
      { return result_type(reinterpret_cast<Tv*>(data(x)), size(x)); }
};

// Real for STRIDE_VECTOR

template <typename T, typename Tv, typename Ti>
struct RealInterface<T, STRIDE_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv const, tagVariable> result_type;
   typedef T const& argument_type;
   result_type operator()(T const& x) const
      { return result_type(reinterpret_cast<Tv const*>(data(x)), size(x), 2 * stride(x)); }
};

template <typename T, typename Tv, typename Ti>
struct RealInterface<T&, STRIDE_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv, tagVariable> result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const
      { return result_type(reinterpret_cast<Tv*>(data(x)), size(x), 2 * stride(x)); }
};

// Imag for CONTIGUOUS_VECTOR

template <typename T, typename Tv, typename Ti>
struct ImagInterface<T, CONTIGUOUS_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv const, boost::mpl::int_<2> > result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
      { return result_type(reinterpret_cast<Tv const*>(data(x))+1, size(x)); }
};

template <typename T, typename Tv, typename Ti>
struct ImagInterface<T&, CONTIGUOUS_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv, boost::mpl::int_<2> > result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const
      { return result_type(reinterpret_cast<Tv*>(data(x))+1, size(x)); }
};

// Imag for STRIDE_VECTOR

template <typename T, typename Tv, typename Ti>
struct ImagInterface<T, STRIDE_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv const, tagVariable> result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
      { return result_type(reinterpret_cast<Tv const*>(data(x))+1, size(x), 2 * stride(x)); }
};

template <typename T, typename Tv, typename Ti>
struct ImagInterface<T&, STRIDE_VECTOR(std::complex<Tv>, Ti)>
{
   typedef VectorMemProxy<Tv, tagVariable> result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const
      { return result_type(reinterpret_cast<Tv*>(data(x))+1, size(x), 2 * stride(x)); }
};

} // namespace LinearAlgebra

#endif
