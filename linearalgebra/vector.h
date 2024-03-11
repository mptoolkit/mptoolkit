// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/vector.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
/* -*- C++ -*- $Id$

  vector.h

  A reference-counted vector class, with VectorRef (reference assignment semantics) and
  Vector (value assignment semantics).  A Vector and a VectorRef can share representations.

*/

#if !defined(VECTOR_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define VECTOR_H_DHJFKHURWEIRY4572893475489YRUI34897

#include "vectorinterface.h"
#include "vectoroperations.h"
#include "datablockreference.h"
#include "vecptriterator.h"
#include "noalias.h"

#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

#include "vectoraddition.h"
#include "vectortransform.h"
#include "crtp_vector.h"

namespace LinearAlgebra
{

// VectorRef

template <typename Scalar>
class VectorRef;

template <typename Scalar>
struct interface<VectorRef<Scalar> >
{
   typedef CONTIGUOUS_VECTOR(Scalar, VectorRef<Scalar>) type;
   typedef Scalar value_type;
};

template <typename Scalar>
class VectorRef : public VectorBase<VectorRef<Scalar> >
{
   public:
      typedef Scalar*       iterator;
      typedef Scalar const* const_iterator;
      typedef Scalar value_type;
      typedef Scalar* pointer;
      typedef Scalar const* const_pointer;
      typedef Scalar& reference;
      typedef Scalar const& const_reference;

      using VectorBase<VectorRef<Scalar> >::operator[];

      VectorRef(VectorRef const& V)
        : Block(V.Block) {}

      VectorRef& operator=(VectorRef const& V)
      {
         assign(*this, V);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorRef<Scalar>&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorRef<Scalar>&>::type
      operator=(NoAliasProxy<U> const& x);

      size_type size() const { return Block.size(); }

      iterator begin() { this->cow(); return Block.get(); }
      iterator end()   { this->cow(); return Block.get() + this->size(); }

      const_iterator begin() const { return Block.get(); }
      const_iterator end() const { return Block.get() + this->size(); }

      Scalar const& operator[](size_type n) const
         { DEBUG_RANGE_CHECK_OPEN(n,0U,this->size()); return *(Block.get() + n); }

      Scalar& operator[](size_type n)
         { DEBUG_RANGE_CHECK_OPEN(n,0U,this->size());
         this->cow(); return *(Block.get() + n); }

      // members defined in StrideVector
      Scalar* data() { this->cow(); return Block.get(); }
      Scalar const* data() const { return Block.get(); }

      difference_type stride() const { return 1; }
      size_type start() const { return 0; }

   protected:
      BlockReference<Scalar> const& block() const { return Block; }
      BlockReference<Scalar>& block() { return Block; }  // danger: this does not call cow()

      VectorRef() {}

      explicit VectorRef(size_type Size_)
        : Block(Size_) {}

      VectorRef(size_type Size_, Scalar const& Fill)
         : Block(Size_, Fill) {}
               //{ using LinearAlgebra::fill; fill(*this, Fill); }

     explicit VectorRef(BlockReference<Scalar> const& Block_) : Block(Block_) {}

     void resize(size_type NewSize);

     void swap(VectorRef& Other) { std::swap(Block, Other.Block); }

   protected:
      void cow() { Block.cow(); }

      BlockReference<Scalar> Block;

   template <typename U> friend class VectorRef;
   template <typename U> friend class Vector;
};

template <typename T>
struct is_const_proxy<VectorRef<T> > : boost::is_const<T> {};

template <typename T>
struct is_const_proxy<VectorRef<T> const> : boost::mpl::true_ {};

template <typename T>
struct is_mutable_proxy<VectorRef<T> > : boost::mpl::not_<boost::is_const<T> > {};

// Size

template <typename T>
struct Size<VectorRef<T>, DENSE_VECTOR(T, VectorRef<T>)>
{
   typedef size_type result_type;
   typedef VectorRef<T> argument_type;

   size_type operator()(VectorRef<T> const& v) const { return v.size(); }
};

// iterators

template <typename T>
struct Iterate<VectorRef<T>&>
{
   typedef VectorRef<T>& argument_type;
   typedef VecPtrIterator<T> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size(), 0);
   }
};

// FIXME: the following specialization probably should not appear - a VectorRef<T>
// is a mutable proxy so it should be covered by the previous one.
// But removing it causes the regression tests to go crazy.  Something for a rainy day?
#if 1
template <typename T>
struct Iterate<VectorRef<T> >
{
   typedef VectorRef<T> argument_type;
   typedef VecPtrIterator<T const> result_type;
   result_type operator()(argument_type const& x) const
   {
      return result_type(x.data(), x.size(), 0);
   }
};
#endif

template <typename T>
struct Iterate<VectorRef<T> const>
{
   typedef VectorRef<T> const& argument_type;
   typedef VecPtrIterator<T const> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size(), 0);
   }
};

//
// Vector
//
// vector with copy-on-write semantics.  This is the 'default' vector class.
//

template <typename Scalar>
class Vector;

template <typename Scalar>
struct interface<Vector<Scalar> >
{
   typedef CONTIGUOUS_VECTOR(Scalar, Vector<Scalar>) type;
   typedef Scalar value_type;
};

template <class Scalar>
class Vector : public VectorRef<Scalar>
{
   public:
      typedef Scalar*       iterator;
      typedef Scalar const* const_iterator;
      typedef Scalar value_type;

      Vector() {}

      explicit Vector(size_type Size_);

      Vector(size_type Size_, Scalar const& fill_);

      template <typename Iter>
      Vector(Iter first, Iter last,
             typename boost::enable_if<
             boost::mpl::not_<boost::is_arithmetic<Iter> > >::type* dummy = 0);

      Vector(Vector<Scalar> const& V)
        : VectorRef<Scalar>(V.block().copy()) {}

      Vector(VectorRef<Scalar> const& V)
        : VectorRef<Scalar>(V.block().copy()) {}

      template <typename U>
      Vector(U const& x, typename boost::enable_if<is_vector<U> >::type* dummy = 0);

      Vector<Scalar>& operator=(Vector<Scalar> const& V);

      template <typename U>
      typename boost::enable_if<is_vector<U>, Vector<Scalar>&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_vector<U>, Vector<Scalar>&>::type
      operator=(NoAliasProxy<U> const& x);

      void resize(size_type NewSize)
        { this->VectorRef<Scalar>::resize(NewSize); }

      void swap(Vector& Other)
         { this->cow(); Other.cow(); this->VectorRef<Scalar>::swap(Other); }
};

// iterators

template <typename T>
struct Iterate<Vector<T>&>
{
   typedef Vector<T>& argument_type;
   typedef VecPtrIterator<T> result_type;
   result_type operator()(argument_type x) const; //__attribute((always_inline));
   //   inline result_type operator()(argument_type x) const
   //   {
   //      return result_type(x.data(), x.size(), 0);
   //   }
};

template <typename T>
inline
typename Iterate<Vector<T>&>::result_type
Iterate<Vector<T>&>::operator()(argument_type x) const
{
   return result_type(x.data(), x.size(), 0);
}

template <typename T>
struct Iterate<Vector<T> >
{
   typedef Vector<T> argument_type;
   typedef VecPtrIterator<T const> result_type;
   inline result_type operator()(argument_type const& x) const
   {
      return result_type(x.data(), x.size(), 0);
   }
};

// resize

template <typename T>
struct Resize<Vector<T>&>
{
   typedef void result_type;
   typedef Vector<T>& first_argument_type;
   typedef size_type second_argument_type;
   void operator()(Vector<T>& v, size_type n) const
   {
      v.resize(n);
   }
};

// swap

template <typename T>
struct Swap<Vector<T>&, Vector<T>&>
{
   typedef void result_type;
   typedef Vector<T>& first_argument_type;
   typedef Vector<T>& second_argument_type;
   result_type operator()(Vector<T>& t, Vector<T>& u) const
   {
      t.swap(u);
   }
};

} // namespace LinearAlgebra

#include "vector.cc"

#endif
