// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/iteratortypes.h
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
  Various interator adaptors for matrix & vector expressions.

  Created 2004-04-05 Ian McCulloch
*/

#if !defined(ITERATORTYPES_H_JHFU8R4Y8Y89YR98RY4WT598YT458Y4)
#define ITERATORTYPES_H_JHFU8R4Y8Y89YR98RY4WT598YT458Y4

#include <cstddef>
#include "common/numerics.h"
#include "common/trace.h"
#include "slice.h"

namespace LinearAlgebra
{

template <typename Iter>
struct IteratorTraits;

template <typename Scalar>
struct IteratorTraits<Scalar*>
{
   typedef Scalar        value_type;
   typedef Scalar&       result_type;
   typedef Scalar const& const_result_type;
};

template <typename Scalar>
struct IteratorTraits<Scalar const*>
{
   typedef Scalar        value_type;
   typedef Scalar const& result_type;
   typedef Scalar const& const_result_type;
};

typedef std::size_t size_type;
typedef std::ptrdiff_t difference_type;

//
// StrideIterator
//

template <typename BaseIterType>
struct StrideIterator;

template <typename BaseIterType>
struct IteratorTraits<StrideIterator<BaseIterType> >
{
   typedef IteratorTraits<BaseIterType>                 base_traits_type;
   typedef typename base_traits_type::value_type        value_type;
   typedef typename base_traits_type::result_type       result_type;
   typedef typename base_traits_type::const_result_type const_result_type;
};

// we never want to instantiate a StrideIterator<StrideIterator<T> >,
// so make it an error to do so
template <typename BaseIterType>
struct IteratorTraits<StrideIterator<StrideIterator<BaseIterType> > >
{
};

template <typename BaseIterType>
class StrideIterator
{
   public:
      typedef BaseIterType base_type;
      typedef IteratorTraits<base_type>               traits_type;
      typedef typename traits_type::value_type        value_type;
      typedef typename traits_type::result_type       result_type;

      StrideIterator() {}
      StrideIterator(BaseIterType I_, difference_type Stride_)
        : I(I_), Stride(Stride_) {}

      StrideIterator(BaseIterType I_, difference_type Offset, difference_type Stride_)
        : I(I_), Stride(Stride_) { I += Offset; }

      template <typename Other>
      StrideIterator(StrideIterator<Other> const& Iter)
        : I(Iter.base()), Stride(Iter.stride()) {}

      StrideIterator& operator++() { I += Stride; return *this; }
      StrideIterator& operator--() { I -= Stride; return *this; }

      StrideIterator& operator+=(difference_type n) { I += Stride * n; return *this; }
      StrideIterator& operator-=(difference_type n) { I -= Stride * n; return *this; }

      result_type operator*() const { return *I; }

      result_type operator[](difference_type n) const { return I[Stride * n]; }

      StrideIterator operator+(difference_type n) const { return StrideIterator(I, Stride * n, Stride); }
      StrideIterator operator-(difference_type n) const { return StrideIterator(I, -Stride * n, Stride); }

      base_type& base() { return I; }
      base_type const& base() const { return I; }
      difference_type stride() const { return Stride; }

   private:
      BaseIterType I;
      difference_type Stride;
};

template <typename B1, typename B2>
inline
bool operator==(StrideIterator<B1> const& i1,
                StrideIterator<B2> const& i2)
{
   DEBUG_PRECONDITION(i1.stride() == i2.stride())(i1.stride())(i2.stride());
   return i1.base() == i2.base();
}

template <typename B1, typename B2>
inline
bool operator!=(StrideIterator<B1> const& i1,
                StrideIterator<B2> const& i2)
{
   DEBUG_PRECONDITION(i1.stride() == i2.stride())(i1.stride())(i2.stride());
   return i1.base() != i2.base();
}

template <typename B1, typename B2>
inline
difference_type operator-(StrideIterator<B1> const& i1,
                          StrideIterator<B2> const& i2)
{
   DEBUG_PRECONDITION(i1.stride() == i2.stride())(i1.stride())(i2.stride());
   DEBUG_CHECK((i1.base() - i2.base()) % i1.stride() == 0);
   return (i1.base() - i2.base()) / i1.stride();
}

// not all iterator types have operator- defined.  But it is useful anyway for those that do.
template <typename Base>
inline
difference_type operator-(StrideIterator<Base> const& i1, StrideIterator<Base> const& i2)
{
   DEBUG_PRECONDITION(i1.stride() == i2.stride());
   DEBUG_CHECK((i1.base() - i2.base()) % i1.stride() == 0);
   return (i1.base() - i2.base()) / i1.stride();
}

//
// TransformIterator
//

template <typename BaseIterType, typename F>
class TransformIterator;

template <typename BaseIterType, typename F>
struct IteratorTraits<TransformIterator<BaseIterType, F> >
{
   typedef F                                     functor_type;
   typedef BaseIterType                          base_type;
   typedef IteratorTraits<base_type>             base_traits_type;
   typedef typename base_traits_type::value_type base_value_type;

   typedef typename functor_type::value_type        value_type;
   typedef typename functor_type::result_type       const_result_type;
};

template <typename BaseIterType, typename F>
class TransformIterator
{
   public:
      typedef BaseIterType                                        base_type;
      typedef F                                                   functor_type;
      typedef IteratorTraits<TransformIterator<BaseIterType, F> > traits_type;
      typedef typename traits_type::value_type                    value_type;
      typedef typename traits_type::const_result_type             const_result_type;

      TransformIterator() {}

      TransformIterator(base_type const& I_, functor_type const& f_ = functor_type())
        : I(I_), f(f_) {}

      TransformIterator& operator++() { ++I; return *this; }
      TransformIterator& operator--() { --I; return *this; }

      TransformIterator& operator+=(difference_type n) { I += n; return *this; }
      TransformIterator& operator-=(difference_type n) { I -= n; return *this; }

      const_result_type operator*() const { return f(*I); }

      const_result_type operator[](difference_type n) const { return f(I[n]); }

      bool operator==(TransformIterator const& Other) const { return I == Other.I; }
      bool operator!=(TransformIterator const& Other) const { return I != Other.I; }

      TransformIterator operator+(difference_type n) const { return TransformIterator(I+n, f); }
      TransformIterator operator-(difference_type n) const { return TransformIterator(I+n, f); }

      // functions for TransformIterator-aware optimizations
      functor_type const& func() const { return f; }
      base_type const& base() const { return I; }

   private:
      base_type I;
      functor_type f;
};

template <typename B1, typename B2, typename F>
inline
difference_type operator-(TransformIterator<B1, F> const& i1, TransformIterator<B2, F> const& i2)
{
   return i1.base() - i2.base();
}

//
// BinaryIterator
//

template <typename BaseIter1, typename BaseIter2, typename F>
class BinaryIterator;

template <typename BaseIter1, typename BaseIter2, typename F>
struct IteratorTraits<BinaryIterator<BaseIter1, BaseIter2, F> >
{
   typedef F                                  functor_type;
   typedef BaseIter1                          base_type1;
   typedef BaseIter2                          base_type2;
   typedef typename functor_type::value_type  value_type;
   typedef typename functor_type::result_type const_result_type;
};

template <typename BaseIter1, typename BaseIter2, typename F>
class BinaryIterator
{
   public:
      typedef BaseIter1                                                base_type1;
      typedef BaseIter2                                                base_type2;
      typedef F                                                        functor_type;
      typedef IteratorTraits<BinaryIterator<BaseIter1, BaseIter2, F> > traits_type;
      typedef typename traits_type::const_result_type                  const_result_type;

      BinaryIterator() {}

      BinaryIterator(base_type1 const& i1_, base_type2 const& i2_,
                     functor_type f_ = functor_type())
        : i1(i1_), i2(i2_), f(f_) {}

      BinaryIterator& operator++() { ++i1; ++i2; return *this; }
      BinaryIterator& operator--() { --i1; --i2; return *this; }

      BinaryIterator& operator+=(difference_type n) { i1 += n; i2 += n; return *this; }
      BinaryIterator& operator-=(difference_type n) { i2 -= n; i2 += n; return *this; }

      const_result_type operator*() const { return f(*i1, *i2); }

      const_result_type operator[](difference_type n) const { return f(i1[n], i2[n]); }

      bool operator==(BinaryIterator const& Other) const
         { return i1 == Other.i1 && i2 == Other.i2; }

      bool operator!=(BinaryIterator const& Other) const
         { return i1 != Other.i1 || i2 != Other.i2; }

      BinaryIterator operator+(difference_type n) const { return BinaryIterator(i1+n, i2+n, f); }
      BinaryIterator operator-(difference_type n) const { return BinaryIterator(i1+n, i2+n, f); }

      // functions for BinaryIterator-aware optimizations
      functor_type const& func() const { return f; }
      base_type1 const& base1() const { return i1; }
      base_type2 const& base2() const { return i2; }

   private:
      base_type1 i1;
      base_type2 i2;
      functor_type f;
};

template <typename A1, typename A2, typename B1, typename B2, typename F>
inline
difference_type operator-(BinaryIterator<A1, A2, F> const& i1, BinaryIterator<B1, B2, F> const& i2)
{
  DEBUG_PRECONDITION(i1.base1() - i2.base1() == i1.base2() - i2.base2())
                    (i1.base1() - i2.base1())
                    (i1.base2() - i2.base2());
   return i1.base1() - i2.base1();
}

//
// IteratorConstructor's
//

} // namespace LinearAlgebra

#endif
