// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/dataops_matrix.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
//
// Overloads in namespace ops:: for matrix operations.
//

#if !defined(DATAOPS_MATRIX_H_HJCDSUI34H89H9H9P8H849HP)
#define DATAOPS_MATRIX_H_HJCDSUI34H89H9H9P8H849HP

#include "matrixiterators.h"
#include "dataops.h"

namespace LinearAlgebra
{

template <typename OuterIter>
class UnwrappedIterator
{
   public:
      typedef OuterIter outer_iterator;
      typedef typename outer_iterator::iterator inner_iterator;
      typedef IteratorTraits<inner_iterator> traits_type;
      typedef typename traits_type::value_type value_type;
      typedef typename traits_type::result_type result_type;

      UnwrappedIterator() {}

      explicit UnwrappedIterator(outer_iterator const& Out)
        : Outer(Out), Inner(Out.begin()), InnerEnd(Out.end()) {}

      UnwrappedIterator& operator++()
        { if (++Inner == InnerEnd) { ++Outer; Inner = Outer.begin(); InnerEnd = Outer.end(); } return *this; }

      UnwrappedIterator& operator--()
        { if (Inner == Outer.begin()) { --Outer; Inner = InnerEnd = Outer.end(); } --Inner; return *this; }

      UnwrappedIterator& operator+=(difference_type n)
        {
           if (n > 0)
           {
              int Remain = InnerEnd - Inner;
              if (n < Remain)
              {
                 Inner += n;
                 return *this;
              }
              n -= Remain;
              ++Outer;
           }
           else
           {
              int Remain = Outer.begin() - Inner;  // negative
              if (n > Remain)
              {
                 Inner += n;
                 return *this;
              }
              n -= Remain;
           }

           numerics::div_result<difference_type> Res = numerics::divp(n, Outer.size());
           Outer += Res.quot;
           Inner = Outer.begin() + Res.rem;
           InnerEnd = Outer.end();
           return *this;
        }

      UnwrappedIterator& operator-=(difference_type n) { return this->operator+=(-n); }

     result_type operator*() const { return *Inner; }

     result_type operator[](difference_type n) const { UnwrappedIterator I(*this); I += n; return *I; }

     inner_iterator const& inner() const { return Inner; }
     outer_iterator const& outer() const { return Outer; }

   private:
      outer_iterator Outer;
      inner_iterator Inner, InnerEnd;
};

template <typename Iter1, typename Iter2>
bool operator==(UnwrappedIterator<Iter1> const& i1, UnwrappedIterator<Iter2> const& i2)
{
   return i1.outer() == i2.outer() && i1.inner() == i2.inner();
}

template <typename Iter1, typename Iter2>
bool operator!=(UnwrappedIterator<Iter1> const& i1, UnwrappedIterator<Iter2> const& i2)
{
   return i1.outer() != i2.outer() || i2.inner() != i2.inner();
}

// If ther Iter types are MatrixOuterIter, we know that it is a 2D view of a single
// container, so it is valid to compare inner iterators directly.
// However we cannot avoid checking the outer iterator because
// there are some circumstances where the inner iterator may be the same
// but the outer is different.  For example, the end() of a transposed iterator
// points to the (1,0) element.  So the best we can do is optimize
// the order for short-circuit evaluation.
template <typename Iter1, typename Iter2, typename F1, typename F2>
bool operator==(UnwrappedIterator<MatrixOuterIterator<Iter1, F1, F2> > const& i1,
                UnwrappedIterator<MatrixOuterIterator<Iter2, F1, F2> > const& i2)
{
   return i1.inner() == i2.inner() && i1.outer() == i2.outer();
}

template <typename Iter1, typename Iter2, typename F1, typename F2>
bool operator!=(UnwrappedIterator<MatrixOuterIterator<Iter1, F1, F2> > const& i1,
                UnwrappedIterator<MatrixOuterIterator<Iter2, F1, F2> > const& i2)
{
   return i1.inner() != i2.inner() || i1.outer() != i2.outer();
}

template <typename OuterIterator>
struct IteratorTraits<UnwrappedIterator<OuterIterator> >
  : public IteratorTraits<typename OuterIterator::inner_iterator>
{ };


} // namespace LinearAlgebra


namespace ops
{

#if !defined(DISABLE_FANCY_OVERLOADS)


template <typename Iter>
struct UnwrapIterator
{
   typedef UnwrappedIterator<Iter> result_type;
   static bool is_unwrap(Iter const& first) { return false; }
   static result_type unwrap(Iter const& x) { PANIC("Attempting to unwrap an iterator that cannot be unwrapped"); return result_type(); }
};

template <typename Iter>
struct UnwrapIterator<MatrixOuterIterator<Iter, Slice, Slice> >
{
   typedef StrideIterator<Iter> result_type;

   static bool is_unwrap(MatrixOuterIterator<Iter, Slice, Slice> const& first)
   {
      int Stride1 = first.func1().stride();
      int Stride2 = first.func2().stride();
      int Size2 = first.func2().size();

      return Stride1 == Size2 * Stride2;
   }

   static result_type unwrap(MatrixOuterIterator<Iter, Slice, Slice> const& first)
   {
      DEBUG_PRECONDITION((UnwrapIterator<MatrixOuterIterator<Iter, Slice, Slice> >::is_unwrap(first)));
      int Stride1 = first.func1().stride();
      int Stride2 = first.func2().stride();
      int Start = first.func1().start() + first.func2().start();
      int FirstOffset = first.index1() * Stride1 + Start;
      return StrideIterator<Iter>(first.base_begin() + FirstOffset, Stride2);
   }
};

template <typename Iter>
struct UnwrapIterator<MatrixOuterIterator<Iter, Slice, Range> >
{
   typedef StrideIterator<Iter> result_type;

   static bool is_unwrap(MatrixOuterIterator<Iter, Slice, Range> const& first)
   {
      int Stride1 = first.func1().stride();
      int Stride2 = first.func2().stride();
      int Size2 = first.func2().size();

      return Stride1 == Size2 * Stride2;
   }

   static result_type unwrap(MatrixOuterIterator<Iter, Slice, Range> const& first)
   {
      DEBUG_PRECONDITION((UnwrapIterator<MatrixOuterIterator<Iter, Slice, Range> >::is_unwrap(first)));
      int Start = first.func1().start() + first.func2().first();
      int Stride1 = first.func1().stride();
      int FirstOffset = first.index1() * Stride1 + Start;
      return StrideIterator<Iter>(first.base_begin() + FirstOffset, 1);
   }
};

// we can never unwrap a (range,slice) ???
template <typename Iter>
struct UnwrapIterator<MatrixOuterIterator<Iter, Range, Slice> >
{
   typedef StrideIterator<Iter> result_type;

   static bool is_unwrap(MatrixOuterIterator<Iter, Range, Slice> const& first)
   {
      int Stride1 = first.func1().stride();
      int Stride2 = first.func2().stride();
      int Size2 = first.func2().size();

      return Stride1 == Size2 * Stride2;
   }

   static result_type unwrap(MatrixOuterIterator<Iter, Range, Slice> const& first)
   {
      DEBUG_PRECONDITION((UnwrapIterator<MatrixOuterIterator<Iter, Range, Slice> >::is_unwrap(first)));
      int Start = first.func1().first() + first.func2().start();
      int FirstOffset = first.index1() + Start;
      int Stride2 = first.func2().stride();
      return StrideIterator<Iter>(first.base_begin() + FirstOffset, Stride2);
   }
};

template <typename Iter, typename F>
typename F::result_type
unwrap_dispatch_unordered(MatrixOuterIterator<Iter, Slice, Slice> const& first,
                          MatrixOuterIterator<Iter, Slice, Slice> const& last,
                          F f)
{
   DEBUG_PRECONDITION(first.base_begin() == last.base_begin());
   DEBUG_PRECONDITION(first.func1() == last.func1());
   DEBUG_PRECONDITION(first.func2() == last.func2());

   int Stride1 = first.func1().stride();
   int Stride2 = first.func2().stride();

   int Size1 = last-first;
   int Size2 = first.func2().size();

   int Start = first.func1().start() + first.func2().start();

   // see if we can linearize the iterator.  This can be done if
   // we can find a fixed stride that will iterate over all elements
   if (Stride1 == Size2 * Stride2)
   {
      int FirstOffset = first.index1() * Stride1 + Start;
      int LastOffset = last.index1() * Stride1 + Start;
      return f(StrideIterator<Iter>(first.base_begin() + FirstOffset, Stride2),
               StrideIterator<Iter>(first.base_begin() + LastOffset, Stride2));
   }

   if (Stride2 == Size1 * Stride1)
   {
      int FirstOffset = first.index1() * Stride1 + Start;
      int LastOffset = FirstOffset + Size2 * Stride2;
      return f(StrideIterator<Iter>(first.base_begin() + FirstOffset, Stride1),
               StrideIterator<Iter>(first.base_begin() + LastOffset, Stride1));
   }

   // else linearization is not possible
   typedef UnwrappedIterator<MatrixOuterIterator<Iter, Slice, Slice> > Unwrapped;
   return f(Unwrapped(first), Unwrapped(last));
}

struct FastNorm2Sq
{
   typedef double result_type;
   template <typename Iter>
   double operator()(Iter const& first, Iter const& last) const { return fast_norm_2_sq(first, last); }
};

template <typename Iter>
inline
double fast_norm_2_sq(MatrixOuterIterator<Iter, Slice, Slice> first,
                      MatrixOuterIterator<Iter, Slice, Slice> last)
{
   return unwrap_dispatch_unordered(first, last, FastNorm2Sq());
}

struct FastFill
{
   typedef void result_type;
   template <typename Iter>
   void operator()(Iter const& first, Iter const& last) const { return fast_fill(first, last); }
};

template <typename Iter>
inline
void fast_fill(MatrixOuterIterator<Iter, Slice, Slice> first,
               MatrixOuterIterator<Iter, Slice, Slice> last)
{
   return unwrap_dispatch_unordered(first, last, FastFill());
}

template <typename Iter, typename F1, typename F2, typename Scalar>
inline
void fast_fill(MatrixOuterIterator<Iter, F1, F2> first,
               MatrixOuterIterator<Iter, F1, F2> last, Scalar const& s)
{
   typedef UnwrapIterator<MatrixOuterIterator<Iter, F1, F2> > UnwrapIter;
   typedef UnwrappedIterator<MatrixOuterIterator<Iter, F1, F2> > UIterType;

   if (UnwrapIter::is_unwrap(first))
   {
      fast_fill(UnwrapIter::unwrap(first), UnwrapIter::unwrap(last), s);
   }
   else
   {
      fast_fill(UIterType(first), UIterType(last), s);
   }
}

template <typename Iter, typename F1, typename F2, typename Scalar>
inline
void fast_multiply(MatrixOuterIterator<Iter, F1, F2> first,
                   MatrixOuterIterator<Iter, F1, F2> last, Scalar const& s)
{
   typedef UnwrapIterator<MatrixOuterIterator<Iter, F1, F2> > UnwrapIter;
   typedef UnwrappedIterator<MatrixOuterIterator<Iter, F1, F2> > UIterType;

   if (UnwrapIter::is_unwrap(first))
   {
      fast_multiply(UnwrapIter::unwrap(first), UnwrapIter::unwrap(last), s);
   }
   else
   {
      fast_multiply(UIterType(first), UIterType(last), s);
   }
}

template <typename Iter1, typename F11, typename F12, typename Iter2, typename F21, typename F22>
void fast_add(MatrixOuterIterator<Iter1, F11, F12> first, MatrixOuterIterator<Iter1, F11, F12> last,
              MatrixOuterIterator<Iter2, F21, F22> dest)
{
   typedef UnwrapIterator<MatrixOuterIterator<Iter1, F11, F12> > UnwrapIter1;
   typedef UnwrapIterator<MatrixOuterIterator<Iter2, F21, F22> > UnwrapIter2;

   typedef UnwrappedIterator<MatrixOuterIterator<Iter1, F11, F12> > UIterType1;
   typedef UnwrappedIterator<MatrixOuterIterator<Iter2, F21, F22> > UIterType2;

   if (UnwrapIter1::is_unwrap(first))
   {
      if (UnwrapIter2::is_unwrap(dest))
      {
         fast_add(UnwrapIter1::unwrap(first), UnwrapIter1::unwrap(last), UnwrapIter2::unwrap(dest));
      }
      else
      {
         fast_add(UnwrapIter1::unwrap(first), UnwrapIter1::unwrap(last), UIterType2(dest));
      }
   }
   else
   {
      if (UnwrapIter2::is_unwrap(dest))
      {
         fast_add(UIterType1(first), UIterType1(last), UnwrapIter2::unwrap(dest));
      }
      else
      {
         fast_add(UIterType1(first), UIterType1(last), UIterType2(dest));
      }
   }
}

template <typename Iter1, typename F11, typename F12, typename Iter2, typename F21, typename F22>
void fast_subtract(MatrixOuterIterator<Iter1, F11, F12> first, MatrixOuterIterator<Iter1, F11, F12> last,
              MatrixOuterIterator<Iter2, F21, F22> dest)
{
   typedef UnwrapIterator<MatrixOuterIterator<Iter1, F11, F12> > UnwrapIter1;
   typedef UnwrapIterator<MatrixOuterIterator<Iter2, F21, F22> > UnwrapIter2;

   typedef UnwrappedIterator<MatrixOuterIterator<Iter1, F11, F12> > UIterType1;
   typedef UnwrappedIterator<MatrixOuterIterator<Iter2, F21, F22> > UIterType2;

   if (UnwrapIter1::is_unwrap(first))
   {
      if (UnwrapIter2::is_unwrap(dest))
      {
         fast_subtract(UnwrapIter1::unwrap(first), UnwrapIter1::unwrap(last), UnwrapIter2::unwrap(dest));
      }
      else
      {
         fast_subtract(UnwrapIter1::unwrap(first), UnwrapIter1::unwrap(last), UIterType2(dest));
      }
   }
   else
   {
      if (UnwrapIter2::is_unwrap(dest))
      {
         fast_subtract(UIterType1(first), UIterType1(last), UnwrapIter2::unwrap(dest));
      }
      else
      {
         fast_subtract(UIterType1(first), UIterType1(last), UIterType2(dest));
      }
   }
}

template <typename Iter1, typename F11, typename F12, typename Iter2, typename F21, typename F22>
void fast_copy(MatrixOuterIterator<Iter1, F11, F12> first, MatrixOuterIterator<Iter1, F11, F12> last,
               MatrixOuterIterator<Iter2, F21, F22> dest)
{
   typedef UnwrapIterator<MatrixOuterIterator<Iter1, F11, F12> > UnwrapIter1;
   typedef UnwrapIterator<MatrixOuterIterator<Iter2, F21, F22> > UnwrapIter2;

   typedef UnwrappedIterator<MatrixOuterIterator<Iter1, F11, F12> > UIterType1;
   typedef UnwrappedIterator<MatrixOuterIterator<Iter2, F21, F22> > UIterType2;

   if (UnwrapIter1::is_unwrap(first))
   {
      if (UnwrapIter2::is_unwrap(dest))
      {
         fast_copy(UnwrapIter1::unwrap(first), UnwrapIter1::unwrap(last), UnwrapIter2::unwrap(dest));
      }
      else
      {
         fast_copy(UnwrapIter1::unwrap(first), UnwrapIter1::unwrap(last), UIterType2(dest));
      }
   }
   else
   {
      if (UnwrapIter2::is_unwrap(dest))
      {
         fast_copy(UIterType1(first), UIterType1(last), UnwrapIter2::unwrap(dest));
      }
      else
      {
         fast_copy(UIterType1(first), UIterType1(last), UIterType2(dest));
      }
   }
}

#endif

} // namespace ops

#endif
