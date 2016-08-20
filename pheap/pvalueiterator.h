// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pvalueiterator.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
//
// An iterator for a pvalue_handle that automatically dereferences.
// Warning: This iterator has non-standard semantics due to the
// copy-on-write properties of the pvalue_handle.  Multiple const_iterators
// work fine, but if an object is modified via a mutable iterator,
// while a second iterator is pointing at the same object (whether const or non-const),
// that second iterator will not see the changes.
// This could be avoided only at the expense of keeping a list of all active
// mutable iterators and where they are pointing to.

#if !defined(PVALUEITERATOR_H_SDUI3489SD)
#define PVALUEITERATOR_H_SDUI3489SD

#include "pvalueptr.h"
#include <iterator>

// metafunction to compute the iterator category.  Random access iterators
// are not supported.
template <typename Category>
struct pvalue_iter_category
{
   typedef Category type;
};

template <>
struct pvalue_iter_category<std::random_access_iterator_tag>
{
   typedef std::bidirectional_iterator_tag type;
};

template <typename BaseIter>
class pvalue_handle_iterator
{
   public:
      typedef BaseIter base_iterator;
      typedef typename std::iterator_traits<BaseIter>::value_type handle_type;
      typedef typename pvalue_iter_category<typename std::iterator_traits<BaseIter>::iterator_category>::type iterator_category;
   //typedef typename std::iterator_traits<BaseIter>::size_type size_type;
      typedef typename std::iterator_traits<BaseIter>::difference_type difference_type;
      typedef typename handle_type::ptr_type ptr_type;
      typedef typename handle_type::element_type value_type;
      typedef value_type& reference;
      typedef value_type* pointer;

      pvalue_handle_iterator() : rPtr(NULL) {}

      explicit pvalue_handle_iterator(base_iterator const& b)
         : Iter(b), rPtr(NULL) {}

      // This is a tricky case, copying an iterator and simultaneously writing
      // to both, will not work properly
      pvalue_handle_iterator(pvalue_handle_iterator const& other)
         : Iter(other.Iter), rPtr(NULL) {}

      pvalue_handle_iterator& operator=(pvalue_handle_iterator const& other)
      {
         Iter = other.Iter;
         if (this != &other)
         {
            Ptr.release();
            rPtr = NULL;
         }
         return *this;
      }

      reference operator*() const
      {
         if (Ptr.is_null() && !Iter->is_null())
         {
            Ptr = Iter->load();
            rPtr = Ptr.mutate();
            *Iter = Ptr;         // re-assign after the copy-on-write
         }
         CHECK(rPtr != NULL)(rPtr);
         return *rPtr;
      }

      pointer operator->() const
      {
         if (Ptr.is_null() && !Iter->is_null())
         {
            Ptr = Iter->load();
            rPtr = Ptr.mutate();
            *Iter = Ptr;         // re-assign after the copy-on-write
         }
         CHECK(rPtr != NULL)(rPtr);
         return rPtr;
      }

      pvalue_handle_iterator& operator++()
      {
         Ptr.release();
         rPtr = NULL;
         ++Iter;
         return *this;
      }

      pvalue_handle_iterator operator++(int)
      {
         pvalue_handle_iterator Temp(*this);
         Ptr.release();
         rPtr = NULL;
         ++Iter;
         return Temp;
      }

      pvalue_handle_iterator operator+=(int n)
      {
         Ptr.release();
         rPtr = NULL;
         Iter += n;
         return *this;
      }

      pvalue_handle_iterator& operator--()
      {
         Ptr.release();
         rPtr = NULL;
         --Iter;
         return *this;
      }

      pvalue_handle_iterator operator--(int)
      {
         pvalue_handle_iterator Temp(*this);
         Ptr.release();
         rPtr = NULL;
         --Iter;
         return Temp;
      }

      pvalue_handle_iterator operator-=(int n)
      {
         Ptr.release();
         rPtr = NULL;
         Iter -= n;
         return *this;
      }

      base_iterator const& base() const { return Iter; }

   private:
      base_iterator    Iter;
      mutable ptr_type Ptr;
      mutable pointer  rPtr;
};

template <typename BaseIter>
inline
bool operator==(pvalue_handle_iterator<BaseIter> const& a,
                pvalue_handle_iterator<BaseIter> const& b)
{
   return a.base() == b.base();
}

template <typename BaseIter>
inline
bool operator!=(pvalue_handle_iterator<BaseIter> const& a,
                pvalue_handle_iterator<BaseIter> const& b)
{
   return a.base() != b.base();
}

template <typename BaseIter>
inline
pvalue_handle_iterator<BaseIter>
operator+(pvalue_handle_iterator<BaseIter> const& a, int n)
{
   pvalue_handle_iterator<BaseIter> Result(a);
   Result += n;
   return Result;
}

template <typename BaseIter>
inline
pvalue_handle_iterator<BaseIter>
operator-(pvalue_handle_iterator<BaseIter> const& a, int n)
{
   pvalue_handle_iterator<BaseIter> Result(a);
   Result -= n;
   return Result;
}

template <typename BaseIter>
class const_pvalue_handle_iterator
{
   public:
      typedef BaseIter base_iterator;
      typedef typename std::iterator_traits<BaseIter>::value_type handle_type;
      typedef typename pvalue_iter_category<typename std::iterator_traits<BaseIter>::iterator_category>::type iterator_category;
   //typedef typename std::iterator_traits<BaseIter>::size_type size_type;
      typedef typename std::iterator_traits<BaseIter>::difference_type difference_type;
      typedef typename handle_type::ptr_type ptr_type;
      typedef typename handle_type::element_type value_type;
      typedef value_type const& reference;
      typedef value_type const* pointer;

      const_pvalue_handle_iterator() {}

      explicit const_pvalue_handle_iterator(base_iterator const& b)
         : Iter(b) {}

      template <typename Other>
      const_pvalue_handle_iterator(pvalue_handle_iterator<Other> const& I)
         : Iter(I.base()) {}

      reference operator*() const
      {
         if (Ptr.is_null() && !Iter->is_null())
            Ptr = Iter->load();
         return *Ptr;
      }

      pointer operator->() const
      {
         if (Ptr.is_null() && !Iter->is_null())
            Ptr = Iter->load();
         return Ptr.get();
      }

      const_pvalue_handle_iterator& operator++()
      {
         Ptr.release();
         ++Iter;
         return *this;
      }

      const_pvalue_handle_iterator operator++(int)
      {
         const_pvalue_handle_iterator Temp(*this);
         Ptr.release();
         ++Iter;
         return Temp;
      }

      const_pvalue_handle_iterator operator+=(int n)
      {
         Ptr.release();
         Iter += n;
         return *this;
      }

      const_pvalue_handle_iterator& operator--()
      {
         Ptr.release();
         --Iter;
         return *this;
      }

      const_pvalue_handle_iterator operator--(int)
      {
         const_pvalue_handle_iterator Temp(*this);
         Ptr.release();
         --Iter;
         return Temp;
      }

      const_pvalue_handle_iterator operator-=(int n)
      {
         Ptr.release();
         Iter -= n;
         return *this;
      }

      base_iterator const& base() const { return Iter; }

   private:
      BaseIter Iter;
      mutable ptr_type Ptr;
};

template <typename BaseIter>
inline
bool operator==(const_pvalue_handle_iterator<BaseIter> const& a,
                const_pvalue_handle_iterator<BaseIter> const& b)
{
   return a.base() == b.base();
}

template <typename BaseIter>
inline
bool operator!=(const_pvalue_handle_iterator<BaseIter> const& a,
                const_pvalue_handle_iterator<BaseIter> const& b)
{
   return a.base() != b.base();
}

// mixed const/non-const comparisons

template <typename BaseIter1, typename BaseIter2>
inline
bool operator==(pvalue_handle_iterator<BaseIter1> const& a,
                const_pvalue_handle_iterator<BaseIter2> const& b)
{
   return a.base() == b.base();
}

template <typename BaseIter1, typename BaseIter2>
inline
bool operator!=(pvalue_handle_iterator<BaseIter1> const& a,
                const_pvalue_handle_iterator<BaseIter2> const& b)
{
   return a.base() != b.base();
}

template <typename BaseIter1, typename BaseIter2>
inline
bool operator==(const_pvalue_handle_iterator<BaseIter1> const& a,
                pvalue_handle_iterator<BaseIter2> const& b)
{
   return a.base() == b.base();
}

template <typename BaseIter1, typename BaseIter2>
inline
bool operator!=(const_pvalue_handle_iterator<BaseIter1> const& a,
                pvalue_handle_iterator<BaseIter2> const& b)
{
   return a.base() != b.base();
}

template <typename BaseIter>
inline
const_pvalue_handle_iterator<BaseIter>
operator+(const_pvalue_handle_iterator<BaseIter> const& a, int n)
{
   const_pvalue_handle_iterator<BaseIter> Result(a);
   Result += n;
   return Result;
}

template <typename BaseIter>
inline
const_pvalue_handle_iterator<BaseIter>
operator-(const_pvalue_handle_iterator<BaseIter> const& a, int n)
{
   const_pvalue_handle_iterator<BaseIter> Result(a);
   Result -= n;
   return Result;
}

#endif
