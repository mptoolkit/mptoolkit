// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/runlengthcompressed.h
//
// Copyright (C) 2006-2016 Ian McCulloch <ian@qusim.net>
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
// A run-length-encoded compressed vector.
// Created 2006-11-22 Ian McCulloch
//
// This tries as much as possible to keep the representation canonical,
// that is, repeats of length 1 are collapsed down to its component part,
// and nested arrays are flattened.  This can be guaranteed, for all
// operations that act on a run_length_compressed object.
// If the user modifies the constituent run_length_repeat or run_length_array
// objects, then it is up to the user to canonicalize the representation.
//
// The const_iterator is a standard-conforming bidirectional iterator.
// To modify a run_length_compressed object, use the visitor pattern
// with a function object that implements operator() for T itself,
// run_length_repeat<T> and run_length_array<T>.  The latter two types
// act as containers of run_length_compressed<T> objects.
// Note that for run_length_array and run_length_repeat, the
// size() function returns the physical size of the container,
// which is NOT the same as the logical size returned by
// run_length_compressed::size().  The repeat and array objects have
// instead a logical_size() function that takes into account the
// complete size of nested objects.
//

#if !defined(MPTOOLKIT_COMMON_RUNLENGTHCOMPRESSED_H)
#define MPTOOLKIT_COMMON_RUNLENGTHCOMPRESSED_H

#if defined(HAVE_CONFIG_H)  // we should include config.h before boost headers
#include "config.h"
#endif
#include <boost/variant.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>
#include <tuple>
#include <deque>
#include <algorithm>
#include <functional>
#include "common/trace.h"
#if defined(USE_PSTREAM)
#include "pstream/variant.h"
#endif

template <typename T>
class run_length_repeat;

template <typename T>
class run_length_array;

template <typename T>
class run_length_compressed
{
   public:
      typedef T value_type;
      class const_iterator;

      run_length_compressed();
      run_length_compressed(run_length_compressed const& x) : Data(x.Data) {}
      run_length_compressed(run_length_compressed&& x) : Data(std::move(x.Data)) {}

      explicit run_length_compressed(T const& x);
      run_length_compressed(int Size, T const& x);
      run_length_compressed(run_length_repeat<T> const& x);
      run_length_compressed(run_length_array<T> const& x);
      run_length_compressed(int Size, run_length_compressed const& x);

      run_length_compressed(T const& x1, T const& x2);
      run_length_compressed(T const& x1, T const& x2, T const& x3);
      run_length_compressed(T const& x1, T const& x2, T const& x3, T const& x4);

      // construction from an existing container
      template <typename FwdIter>
      run_length_compressed(FwdIter start, FwdIter finish);

      bool empty() const;

      int size() const;

      // returns the number of actual T objects stored in the container.
      int leaf_count() const;

      // returns the number of non-terminal nodes in the container.
      // The compression ratio is roughly size() / (leaf_count() + node_count())
      int node_count() const;

      T const& front() const;
      T const& back() const;

      // returns the n'th element of the container, equivalent to
      // *std::advance(this->begin(), n);
      T const& find(int n) const;

      void push_front(T const& x);
      void push_front(run_length_repeat<T> const& x);
      void push_front(run_length_array<T> const& x);
      void push_front(run_length_compressed<T> const& x);

      void push_back(T const& x);
      void push_back(run_length_repeat<T> const& x);
      void push_back(run_length_array<T> const& x);
      void push_back(run_length_compressed<T> const& x);

      template <typename FwdIter>
      void append(FwdIter start, FwdIter finish);

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v) const
      {
         return boost::apply_visitor(v, Data);
      }

      // should the Visitor be non-const here?  I don't see why not,
      // but check boost docs!
      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v)
      {
         return boost::apply_visitor(v, Data);
      }

      const_iterator begin() const;
      const_iterator end() const;

      typedef boost::variant<T,
         boost::recursive_wrapper<run_length_repeat<T> >,
         boost::recursive_wrapper<run_length_array<T> > > data_type;

      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

   private:
      data_type Data;
};

// Replicates x by RepeatCount
template <typename T>
run_length_compressed<T>
repeat(run_length_compressed<T> const& x, int RepeatCount);

// Joins two items together
template <typename T>
run_length_compressed<T>
join(run_length_compressed<T> const& x, run_length_compressed<T> const& y);

// Splits a run_length_compressed<T> container into two
// parts, having sizes s and size()-s.
template <typename T>
std::pair<run_length_compressed<T>, run_length_compressed<T> >
split(run_length_compressed<T> const& x, int Loc);

// Variant of split, breaks a run into 3 components,
// the left array, the element itself, and the right array.
// This is currently implemented in terms of split(), but
// it might be that this is the more useful primitive.
template <typename T>
std::tuple<run_length_compressed<T>, T, run_length_compressed<T> >
split_lcr(run_length_compressed<T> const& x, int Loc);

// Forwards to boost::apply_visitor(v, x.data())
template <typename T, typename Visitor>
typename Visitor::result_type
apply_visitor(Visitor const& v, run_length_compressed<T> const& x);

template <typename T, typename Visitor>
typename Visitor::result_type
apply_visitor(Visitor const& v, run_length_compressed<T>& x);

template <typename T, typename Visitor>
typename Visitor::result_type
apply_visitor(Visitor const& v,
              run_length_compressed<T> const& x,
              run_length_compressed<T> const& y);

//
// run_length_repeat
//

template <typename T>
class run_length_repeat
{
   public:
      typedef run_length_compressed<T> value_type;

      struct const_iterator
      {
         typedef std::bidirectional_iterator_tag iterator_category;
         typedef run_length_compressed<T> value_type;
         typedef int difference_type;
         typedef value_type const& reference;
         typedef value_type const* pointer;
         typedef int size_type;

         const_iterator() {}

         const_iterator(run_length_compressed<T> const* Value_, int Loc_)
            : Value(Value_), Loc(Loc_) { DEBUG_CHECK(Loc >= 0); }

         reference operator*() const { return *Value; }

         pointer operator->() const { return Value; }

         bool operator==(const_iterator i) const
         {
            DEBUG_CHECK(Value == i.Value);
            return Loc == i.Loc;
         }

         bool operator!=(const_iterator i) const
         {
            DEBUG_CHECK(Value == i.Value);
            return Loc != i.Loc;
         }

         const_iterator& operator++()
         {
            ++Loc;
            return *this;
         }

         const_iterator operator++(int)
         {
            return const_iterator(Value, Loc++);
         }

         const_iterator& operator--()
         {
            --Loc;
            DEBUG_CHECK(Loc >= 0)(Loc);
            return *this;
         }

         const_iterator operator--(int)
         {
            return const_iterator(Value, Loc--);
         }

         run_length_compressed<T> const* Value;
         int Loc;
      };

      run_length_repeat() : Count(0) {}
      run_length_repeat(int RepeatCount, run_length_compressed<T> const& x);

      const_iterator begin() const { return const_iterator(&Value, 0); }
      const_iterator end() const { return const_iterator(&Value, Count); }


      int size() const { return Count; }

      bool empty() const { return Count == 0; }

      int logical_size() const { return Count * Value.size(); }

      value_type const& nested() const { return Value; }

      run_length_compressed<T> const& front() const;
      run_length_compressed<T> const& back() const;

      value_type& nested() { return Value; }

   //   private:
      value_type Value;
      int Count;
};

//
// run_length_array
//

template <typename T>
class run_length_array
{
   public:
      typedef run_length_compressed<T> value_type;
      typedef std::deque<value_type> data_type;

      typedef typename data_type::iterator iterator;
      typedef typename data_type::const_iterator const_iterator;
      typedef typename data_type::reference reference;
      typedef typename data_type::const_reference const_reference;

      run_length_array() {}
      run_length_array(T const& x);
      run_length_array(run_length_compressed<T> const& x);
      run_length_array(run_length_repeat<T> const& x);

      template <typename FwdIter>
      run_length_array(FwdIter first, FwdIter last);

      bool empty() const { return Data.empty(); }

      int size() const { return Data.size(); }

      int logical_size() const;

      run_length_compressed<T> const& front() const;
      run_length_compressed<T> const& back() const;

      void push_front(T const& x);
      void push_front(run_length_repeat<T> const& x);
      void push_front(run_length_array<T> const& x);
      void push_front(run_length_compressed<T> const& x);

      void push_back(T const& x);
      void push_back(run_length_repeat<T> const& x);
      void push_back(run_length_array<T> const& x);
      void push_back(run_length_compressed<T> const& x);

      template <typename FwdIter>
      void append(FwdIter start, FwdIter finish);

      iterator begin() { return Data.begin(); }
      iterator end() { return Data.end(); }

      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      // direct access to the representation, internal use only
      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

   private:
      data_type Data;
};

#if defined(USE_PSTREAM)
template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, run_length_compressed<T> const& x);
template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, run_length_repeat<T> const& x);
template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, run_length_array<T> const& x);

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, run_length_compressed<T>& x);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, run_length_repeat<T>& x);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, run_length_array<T>& x);
#endif

#include "runlengthcompressed.cc"

#endif
