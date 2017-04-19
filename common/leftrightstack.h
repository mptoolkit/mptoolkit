// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/leftrightstack.h
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
//
// A dual stack based around a stack of pvalue_handle,
// with direct access to the stack top.

#if !defined(MPTOOLKIT_COMMON_LEFTRIGHTSTACK_H)
#define MPTOOLKIT_COMMON_LEFTRIGHTSTACK_H

#include <deque>
#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"

template <typename T>
class left_right_stack;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, left_right_stack<T> const& s);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, left_right_stack<T>& s);

template <typename T>
class left_right_stack
{
   public:
      typedef T value_type;

      left_right_stack();

      value_type const& left() const { return *leftTop; }
      value_type const& right() const { return *rightTop; }

      value_type& left() { return *leftTop.mutate(); }
      value_type& right() { return *rightTop.mutate(); }

      void push_left(value_type const& L);
      void push_right(value_type const& R);

      void pop_left();
      void pop_right();

      int size_left() const;
      int size_right() const;

   private:
      typedef pvalue_handle<value_type> handle_type;
      typedef pvalue_ptr<value_type> ptr_type;
      std::deque<handle_type> leftStack, rightStack;
      ptr_type leftTop, rightTop;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, left_right_stack const& s);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, left_right_stack& s);
};

template <typename T>
inline
left_right_stack<T>::left_right_stack()
{
}

template <typename T>
inline
void left_right_stack<T>::push_left(value_type const& L)
{
   if (!leftTop.is_null()) leftStack.push_back(leftTop);
   leftTop = ptr_type(new value_type(L));
}

template <typename T>
inline
void left_right_stack<T>::push_right(value_type const& R)
{
   if (!rightTop.is_null()) rightStack.push_back(rightTop);
   rightTop = ptr_type(new value_type(R));
}

template <typename T>
inline
void left_right_stack<T>::pop_left()
{
   PRECONDITION(!leftTop.is_null());
   if (leftStack.empty())
   {
      leftTop = ptr_type();
   }
   else
   {
      leftTop = leftStack.back().load();
      leftStack.pop_back();
   }
}

template <typename T>
inline
void left_right_stack<T>::pop_right()
{
   PRECONDITION(!rightTop.is_null());
   if (rightStack.empty())
   {
      rightTop = ptr_type();
   }
   else
   {
      rightTop = rightStack.back().load();
      rightStack.pop_back();
   }
}

template <typename T>
inline
int left_right_stack<T>::size_left() const
{
   if (leftStack.empty())
   {
      return leftTop.is_null() ? 0 : 1;
   }
   return leftStack.size() + 1;  // +1 for leftTop
}

template <typename T>
inline
int left_right_stack<T>::size_right() const
{
   if (rightStack.empty())
   {
      return rightTop.is_null() ? 0 : 1;
   }
   return rightStack.size() + 1;  // +1 for rightTop
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, left_right_stack<T> const& s)
{
   return out << s.leftStack
              << s.rightStack
              << s.leftTop
              << s.rightTop;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, left_right_stack<T>& s)
{
   return in >> s.leftStack
             >> s.rightStack
             >> s.leftTop
             >> s.rightTop;
}

#endif
