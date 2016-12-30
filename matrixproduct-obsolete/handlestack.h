// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/handlestack.h
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

#if !defined(HANDLESTACK_H_JHCHIUHFIULSHDLCIWEHLH)
#define HANDLESTACK_H_JHCHIUHFIULSHDLCIWEHLH

#include <deque>
#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"

template <typename T>
class HandleStack;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, HandleStack<T> const& s);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, HandleStack<T>& s);

template <typename T>
class HandleStack
{
   public:
      typedef T value_type;

      HandleStack();

      value_type const& top() const { return *Top; }
      value_type& top() { return *Top.mutate(); }

      void push(value_type const& x);

      void pop();

   private:
      typedef pvalue_handle<value_type> handle_type;
      typedef pvalue_ptr<value_type> ptr_type;
      std::deque<handle_type> Stack;
      ptr_type Top;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, HandleStack const& s);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, HandleStack& s);
};

template <typename T>
inline
HandleStack<T>::HandleStack()
{
}

template <typename T>
inline
void HandleStack<T>::push(value_type const& L)
{
   if (!Top.is_null()) Stack.push_back(Top);
   Top = ptr_type(new value_type(L));
}

template <typename T>
inline
void HandleStack<T>::pop()
{
   PRECONDITION(!Top.is_null());
   if (Stack.empty())
   {
      Top = ptr_type();
   }
   else
   {
      Top = Stack.back().load();
      Stack.pop_back();
   }
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, HandleStack<T> const& s)
{
   return out << s.Stack
              << s.Top;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, HandleStack<T>& s)
{
   return in >> s.Stack
             >> s.Top;
}

#endif
