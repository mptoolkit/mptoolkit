// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/leftrightstack.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(LEFTRIGHTSTACK_H_JHCHIUHFIULSHDLCIWEHLH)
#define LEFTRIGHTSTACK_H_JHCHIUHFIULSHDLCIWEHLH

#include <deque>
#include "quantumnumbers/symmetrylist.h"
#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"

template <typename T>
class LeftRightStack;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, LeftRightStack<T> const& s);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, LeftRightStack<T>& s);

template <typename T>
class LeftRightStack
{
   public:
      typedef T value_type;

      LeftRightStack();

      value_type const& Left() const { return *LeftTop; }
      value_type const& Right() const { return *RightTop; }

      value_type& Left() { return *LeftTop.mutate(); }
      value_type& Right() { return *RightTop.mutate(); }

      void PushLeft(value_type const& L);
      void PushRight(value_type const& R);

      void PopLeft();
      void PopRight();

   private:
      typedef pvalue_handle<value_type> handle_type;
      typedef pvalue_ptr<value_type> ptr_type;
      std::deque<handle_type> LeftStack, RightStack;
      ptr_type LeftTop, RightTop;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, LeftRightStack const& s);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, LeftRightStack& s);
};

template <typename T>
inline
LeftRightStack<T>::LeftRightStack()
{
}

template <typename T>
inline
void LeftRightStack<T>::PushLeft(value_type const& L)
{
   LeftStack.push_back(LeftTop);
   LeftTop = ptr_type(new value_type(L));
}

template <typename T>
inline
void LeftRightStack<T>::PushRight(value_type const& R)
{
   RightStack.push_back(RightTop);
   RightTop = ptr_type(new value_type(R));
}

template <typename T>
inline
void LeftRightStack<T>::PopLeft()
{
   if (LeftStack.empty())
   {
      LeftTop = ptr_type();
   }
   else
   {
      LeftTop = LeftStack.back().load();
      LeftStack.pop_back();
   }
}

template <typename T>
inline
void LeftRightStack<T>::PopRight()
{
   PRECONDITION(!RightTop.is_null());
   if (RightStack.empty())
   {
      RightTop = ptr_type();
   }
   else
   {
      RightTop = RightStack.back().load();
      RightStack.pop_back();
   }
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, LeftRightStack<T> const& s)
{
   return out << s.LeftStack
              << s.RightStack
              << s.LeftTop
              << s.RightTop;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, LeftRightStack<T>& s)
{
   return in >> s.LeftStack
             >> s.RightStack
             >> s.LeftTop
             >> s.RightTop;
}

#endif
