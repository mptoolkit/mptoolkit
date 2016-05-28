// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/operatorstack.h
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

#if !defined(OPERATORSTACK_H_JHCHIUHFIULSHDLCIWEHLH)
#define OPERATORSTACK_H_JHCHIUHFIULSHDLCIWEHLH

#include <deque>
#include "quantumnumbers/symmetrylist.h"
#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"

template <typename T>
class OperatorStack;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, OperatorStack<T> const& s);
template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, OperatorStack<T>& s);

template <typename T>
class OperatorStack
{
   public:
      typedef T value_type;

      OperatorStack();

      SymmetryList GetSymmetryList() const;

      value_type const& Left() const { return *LeftTop; }
      value_type const& Right() const { return *RightTop; }

      value_type& Left() { return *LeftTop.mutate(); }
      value_type& Right() { return *RightTop.mutate(); }

      void PushLeft(value_type const& L);
      void PushRight(value_type const& R);

      void PopLeft();
      void PopRight();

      int LeftSize() const;
      int RightSize() const;

   private:
      typedef pvalue_handle<value_type> handle_type;
      typedef pvalue_ptr<value_type> ptr_type;
      std::deque<handle_type> LeftStack, RightStack;
      ptr_type LeftTop, RightTop;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, OperatorStack const& s);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, OperatorStack& s);
};

template <typename T>
inline
OperatorStack<T>::OperatorStack()
{
}

template <typename T>
inline
SymmetryList OperatorStack<T>::GetSymmetryList() const
{
   if (this->Left().is_null())
      return this->Right().GetSymmetryList();
   else
      return this->Left().GetSymmetryList();
}

template <typename T>
inline
void OperatorStack<T>::PushLeft(value_type const& L)
{
   PRECONDITION(LeftTop.is_null() || LeftTop->GetSymmetryList() == L.GetSymmetryList());
   PRECONDITION(RightTop.is_null() || RightTop->GetSymmetryList() == L.GetSymmetryList());
   if (!LeftTop.is_null()) LeftStack.push_back(LeftTop);
   LeftTop = ptr_type(new value_type(L));
}

template <typename T>
inline
void OperatorStack<T>::PushRight(value_type const& R)
{
   PRECONDITION(LeftTop.is_null() || LeftTop->GetSymmetryList() == R.GetSymmetryList());
   PRECONDITION(RightTop.is_null() || RightTop->GetSymmetryList() == R.GetSymmetryList());
   if (!RightTop.is_null()) RightStack.push_back(RightTop);
   RightTop = ptr_type(new value_type(R));
}

template <typename T>
inline
void OperatorStack<T>::PopLeft()
{
   PRECONDITION(!LeftTop.is_null());
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
void OperatorStack<T>::PopRight()
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
inline
int OperatorStack<T>::LeftSize() const
{
   if (LeftStack.empty())
   {
      return LeftTop.is_null() ? 0 : 1;
   }
   return LeftStack.size() + 1;  // +1 for LeftTop
}

template <typename T>
inline
int OperatorStack<T>::RightSize() const
{
   if (RightStack.empty())
   {
      return RightTop.is_null() ? 0 : 1;
   }
   return RightStack.size() + 1;  // +1 for RightTop
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, OperatorStack<T> const& s)
{
   return out << s.LeftStack
              << s.RightStack
              << s.LeftTop
              << s.RightTop;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, OperatorStack<T>& s)
{
   return in >> s.LeftStack
             >> s.RightStack
             >> s.LeftTop
             >> s.RightTop;
}

#endif
