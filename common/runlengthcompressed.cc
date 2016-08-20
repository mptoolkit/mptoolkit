// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/runlengthcompressed.cc
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

#include <stack>
#include <list>
#include "common/trace.h"
#include <tuple>

// run_length_compressed<T>::const_iterator

template <typename T>
class run_length_compressed<T>::const_iterator
{
   public:
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef T value_type;
      typedef int difference_type;
      typedef T const* pointer;
      typedef T const& reference;

      const_iterator() {}

      const_iterator& operator++();
      const_iterator operator++(int);

      const_iterator& operator--();
      const_iterator operator--(int);

      bool operator==(const_iterator const& x)
      {
         // We need to check that the nodes are the same before checking the iterators,
         // otherwise we would be comparing iterators from different containers
         return Item == x.Item && (Nodes.empty() || (Nodes == x.Nodes && Current == x.Current));
      }

      bool operator!=(const_iterator const& x)
      {
         return !this->operator==(x);
      }

      pointer operator->() const;
      reference operator*() const;

   private:
      typedef typename run_length_repeat<T>::const_iterator repeat_iterator;
      typedef typename run_length_array<T>::const_iterator array_iterator;
      typedef boost::variant<repeat_iterator, array_iterator> iter_base;
      typedef std::stack<iter_base> IterStackType;
      typedef std::stack<run_length_compressed<T> const*> NodeStackType;

      NodeStackType Nodes;
      IterStackType Current;
      T const* Item;
      run_length_compressed<T> const* NodeTop;

      struct BeginTag {};
      struct EndTag {};

      const_iterator(run_length_compressed<T> const& c, BeginTag)
         : NodeTop(&c)
      {
         if (c.empty())
         {
            Item = NULL;
            return;
         }
         Nodes.push(NodeTop);
         boost::apply_visitor(IterBegin(&Nodes, &Current, &Item), c.data());
      }

      const_iterator(run_length_compressed<T> const& c, EndTag)
         : Item(NULL), NodeTop(&c)
      {
      }

      struct GetIterBegin;
      struct GetIterEnd;

      struct IterBegin;
      struct IterEndDec;

   friend class run_length_compressed<T>;
};

template <typename T>
struct run_length_compressed<T>::const_iterator::GetIterEnd
   : boost::static_visitor<typename run_length_compressed<T>::const_iterator::iter_base>
{
   iter_base operator()(T const&) const
   {
      PANIC("invalid iterator");
      return iter_base();
   }
   template <typename U>
   iter_base operator()(U const& x) const
   {
      return x.end();
   }
};

template <typename T>
struct run_length_compressed<T>::const_iterator::GetIterBegin
   : boost::static_visitor<typename run_length_compressed<T>::const_iterator::iter_base>
{
   iter_base operator()(T const&) const
   {
      PANIC("invalid iterator");
      return iter_base();
   }
   template <typename U>
   iter_base operator()(U const& x) const
   {
      return x.begin();
   }
};

template <typename T>
struct IterArrow : public boost::static_visitor<T*>
{
   template <typename U>
   T* operator()(U const& x) const
   {
      return x.operator->();
   }
};

template <typename T>
struct run_length_compressed<T>::const_iterator::IterBegin : public boost::static_visitor<>
{
   IterBegin(NodeStackType* Nodes_, IterStackType* Current_, T const** Item_)
      : Nodes(Nodes_), Current(Current_), Item(Item_) {}

   void operator()(T const& x) const
   {
      Nodes->pop();
      *Item = &x;
   }
   void operator()(run_length_repeat<T> const& x) const
   {
      Current->push(x.begin());
      Nodes->push(boost::apply_visitor(IterArrow<run_length_compressed<T> const>(), Current->top()));
      boost::apply_visitor(*this, x.nested());
   }
   void operator()(run_length_array<T> const& x) const
   {
      Current->push(x.begin());
      Nodes->push(boost::apply_visitor(IterArrow<run_length_compressed<T> const>(), Current->top()));
      boost::apply_visitor(*this, x.front());
   }

   // these overloads also allow this to work on the Current.top() iterator
   void operator()(repeat_iterator const& x) const
   {
      boost::apply_visitor(*this, *x);
   }
   void operator()(array_iterator const& x) const
   {
      boost::apply_visitor(*this, *x);
   }

   NodeStackType* Nodes;
   IterStackType* Current;
   T const** Item;
};

template <typename T>
struct run_length_compressed<T>::const_iterator::IterEndDec
   : public boost::static_visitor<>
{
   IterEndDec(NodeStackType* Nodes_, IterStackType* Current_, T const** Item_)
      : Nodes(Nodes_), Current(Current_), Item(Item_) {}

   void operator()(T const& x) const
   {
      Nodes->pop();
      *Item = &x;
   }
   void operator()(run_length_repeat<T> const& x) const
   {
      repeat_iterator i = x.end();
      --i;
      Current->push(i);
      Nodes->push(boost::apply_visitor(IterArrow<run_length_compressed<T> const>(), Current->top()));
      boost::apply_visitor(*this, *i);
   }
   void operator()(run_length_array<T> const& x) const
   {
      array_iterator i = x.end();
      --i;
      Current->push(i);
      Nodes->push(boost::apply_visitor(IterArrow<run_length_compressed<T> const>(), Current->top()));
      boost::apply_visitor(*this, *i);
   }

   // these overloads also allow this to work on the Current.top() iterator
   void operator()(repeat_iterator& x) const
   {
      --x;
      Nodes->push(&*x);
      boost::apply_visitor(*this, *x);
   }
   void operator()(array_iterator& x) const
   {
      --x;
      Nodes->push(&*x);
      boost::apply_visitor(*this, *x);
   }

   NodeStackType* Nodes;
   IterStackType* Current;
   T const** Item;
};

struct IterIncrement : public boost::static_visitor<>
{
   template <typename U>
   void operator()(U& x) const
   {
      ++x;
   }
};

struct IterDecrement : public boost::static_visitor<>
{
   template <typename U>
   void operator()(U& x) const
   {
      --x;
   }
};

template <typename T>
typename run_length_compressed<T>::const_iterator&
run_length_compressed<T>::const_iterator::operator++()
{
   if (!Current.empty())
   {
      boost::apply_visitor(IterIncrement(), Current.top());
      while (!Current.empty() && Current.top() == boost::apply_visitor(GetIterEnd(), Nodes.top()->data()))
      {
         Current.pop();
         Nodes.pop();
         if (!Current.empty())
            boost::apply_visitor(IterIncrement(), Current.top());
      }
   }
   if (Current.empty()) // iterator over a single item is a trivial iterator
      Item = NULL;
   else
   {
      // drill down until we get a leaf node
      Nodes.push(boost::apply_visitor(IterArrow<run_length_compressed<T> const>(), Current.top()));
      boost::apply_visitor(IterBegin(&Nodes, &Current, &Item), Current.top());
   }
   return *this;
}

template <typename T>
typename run_length_compressed<T>::const_iterator&
run_length_compressed<T>::const_iterator::operator--()
{
   if (!Item)
   {
      DEBUG_CHECK(Current.empty());
      Nodes.push(NodeTop);
      boost::apply_visitor(IterEndDec(&Nodes, &Current, &Item), *NodeTop);
      return *this;
   }
   if (Current.empty())
   {
      PANIC("Invalid iterator");
      return *this;
   }
   // else while we are at the beginning of a node, move back up the heirachy to an iterator we can safely decrement
   while (!Current.empty() && Current.top() == boost::apply_visitor(GetIterBegin(), Nodes.top()->data()))
   {
      Current.pop();
      Nodes.pop();
   }
   if (Current.empty())
   {
      PANIC("Attempt to decrement iterator off the begin of the container");
   }
   // Now Current.top() is an iterator we can decrement, then drill down until we get a leaf node
   boost::apply_visitor(IterEndDec(&Nodes, &Current, &Item), Current.top());
   return *this;
}

template <typename T>
typename run_length_compressed<T>::const_iterator
run_length_compressed<T>::const_iterator::operator++(int)
{
   const_iterator Temp(*this);
   this->operator++();
   return Temp;
}

template <typename T>
T const*
run_length_compressed<T>::const_iterator::operator->() const
{
   return Item;
}

template <typename T>
T const&
run_length_compressed<T>::const_iterator::operator*() const
{
   return *Item;
}

// run_length_array

template <typename T>
inline
run_length_array<T>::run_length_array(T const& x)
   : Data(run_length_compressed<T>(x))
{
}

template <typename T>
run_length_array<T>::run_length_array(run_length_compressed<T> const& x)
{
   if (run_length_array<T> const* a = boost::get<run_length_array<T> >(&x.data()))
      Data = a->Data;
   else
   {
      Data.clear();
      Data.push_back(x);
      //Data = data_type(x);
   }
}

template <typename T>
run_length_array<T>::run_length_array(run_length_repeat<T> const& x)
{
   if (x.size() > 1)
      Data = data_type(1, run_length_compressed<T>(x));
   else if (x.size() == 1)
      Data = data_type(1, x.nested());
}

template <typename T>
template <typename FwdIter>
run_length_array<T>::run_length_array(FwdIter first, FwdIter last)
   : Data(first, last)
{
   // TODO:   if (Data.size() == 1)
}

template <typename T>
int
run_length_array<T>::logical_size() const
{
   int Result = 0;
   for (const_iterator I = this->begin(); I != this->end(); ++I)
   {
      Result += I->size();
   }
   return Result;
}

template <typename T>
inline
run_length_compressed<T> const& run_length_array<T>::front() const
{
   DEBUG_CHECK(!Data.empty());
   return Data.front();
}

template <typename T>
inline
run_length_compressed<T> const& run_length_array<T>::back() const
{
   DEBUG_CHECK(!Data.empty());
   return Data.back();
}

template <typename T>
inline
void run_length_array<T>::push_front(T const& x)
{
   Data.push_front(run_length_compressed<T>(x));
}

template <typename T>
void run_length_array<T>::push_front(run_length_repeat<T> const& x)
{
   if (x.size() > 0)
   {
      if (x.size() == 1)
         Data.push_front(x.nested());
      else
          Data.push_front(run_length_compressed<T>(x));
   }
}

template <typename T>
inline
void run_length_array<T>::push_front(run_length_array<T> const& x)
{
   Data.insert(Data.begin(), x.begin(), x.end());
}

template <typename T>
void run_length_array<T>::push_front(run_length_compressed<T> const& x)
{
   if (run_length_array<T> const* a = boost::get<run_length_array<T> >(&x.data()))
      Data.insert(Data.begin(), a->begin(), a->end());
   else
      Data.push_front(x);
}

template <typename T>
inline
void run_length_array<T>::push_back(T const& x)
{
   Data.push_back(run_length_compressed<T>(x));
}

template <typename T>
void run_length_array<T>::push_back(run_length_repeat<T> const& x)
{
   if (x.size() > 0)
   {
      if (x.size() == 1)
         Data.push_back(x.nested());
      else
          Data.push_back(run_length_compressed<T>(x));
   }
}

template <typename T>
inline
void run_length_array<T>::push_back(run_length_array<T> const& x)
{
   Data.insert(Data.end(), x.begin(), x.end());
}

template <typename T>
void run_length_array<T>::push_back(run_length_compressed<T> const& x)
{
   if (run_length_array<T> const* a = boost::get<run_length_array<T> >(&x.data()))
      Data.insert(Data.end(), a->begin(), a->end());
   else
      Data.push_back(x);
}

template <typename T>
template <typename FwdIter>
inline
void
run_length_array<T>::append(FwdIter start, FwdIter finish)
{
   Data.insert(Data.end(), start, finish);
}

// run_length_repeat

template <typename T>
inline
run_length_repeat<T>::run_length_repeat(int RepeatCount, run_length_compressed<T> const& x)
   : Value(x), Count(RepeatCount)
{
   // if x is itself a run_length_repeat, then collapse it
   if (run_length_repeat<T> const* a = boost::get<run_length_repeat<T> >(&x.data()))
   {
      Value = a->nested();
      Count = RepeatCount * a->size();
   }
}

template <typename T>
inline
run_length_compressed<T> const& run_length_repeat<T>::front() const
{
   DEBUG_CHECK(Count != 0);
   return Value;
}

template <typename T>
inline
run_length_compressed<T> const& run_length_repeat<T>::back() const
{
   DEBUG_CHECK(Count != 0);
   return Value;
}

// run_length_compressed

template <typename T>
inline
run_length_compressed<T>::run_length_compressed()
   : Data(run_length_array<T>())
{
}

template <typename T>
inline
run_length_compressed<T>::run_length_compressed(T const& x)
   : Data(x)
{
}

template <typename T>
inline
run_length_compressed<T>::run_length_compressed(int Size, T const& x)
   : Data(run_length_repeat<T>(Size, x))
{
   // corner case: repeat count is 1
   if (Size == 1)
      Data = x;
}

template <typename T>
inline
run_length_compressed<T>::run_length_compressed(run_length_repeat<T> const& x)
   : Data(x)
{
   // corner case: repeat count is 1
   if (x.size() == 1)
      Data = x.nested().data();
}

template <typename T>
inline
run_length_compressed<T>::run_length_compressed(run_length_array<T> const& x)
   : Data(x)
{
   // corner case: array size is 1
   if (x.size() == 1)
      Data = x.front().data();
}

template <typename T>
inline
run_length_compressed<T>::run_length_compressed(int Size, run_length_compressed<T> const& x)
   : Data(run_length_repeat<T>(Size, x))
{
   // corner case: repeat count is 1
   if (Size == 1)
      Data = x.data();
}

template <typename T>
template <typename FwdIter>
run_length_compressed<T>::run_length_compressed(FwdIter start, FwdIter finish)
   : Data(run_length_array<T>(start, finish))
{
   // corner case: array length is 1
   if (std::distance(start, finish) == 1)
      Data = *start;
}

template <typename T>
run_length_compressed<T>::run_length_compressed(T const& x1, T const& x2)
{
   std::list<T> ls(1, x1);
   ls.push_back(x2);
   Data = run_length_array<T>(ls.begin(), ls.end());
}

template <typename T>
run_length_compressed<T>::run_length_compressed(T const& x1, T const& x2, T const& x3)
{
   std::list<T> ls(1, x1);
   ls.push_back(x2);
   ls.push_back(x3);
   Data = run_length_array<T>(ls.begin(), ls.end());
}

template <typename T>
run_length_compressed<T>::run_length_compressed(T const& x1, T const& x2,
                                                T const& x3, T const& x4)
{
   std::list<T> ls(1, x1);
   ls.push_back(x2);
   ls.push_back(x3);
   ls.push_back(x4);
   Data = run_length_array<T>(ls.begin(), ls.end());
}

template <typename T>
typename run_length_compressed<T>::const_iterator
run_length_compressed<T>::begin() const
{
   return const_iterator(*this, typename const_iterator::BeginTag());
}

template <typename T>
typename run_length_compressed<T>::const_iterator
run_length_compressed<T>::end() const
{
   return const_iterator(*this, typename const_iterator::EndTag());
}

template <typename T>
struct RLE_size : public boost::static_visitor<int>
{
   int operator()(T const&) const
   {
      return 1;
   }

   int operator()(run_length_repeat<T> const& x) const
   {
      return x.logical_size();
   }

   int operator()(run_length_array<T> const& x) const
   {
      return x.logical_size();
   }
};

template <typename T>
int
run_length_compressed<T>::size() const
{
   return boost::apply_visitor(RLE_size<T>(), Data);
}

template <typename T>
struct RLE_empty : public boost::static_visitor<bool>
{
   bool operator()(T const&) const
   {
      return false;
   }
   template <typename U>
   bool operator()(U const& x) const
   {
      return x.empty();
   }
};

template <typename T>
bool
run_length_compressed<T>::empty() const
{
   return boost::apply_visitor(RLE_empty<T>(), Data);
}

template <typename T>
struct RLE_front : public boost::static_visitor<T const&>
{
   T const& operator()(T const& x) const
   {
      return x;
   }
   template <typename U>
   T const& operator()(U const& x) const
   {
      return  boost::apply_visitor(RLE_front<T>(), x.front());
   }
};

template <typename T>
struct RLE_back : public boost::static_visitor<T const&>
{
   T const& operator()(T const& x) const
   {
      return x;
   }
   template <typename U>
   T const& operator()(U const& x) const
   {
      return boost::apply_visitor(RLE_back<T>(), x.back());
   }
};

template <typename T>
T const& run_length_compressed<T>::front() const
{
   return boost::apply_visitor(RLE_front<T>(), Data);
}

template <typename T>
T const& run_length_compressed<T>::back() const
{
   return boost::apply_visitor(RLE_back<T>(), Data);
}

template <typename T>
void run_length_compressed<T>::push_back(T const& x)
{
   if (this->empty())
      Data = x;
   else
   {
      run_length_array<T>* a = boost::get<run_length_array<T> >(&Data);
      if (!a)
      {
         Data = run_length_array<T>(*this);
         a = boost::get<run_length_array<T> >(&Data);
      }
      a->push_back(x);
   }
}

template <typename T>
void run_length_compressed<T>::push_back(run_length_array<T> const& x)
{
   if (x.empty())
      return;
   if (this->empty())
   {
      Data = x;
      return;
   }
   if (x.size() == 1)
   {
      this->push_back(x.front());
      return;
   }
   // else
   run_length_array<T>* a = boost::get<run_length_array<T> >(&Data);
   if (!a)
   {
      Data = run_length_array<T>(*this);
      a = boost::get<run_length_array<T> >(&Data);
   }
   a->push_back(x);
}

template <typename T>
void run_length_compressed<T>::push_back(run_length_repeat<T> const& x)
{
   if (x.empty())
      return;
   if (this->empty())
   {
      Data = x;
      return;
   }
   if (x.size() == 1)
   {
      this->push_back(x.front());
      return;
   }
   // else
   run_length_array<T>* a = boost::get<run_length_array<T> >(&Data);
   if (!a)
   {
      Data = run_length_array<T>(*this);
      a = boost::get<run_length_array<T> >(&Data);
   }
   a->push_back(x);
}

template <typename T>
struct RLE_push_back : boost::static_visitor<>
{
   RLE_push_back(run_length_compressed<T>* Data_) : Data(Data_) {}
   template <typename U>
   void operator()(U const& x) const
   {
      Data->push_back(x);
   }
   run_length_compressed<T>* Data;
};

template <typename T>
inline
void run_length_compressed<T>::push_back(run_length_compressed<T> const& x)
{
   boost::apply_visitor(RLE_push_back<T>(this), x.data());
}

template <typename T>
void run_length_compressed<T>::push_front(T const& x)
{
   if (this->empty())
      Data = x;
   else
   {
      run_length_array<T>* a = boost::get<run_length_array<T> >(&Data);
      if (!a)
      {
         Data = run_length_array<T>(*this);
         a = boost::get<run_length_array<T> >(&Data);
      }
      a->push_front(x);
   }
}

template <typename T>
void run_length_compressed<T>::push_front(run_length_array<T> const& x)
{
   if (x.empty())
      return;
   if (this->empty())
   {
      Data = x;
      return;
   }
   if (x.size() == 1)
   {
      this->push_front(x.front());
      return;
   }
   // else
   run_length_array<T>* a = boost::get<run_length_array<T> >(&Data);
   if (!a)
   {
      Data = run_length_array<T>(*this);
      a = boost::get<run_length_array<T> >(&Data);
   }
   a->push_front(x);
}

template <typename T>
void run_length_compressed<T>::push_front(run_length_repeat<T> const& x)
{
   if (x.empty())
      return;
   if (this->empty())
   {
      Data = x;
      return;
   }
   if (x.size() == 1)
   {
      this->push_front(x.front());
      return;
   }
   // else
   run_length_array<T>* a = boost::get<run_length_array<T> >(&Data);
   if (!a)
   {
      Data = run_length_array<T>(*this);
      a = boost::get<run_length_array<T> >(&Data);
   }
   a->push_front(x);
}

template <typename T>
struct RLE_push_front : boost::static_visitor<>
{
   RLE_push_front(run_length_compressed<T>* Data_) : Data(Data_) {}
   template <typename U>
   void operator()(U const& x) const
   {
      Data->push_front(x);
   }
   run_length_compressed<T>* Data;
};

template <typename T>
inline
void run_length_compressed<T>::push_front(run_length_compressed<T> const& x)
{
   boost::apply_visitor(RLE_push_front<T>(this), x.data());
}

template <typename T>
template <typename FwdIter>
void run_length_compressed<T>::append(FwdIter first, FwdIter last)
{
   run_length_array<T>* a = boost::get<run_length_array<T> >(&Data);
   if (!a)
   {
      Data = run_length_array<T>(*this);
      a = boost::get<run_length_array<T> >(&Data);
   }
   if (a->empty())
      Data = run_length_array<T>(first, last);
   else
      a->append(first, last);
}

template <typename T>
void canonicalize(run_length_compressed<T>& x)
{
   if (run_length_array<T>* a = boost::get<run_length_array<T> >(&x.data()))
   {
      typename run_length_array<T>::iterator I = a->begin();
      while (I != a->end())
      {
         canonicalize(*I);
         // if it is an empty array/repeat, remove it
         if (I->size() == 0)
            I = a->erase(I);
         else if (run_length_array<T>* b = boost::get<run_length_array<T> >(*I))
         {
            // otherwise, if it is a nested array, flatten it
            int Loc = I - a->begin();
            std::deque<run_length_compressed<T> > Temp;
            b->swap(Temp);
            I = a->erase(I);  // erase the element pointed to by b
            a->insert(I, Temp.begin(), Temp.end());
            I = a->begin() + Loc; // and our new I is at the start of the inserted elements
         }
         else
            ++I;
      }
   }
   else if (run_length_repeat<T>* a = boost::get<run_length_repeat<T> >(&x.data()))
   {
      // first, check for zero
      if (a->size() == 0 || a->nested_size() == 0)
      {
         x = run_length_array<T>();
         return;
      }
      // else
      canonicalize(a->nested());
      if (a->size() == 1)
         x = a->nested();
   }
   // else it is a T itself, no need to do anything
}

template <typename T>
run_length_compressed<T>
repeat(run_length_compressed<T> const& x, int RepeatCount)
{
   if (RepeatCount == 0)
      return run_length_array<T>();
   // else
   if (RepeatCount == 1)
      return x;
   // else
   return run_length_repeat<T>(RepeatCount, x);
}

template <typename T>
run_length_compressed<T>
join(run_length_compressed<T> const& x, run_length_compressed<T> const& y)
{
   if (x.empty())
      return y;
   // else
   if (y.empty())
      return x;
   run_length_array<T> r(x);
   r.push_back(y);
   return r;
}

template <typename T>
struct RLE_leaf_count : public boost::static_visitor<int>
{
   int operator()(T const&) const
   {
      return 1;
   }
   int operator()(run_length_repeat<T> const& x) const
   {
      return x.nested().leaf_count();
   }
   int operator()(run_length_array<T> const& x) const
   {
      int Result = 0;
      for (typename run_length_array<T>::const_iterator I = x.begin(); I != x.end(); ++I)
         Result += I->leaf_count();
      return Result;
   }
};

template <typename T>
int
run_length_compressed<T>::leaf_count() const
{
   return boost::apply_visitor(RLE_leaf_count<T>(), Data);
}

template <typename T>
struct RLE_node_count : public boost::static_visitor<int>
{
   int operator()(T const&) const
   {
      return 0;
   }
   int operator()(run_length_repeat<T> const& x) const
   {
      return x.nested().node_count() + 1;
   }
   int operator()(run_length_array<T> const& x) const
   {
      int Result = 1;  // 1 for the current node
      for (typename run_length_array<T>::const_iterator I = x.begin(); I != x.end(); ++I)
         Result += I->node_count();
      return Result;
   }
};

template <typename T>
int
run_length_compressed<T>::node_count() const
{
   return boost::apply_visitor(RLE_node_count<T>(), Data);
}

template <typename T, typename Visitor>
inline
typename Visitor::result_type
apply_visitor(Visitor const& v, run_length_compressed<T> const& x)
{
   return boost::apply_visitor(v, x.data());
}

template <typename T, typename Visitor>
inline
typename Visitor::result_type
apply_visitor(Visitor const& v, run_length_compressed<T>& x)
{
   return boost::apply_visitor(v, x.data());
}

template <typename T, typename Visitor>
inline
typename Visitor::result_type
apply_visitor(Visitor const& v,
              run_length_compressed<T> const& x,
              run_length_compressed<T> const& y)
{
   return boost::apply_visitor(v, x.data(), y.data());
}

//
// split
//

// Visitor functor to implement the split() function
template <typename T>
struct DoSplit
   : public boost::static_visitor<std::pair<run_length_compressed<T>, run_length_compressed<T> > >
{
   typedef run_length_compressed<T> value_type;
   typedef run_length_repeat<T> repeat_type;
   typedef run_length_array<T> array_type;
   typedef std::pair<value_type, value_type> result_type;

   DoSplit(int Loc_) : Loc(Loc_) {}

   result_type operator()(T const& x) const
   {
      PANIC("DoSplit: attempt to split a single component")(Loc);
      return result_type();
   }

   result_type operator()(repeat_type const& x) const
   {
      // the easy case is where Loc is a multiple of nested.size()
      int LocCount = Loc / x.nested().size();
      DEBUG_CHECK(LocCount < x.size())(Loc)(LocCount)(x.size());
      int LocRem = Loc % x.nested().size();
      if (LocRem == 0)
      {
         return result_type(value_type(LocCount, x.nested()),
                            value_type( x.size()-LocCount, x.nested()));
      }
      // else
      // take the closest integral number of periods and make an array from that
      // plus the remainder
      array_type ResultL(repeat_type(LocCount, x.nested()));
      std::pair<value_type, value_type> JoinPart = x.nested().apply_visitor(DoSplit(LocRem));
      ResultL.push_back(JoinPart.first);

      int RemCount = x.size() - LocCount - 1;  // -1 for the part we just split
      if (RemCount == 0) // if we split the right most part, we have finished
         return result_type(ResultL, JoinPart.second);
      else
      {
         // otherwise add the rest of the array
         array_type ResultR(JoinPart.second);
         ResultR.push_back(repeat_type(RemCount, x.nested()));
         return result_type(ResultL, ResultR);
      }
   }

   result_type operator()(array_type const& a) const
   {
      // Firstly, get an iterator into the array that points to the node containing Loc
      int CurrentLoc = 0;
      typename array_type::const_iterator I = a.begin();
      while (CurrentLoc < Loc)
      {
         DEBUG_CHECK(I != a.end())(CurrentLoc)(Loc)(a.size())(a.logical_size());
         CurrentLoc += I->size();
         ++I;
      }
      // The left-half
      // If the split is exactly on the boundary of array elements, we are in luck
      if (CurrentLoc == Loc)
         return result_type(value_type(a.begin(), I), value_type(I, a.end()));
      // else we have to split the node itself
      --I;
      value_type ResultL(a.begin(), I);
      value_type Temp, ResultR;
      std::tie(Temp, ResultR) = split(*I, Loc - CurrentLoc + I->size());
      ResultL.push_back(Temp);
      ResultR.append(++I, a.end());
      DEBUG_CHECK_EQUAL(ResultL.size(), Loc);
      return result_type(ResultL, ResultR);
   }

   int const Loc;
};

template <typename T>
std::pair<run_length_compressed<T>, run_length_compressed<T> >
split(run_length_compressed<T> const& x, int Loc)
{
   CHECK(Loc >= 0 && Loc <= x.size())(Loc)(x.size());
   if (Loc == 0)
      return std::pair<run_length_compressed<T>,
         run_length_compressed<T> >(run_length_compressed<T>(), x);
   if (Loc == x.size())
      return std::pair<run_length_compressed<T>,
         run_length_compressed<T> >(x, run_length_compressed<T>());

   return x.apply_visitor(DoSplit<T>(Loc));
}

template <typename T>
std::tuple<run_length_compressed<T>, T, run_length_compressed<T> >
split_lcr(run_length_compressed<T> const& x, int Loc)
{
   CHECK(Loc >= 0 && Loc < x.size())(Loc);
   std::pair<run_length_compressed<T>, run_length_compressed<T> > l_cr = split(x, Loc);
   std::pair<run_length_compressed<T>, run_length_compressed<T> > cr = split(l_cr.second, 1);
   return std::make_tuple(l_cr.first, cr.first.front(), cr.second);
}

#if defined(USE_PSTREAM)
template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, run_length_compressed<T> const& x)
{
   return out << x.data();
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, run_length_repeat<T> const& x)
{
   int sz = x.size();
   return out << sz << x.nested();
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, run_length_array<T> const& x)
{
   return out << x.data();
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, run_length_compressed<T>& x)
{
   return in >> x.data();
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, run_length_repeat<T>& x)
{
   int sz;
   run_length_compressed<T> Nested;
   in >> sz >> Nested;
   x = run_length_repeat<T>(sz, Nested);
   return in;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, run_length_array<T>& x)
{
   return in >> x.data();
}

#endif
