// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/proxycontainer.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ian@qusim.net>
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
/*
  proxycontainer.h

  A simple pseudo-container for an iterator range.

  Created 2002-11-13 Ian McCulloch

  A proxy_container represents a fixed size iterator range.  While it is surely useful in other contexts,
  it is primarily designed for streaming.  a proxy_container knows what size it is,
  so the size is NOT explicitly inserted/extracted.

  The copy ctor and assignment use lightweight semantics, ie we just assign begin/end iterators rather than
  copy all the elements.  This is the correct semantics if we view proxy_container as simply an iterator range.
  However if we view proxy_container as a reference to another container, then the semantics should be
  deep-copy; in that case the copy constructor should not be defined and assignment should only work if
  the assignee is the same length as the original, in which case we copy all the elements.
  This also affects swap().

  Which to use?

  Currently there is no default constructor.  If we decide on lightweight semantics, then we should define one.
  It doesn't make sense to have a default ctor if we use deep-copy semantics.
*/

#if !defined(PROXYCONTAINER_H_CSH3784YR78FH785TYA7E5YT89P7YWE)
#define PROXYCONTAINER_H_CSH3784YR78FH785TYA7E5YT89P7YWE

#include <algorithm>

//
// AddConstIter
//
// A proxy iterator that turns a non-const iterator into a const iterator.
//

template <class Iter>
class AddConstIter;

template <class Iter>
bool operator==(AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
template <class Iter>
bool operator<(AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
template <class Iter>
bool operator!=(AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
template <class Iter>
bool operator>(AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
template <class Iter>
bool operator>=(AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
template <class Iter>
inline bool operator<=(AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
template <class Iter>
typename AddConstIter<Iter>::difference_type operator-(AddConstIter<Iter> const& x, AddConstIter<Iter> const& y);
template <class Iter>
AddConstIter<Iter> operator+(typename AddConstIter<Iter>::difference_type n, AddConstIter<Iter> const& x);

template <class Iter>
class AddConstIter
{
   public:
      typedef std::iterator_traits<Iter>::value_type        value_type;
      typedef std::iterator_traits<Iter>::difference_type   difference_type;
      typedef std::iterator_traits<Iter>::iterator_category iterator_category;
      typedef value_type const*                             pointer;
      typedef value_type const&                             reference;

      AddConstIter() {}

      AddConstIter(Iter const& i_) : i(i_) {}

      AddConstIter& operator=(Iter const& i_) { i = i_; return *this; }

      reference operator*() const { return *i; }
      pointer operator->() const { return i; }

      AddConstIterator& operator++() const { ++i; return *this; }
      AddConstIterator operator++(int) { AddConstIterator Temp = *this; ++i; return Temp; }

      AddConstIterator& operator--() const { --i; return *this; }
      AddConstIterator operator--(int) { AddConstIterator Temp = *this; --i; return Temp; }

      AddConstIterator operator+(difference_type n) const { return AddConstIterator(i + n); }
      AddConstIterator operator-(difference_type n) const { return AddConstIterator(i - n); }

      AddConstIterator operator+=(difference_type n) { i += n; return *this; }
      AddConstIterator operator-=(difference_type n) { i -= n; return *this; }

      reference operator[](difference_type n) const { return i[n+i]; }

   private:
      Iter i;

   friend bool operator==(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
   friend bool operator!=(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
   friend bool operator<(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
   friend bool operator<=(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
   friend bool operator>(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
   friend bool operator>=(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y);
   friend difference_type operator-(AddConstIter<Iter> const& x, AddConstIter<Iter> const& y);
   friend AddConstIter<Iter> operator+(difference_type n, AddConstIter<Iter> const& x);
};

template <class Iter>
inline bool operator==(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y){ return y.i == x.i; }
template <class Iter>
inline bool operator<(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y){ return y.i < x.i; }
template <class Iter>
inline bool operator!=(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y){ return y.i != x.i; }
template <class Iter>
inline bool operator>(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y){ return y.i > x.i; }
template <class Iter>
inline bool operator>=(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y){ return y.i >= x.i; }
template <class Iter>
inline bool operator<=(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y){ return y.i <= x.i; }

template <class Iter>
inline
typename AddConstIter<Iter>::difference_type operator-(const AddConstIter<Iter>& x, const AddConstIter<Iter>& y)
{
   return y.i - x.i;
}

template <class Iter>
inline
AddConstIter<Iter> operator+(typename AddConstIter<Iter>::difference_type n, AddConstIter<Iter>& x)
{
   return AddConstIter<Iter>(x.i + n);
}

template <class Iter>
class proxy_container
{
   public:
      typedef Iter iterator;
      typedef std::iterator_traits<iterator>::value_type                             value_type;
      typedef std::iterator_traits<iterator>::pointer                                pointer;
      typedef std::iterator_traits<iterator>::reference                              reference;
      typedef std::iterator_traits<iterator>::difference_type                        difference_type;
      typedef size_t                                                                 size_type;

      // As an optimization, if the pointer type is (value_type const*), then the iterator is already const.
      // Otherwise, make a AddConstIter out of it.
      typedef typename select_type<compare_types<pointer, value_type const*>::equal,
                                   iterator,
                                   AddConstIter<iterator> >::type                    const_iterator;

      typedef std::iterator_traits<const_iterator>::pointer                          const_pointer;
      typedef std::iterator_traits<const_iterator>::reference                        const_reference;

      proxy_container(Iter start_, Iter finish_) : start(start_), finish(finish_) {}

      // use compiler defined copy ctor, copy assignment and dtor.

      iterator begin() { return start; }
      iterator end() { return finish; }

      const_iterator begin() const { return start; }
      const_iterator end() const { return finish; }

      size_type size() const { return finish - start; }
      size_type max_size() const { return size(); }

      bool empty() const { return start == finish; }

      reference operator[](size_t n) { return start[n]; };
      const_reference operator[](size_t n) const { return start[n]; };

      void swap(proxy_container<Iter>& other)
      {
         using std::swap;
         swap(start, other.start);
         swap(finish, other.finish);
      }

   private:
      iterator start, finish;
};

template <class Iter>
bool operator==(proxy_container<Iter> const& x, proxy_container<Iter> const& y)
{
    return x.size() == y.size() && equal(x.begin(), x.end(), y.begin());
}

template <class Iter>
inline bool operator<(proxy_container<Iter> const& x, proxy_container<Iter> const& y)
{
    return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
}

template <class Iter>
inline bool operator!=(proxy_container<Iter> const& x, proxy_container<Iter> const& y) { return !(x == y); }
template <class Iter>
inline bool operator>(proxy_container<Iter> const& x, proxy_container<Iter> const& y) { return y < x; }
template <class Iter>
inline bool operator>=(proxy_container<Iter> const& x, proxy_container<Iter> const& y) { return !(x < y); }
template <class Iter>
inline bool operator<=(proxy_container<Iter> const& x, proxy_container<Iter> const& y) { return !(y < x); }

#if 0

template <class Container>
inline
proxy_container<typename Container::iterator> make_proxy_container(Container& c)
{
   return proxy_container<typename Container::iterator>(c.begin(), c.end());
}

template <class Container>
inline
proxy_container<typename Container::const_iterator> make_proxy_container(Container const& c)
{
   return proxy_container<typename Container::const_iterator>(c.begin(), c.end());
}

#else

template <class Container>
struct HackIntelBug
{
   typedef Container::iterator IterType;
};

template <class Container>
struct HackIntelBug<Container const>
{
   typedef Container::const_iterator IterType;
};

template <class Container>
inline
proxy_container<typename HackIntelBug<Container>::IterType> make_proxy_container(Container& c)
{
   return proxy_container<typename HackIntelBug<Container>::IterType>(c.begin(), c.end());
}

#endif

#endif
