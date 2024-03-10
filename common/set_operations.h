// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/set_operations.h
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

// Some set operations to construct set unions and intersection iterators on the fly.

// currently, only intersections are defined.

#if !defined(SET_OPERATIONS_H_HSDFJKFH894578934789)
#define SET_OPERATIONS_H_HSDFJKFH894578934789

#include <set>

template <typename T>
struct intersection_iterator
{
   typedef std::forward_iterator_tag              iterator_category;
   typedef typename T::const_iterator             iterator;
   typedef typename iterator::value_type          value_type;
   typedef typename iterator::difference_type     difference_type;
   typedef typename iterator::pointer             pointer;
   typedef typename iterator::reference           reference;

   intersection_iterator() {}

   intersection_iterator(iterator x_, iterator y_, iterator xEnd_, iterator yEnd_)
      : x(x_), y(y_), xEnd(xEnd_), yEnd(yEnd_) { this->regularize(); }

   bool operator==(intersection_iterator<T> const& other) const
   {
      // It is probably sufficient to have just x == other.x here
      return x == other.x;
   }

   bool operator!=(intersection_iterator<T> const& other) const
   {
      // It is probably sufficient to have just x != other.x here
      return x != other.x;
   }

   intersection_iterator& operator++()
   {
      ++x; ++y;
      this->regularize();
      return *this;
   }

   intersection_iterator operator++(int)
   {
      intersection_iterator temp = *this;
      ++*this;
      return temp;
   }

   reference operator*() const
   {
      return *x;
   }

   pointer operator->() const
   {
      return x.operator->();
   }

   private:

   void regularize()
   {
      if (x == xEnd)
      {
         y = yEnd;
         return;
      }
      else if (y == yEnd)
      {
         x = xEnd;
         return;
      }
      // else
      // neither x nor y is at the end
      while (*x < *y)
      {
         if (++x == xEnd)
         {
            y = yEnd;
            return;
         }
      }

      while (*y < *x)
      {
         if (++y == yEnd)
         {
            x = xEnd;
            return;
         }
      }
   }

   iterator x,y;
   iterator xEnd, yEnd;
};

template <typename T, typename Compare>
inline
intersection_iterator<std::set<T, Compare> >
set_intersection_begin(std::set<T, Compare> const& x,
                       std::set<T, Compare> const& y)
{
   return intersection_iterator<std::set<T, Compare> >(x.begin(), y.begin(), x.end(), y.end());
}

template <typename T, typename Compare>
inline
intersection_iterator<std::set<T, Compare> >
set_intersection_end(std::set<T, Compare> const& x,
                     std::set<T, Compare> const& y)
{
   return intersection_iterator<std::set<T, Compare> >(x.end(), y.end(), x.end(), y.end());
}

#endif
