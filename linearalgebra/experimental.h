// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/experimental.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

//
// FlattenedIterator
//

template <typename OuterIterType>
class FattenedIterator
{
   public:
      typedef typename OuterIterType                               outer_type;
      typedef typename outer_type::iterator                        iterator;
      typedef typename outer_type::const_iterator                  const_iterator;
      typedef typename IteratorTraits<iterator>::value_type        value_type;
      typedef typename IteratorTraits<iterator>::result_type       result_type;
      typedef typename IteratorTraits<iterator>::const_result_type const_result_type;

      result_type operator*() { return *Inner; }

      result_type operator[](difference_type i) { FlattenedIterator Temp(*this); Temp += i; return *Temp; }

      FlattenedIterator& operator+=(difference_type i) { 

   private:
     outer_type Outer;
     iterator Inner;
};
