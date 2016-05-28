// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/sliceproxy.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  sliceproxy.h

  A proxy class for taking the slice of a vector

  Created 2005-01-12 Ian McCulloch
*/

#if !defined(SLICEPROXY_H_JCHDJKHGUIHRIHEURIEO)
#define SLICEPROXY_H_JCHDJKHGUIHRIHEURIEO

#include "slice.h"
#include "vectoroperations.h"

namespace LinearAlgebra
{

template <typename BaseReference, typename SliceType = Slice>
class SliceProxy
{
   public:

   private:
      SliceType Slice_;
      BaseReference Base_;
};


} // namespace LinearAlgebra

#endif
