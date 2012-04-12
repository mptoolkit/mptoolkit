/* -*- C++ -*- $Id$

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
