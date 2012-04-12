/* -*- C++ -*- $Id$

  Pack and unpack of intrinsic types to/from a buffer according to the specified format.
  
  Created 2004-01-17 Ian McCulloch
*/

#if !defined(PACKFORMAT_H_DHFHFJ43YR43YY4YR4Y)
#define PACKFORMAT_H_DHFHFJ43YR43YY4YR4Y

#include "binaryformats.h"

template <typename Format, typename Buffer, typename T>
inline
void pack_format(Buffer& B, T x);

template <typename Format, typename Buffer, typename T>
inline
void unpack_format(Buffer& B, T& x);

#include "packformat.cc"

#endif
