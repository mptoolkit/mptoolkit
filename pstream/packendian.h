// -*- C++ -*- $Id$

#if !defined(PACKENDIAN_H_DSHFIUYT437Y7LIUHYIURH8743Y7H)
#define PACKENDIAN_H_DSHFIUYT437Y7LIUHYIURH8743Y7H

#include "endian.h"
#include "currentformat.h"

template <EndianType Endian, typename Buffer, typename T>
inline
void pack_endian(Buffer& B, T const& x);

template <EndianType Endian, typename Buffer, typename T>
inline
void unpack_endian(Buffer& B, T& x);

#include "packendian.cc"

#endif
