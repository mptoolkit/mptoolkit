// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/hash.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  hash.h

  Portable hashing.

  Created 2002-11-17 Ian McCulloch
*/

#if !defined(HASH_H_HDSFJKHUIRY879YHF8HQ38HF8OH389UWEO)
#define HASH_H_HDSFJKHUIRY879YHF8HQ38HF8OH389UWEO

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
#include "common/inttype.h"

namespace ext
{

// Returns a hash constructed such that the return value is unique with a probability
// of 1 in 2^32.  The data hashed includes the current time from gettimeofday(),
// the process id, the hostname, and the value from the previous call.
uint32_t get_unique();

} // namespace ext

#endif
