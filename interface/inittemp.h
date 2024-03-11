// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// interface/inittemp.h
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
//
// Helper functions for initialization of the persistent heap
//

#if !defined(INITTEMP_H_ASDCKHWEIOUFYWIUOUY8957893476789Y89)
#define INITTEMP_H_ASDCKHWEIOUFYWIUOUY8957893476789Y89

#if defined(HAVE_CONFIG_H)
# include "config.h"
#else
# if !defined(DEFAULT_PAGE_SIZE)
#  define DEFAULT_PAGE_SIZE 65546
# endif
# if !defined(DEFAULT_PAGE_CACHE_SIZE)
#  define DEFAULT_PAGE_CACHE_SIZE 67108864
# endif
#endif

namespace mp_pheap
{

// returns the page size to use, either MP_PAGE_SIZE environment or DEFAULT_PAGE_SIZE
int PageSize();

// returns the page cache size to use, either MP_PAGE_CACHE_SIZE environment
// or DEFAULT_PAGE_CACHE_SIZE
long CacheSize();

void InitializeTempPHeap(bool Verbose = false);

} // namespace mp_pheap

#endif
