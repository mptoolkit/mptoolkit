// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/valgrind-wrapper.h
//
// Copyright (C) 2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// valgrind wrapper - if we have valgrind.h then use that, otherwise
// declare some dummy macros

#if !defined(MPTOOLKIT_COMMON_VALGRIND_WRAPPER_H)
#define MPTOOLKIT_COMMON_VALGRIND_WRAPPER_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if defined(HAVE_VALGRIND_MEMCHECK_H)

#include <valgrind/memcheck.h>

#else

#define VALGRIND_MAKE_MEM_NOACCESS(a,b)
#define VALGRIND_MAKE_MEM_UNDEFINED(a,b)
#define VALGRIND_MAKE_MEM_DEFINED(a,b)
#define VALGRIND_MAKE_MEM_DEFINED_IF_ADDRESSABLE(a,b)
#define VALGRIND_CREATE_BLOCK(a,b,c)
#define VALGRIND_DISCARD(a)
#define VALGRIND_CHECK_MEM_IS_ADDRESSABLE(a,b)
#define VALGRIND_CHECK_MEM_IS_DEFINED(a,b)
#define VALGRIND_CHECK_VALUE_IS_DEFINED(a)

#endif // defined(HAVE_VALGRIND_MEMCHECK_H)

#endif
