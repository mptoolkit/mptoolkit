// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/prog_options.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
//
// Include wrapper for boost program_options library, adds FormatDefault() function
//

#if !defined(PROG_OPTIONS_H_SDCH348957Y23489J89A)
#define PROG_OPTIONS_H_SDCH348957Y23489J89A

#include "common/prog_opt_accum.h"
#include <string>
#include <boost/lexical_cast.hpp>
#include <sstream>

template <typename T>
std::string FormatDefault(std::string const& Text, T const& Value)
{
   std::ostringstream s;
   s.precision(6);
   s << Value << std::flush;
   return Text + " [default " + s.str() + ']';
}

#endif
