// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// parser/number-parser.h
//
// Copyright (C) 2016-2017 Ian McCulloch <ian@qusim.net>
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

// Simple parser for complex numbers

#if !defined(MPTOOLKIT_NUMBER_PARSER_H)
#define MPTOOLKIT_NUMBER_PARSER_H

#include "mps/state_component.h"
#include "lattice/function.h"
#include <string>
#include <map>

std::complex<double>
ParseNumber(std::string const& Str, Function::ArgumentList const& Args);

inline
std::complex<double>
ParseNumber(std::string const& Str)
{
   return ParseNumber(Str, Function::ArgumentList());
}

#endif
