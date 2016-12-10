// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/copyright.h
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

#if !defined(MPTOOLKIT_MP_COPYRIGHT_H)
#define MPTOOLKIT_MP_COPYRIGHT_H

#include <iostream>
#include <string>

// Display copyright information for the GPL-3
void
print_copyright(std::ostream& out);

// Display warranty information, as recommended for GPL-3
void
print_warranty(std::ostream& out);

// Display copying information, as recommended for GPL-3
void
print_copying(std::ostream& out);

// Display citation list
void
print_citations(std::ostream& out);

// Given the program name, return a wikified version of a URL to the documentation.
std::string
Wikify(std::string const& x, bool itool = true);

// 
void print_copyright(std::ostream& out, std::string const& Category, std::string const& Name);

// The basename() function is useful in help messages for printing the program name
std::string basename(std::string const& FName);

// Escape an argument for bash (eg, to print the command-line arguments)
std::string
EscapeArgument(std::string const& s);

std::string
EscapeCommandline(int argc, char** argv);

// Print a useful preamble, consisting of the
// program name and arguments, and the date.
void print_preamble(std::ostream& out, int argc, char** argv);

#endif
