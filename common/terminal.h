// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/terminal.h
//
// Copyright (C) 2006-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
   terminal.h

   Some simple terminal control functions.

   Created 2006-05-19 Ian McCulloch
*/

#if !defined(TERMINAL_H_DSJFHJFH43Y439Y437YYH9YREW)
#define TERMINAL_H_DSJFHJFH43Y439Y437YYH9YREW

#include <utility>

namespace terminal
{

// return the size of the output terminal, as (rows, cols)
std::pair<int, int> size();

// return the number of rows of the output terminal
int rows();

// return the number of columns in the output terminal
int columns();

} // namespace terminal


#endif
