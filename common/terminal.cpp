// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/terminal.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "terminal.h"
#include <stdlib.h>
#include <unistd.h>
#include <string>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_TERMIOS_H
#include <termios.h>
#endif

#if HAVE_SYS_IOCTL_H
#include <sys/ioctl.h>
#endif

namespace terminal
{

int rows()
{
   // primary option: environment string LINES
   char const* rows_string = getenv("LINES");
   if (rows_string)
   {
      int n_rows = atoi(rows_string);
      if (n_rows > 0)
         return n_rows;
   }
   // fall back to ioctl
#ifdef TIOCGWINSZ
   winsize ws;
   if (!ioctl(1, TIOCGWINSZ, &ws) && ws.ws_row)
      return ws.ws_row;
#endif
   // last resort
   return 25;
}

int columns()
{
   // primary option: environment string COLUMNS
   char const* col_string = getenv("COLUMNS");
   if (col_string)
   {
      int n_cols = atoi(col_string);
      if (n_cols > 0)
         return n_cols;
   }
   // fall back to ioctl
#ifdef TIOCGWINSZ
   winsize ws;
   if (!ioctl(1, TIOCGWINSZ, &ws) && ws.ws_col)
      return ws.ws_col;
#endif
   // last resort
   return 80;
}

std::pair<int, int> tsize()
{
   return std::make_pair(rows(), columns());
}

bool is_cout_terminal()
{
   return isatty(1);
}

std::string color_code(color c)
{
   return "\e[" + std::to_string(int(c)) + 'm';
}

std::string color_text(std::string s, color c)
{
   return color_code(c) + s + color_code(color::Reset);
}

std::string color_text(std::string s, color c1, color c2)
{
   return color_code(c1) + color_code(c2) + s + color_code(color::Reset);
}

} // namespace terminal
