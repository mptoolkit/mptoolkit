// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// parser/color-output.h
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

#if !defined(MPTOOLKIT_PARSER_COLOR_OUTPUT_H)
#define MPTOOLKIT_PARSER_COLOR_OUTPUT_H

#include <string>
#include <boost/lexical_cast.hpp>

// http://misc.flogisoft.com/bash/tip_colors_and_formatting

// returns true if std::cout is an actual terminal, false otherwise (eg redirected to a file)
bool is_cout_terminal();

enum class TerminalColor
{
   Reset = 0,
   Bold = 1,
   Dim = 2,
   Underline = 4,
   Default = 39,
   Black = 30,
   Red = 31,
   Green = 32,
   Yellow = 33,
   Blue = 34,
   Magenta = 35,
   Cyan = 36,
   LightGray = 37,
   DarkGray = 90,
   LightRed = 91,
   LightGreen = 92,
   LightYellow = 93,
   LightBlue = 94,
   LightMagenta = 95,
   LightCyan = 96,
   White = 97 
};

inline
std::string ColorCode(TerminalColor c)
{
   return "\e[" + boost::lexical_cast<std::string>(int(c)) + 'm';
}

inline
std::string ColorText(std::string s, TerminalColor c)
{
   return ColorCode(c) + s + ColorCode(TerminalColor::Reset);
}

inline
std::string ColorText(std::string s, TerminalColor c1, TerminalColor c2)
{
   return ColorCode(c1) + ColorCode(c2) + s + ColorCode(TerminalColor::Reset);
}


#endif
