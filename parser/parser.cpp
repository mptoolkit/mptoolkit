// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// parser/parser.cpp
//
// Copyright (C) 2016-2022 Ian McCulloch <ian@qusim.net>
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

#include "parser.h"
#include "common/terminal.h"
#include <fstream>

namespace Parser
{

//
// Color support
//

bool Color = terminal::is_cout_terminal();

void SetShowColors(ShowColors c)
{
   if (c == ColorNever)
      Color = false;
   else if (c == ColorAuto)
      Color = terminal::is_cout_terminal();
   else if (c == ColorAlways)
      Color = true;
   else
      Color = false;
}

std::string ConditionalColorText(std::string s, terminal::color c)
{
   return Color ? terminal::color_text(s, c) : s;
}

std::string ConditionalColorText(std::string s, terminal::color c1, terminal::color c2)
{
   return Color ? terminal::color_text(s, c1, c2) : s;
}

std::string ColorHighlight(std::string s)
{
   return ConditionalColorText(s, terminal::color::Bold);
}

std::string ColorHighlight(char c)
{
   return ConditionalColorText(std::string(1,c), terminal::color::Bold);
}

std::string ColorError(std::string s)
{
   return ConditionalColorText(s, terminal::color::Red, terminal::color::Bold);
}

std::string ColorWarning(std::string s)
{
   return ConditionalColorText(s, terminal::color::Magenta, terminal::color::Bold);
}

std::string ColorNote(std::string s)
{
   return ConditionalColorText(s, terminal::color::Cyan, terminal::color::Bold);
}

std::string ColorHint(std::string s)
{
   return ConditionalColorText(s, terminal::color::Blue, terminal::color::Bold);
}

std::string ColorPrompt(std::string s)
{
   return ConditionalColorText(s, terminal::color::Green, terminal::color::Bold);
}

std::string ColorQuote(std::string s)
{
   return "'" + ColorHighlight(s) + "'";
}

std::string Spaces(int Size)
{
   return std::string(Size, ' ');
}

//
// ParserError
//

ParserError::ParserError(std::string const& Why)
   : CallStack(1, Why), Pos(NULL), End(NULL)
{
   this->AssembleMessage();
}

ParserError::ParserError(ParserError const& Prev, std::string const& Why)
   : CallStack(Prev.CallStack), Hint(Prev.Hint), Pos(Prev.Pos), End(Prev.End)
{
   CallStack.push_back(Why);
   this->AssembleMessage();
}

ParserError::ParserError(std::exception const& Prev, std::string const& Why)
   : CallStack(1, Prev.what()), Pos(NULL), End(NULL)
{
   CallStack.push_back(Why);
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_, char const* Position,
                         std::string const& Hint_)
   : CallStack(CallStack_), Hint(Hint_), Pos(Position), End(NULL)
{
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_, char const* Position,
                         char const* End_, std::string const& Hint_)
   : CallStack(CallStack_), Hint(Hint_), Pos(Position), End(End_)
{
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_,
                         std::string const& Why, char const* Position, char const* End_,
                         char const* beg, char const* end,
                         std::string const& Hint_)
   : CallStack(CallStack_), Hint(Hint_), Pos(NULL), End(NULL)
{
   std::string Next = Why + "\n" + std::string(beg, end);
   if (Position)
   {
      int Count = End_ ? std::distance(Position, End_) : 1;
      Next = Next + "\n" + std::string(std::distance(beg, Position), ' ')
         + ColorPrompt(std::string(Count, '^'));
   }
   CallStack.push_back(Next);
   this->AssembleMessage();
}

ParserError
ParserError::AtPosition(std::string const& Why, char const* Position)
{
   return ParserError(std::list<std::string>(1, Why), Position);
}

ParserError
ParserError::AtPosition(ParserError const& Prev, char const* Position)
{
   return ParserError(Prev.CallStack, Position, Prev.hint());
}

ParserError
ParserError::AtPosition(std::exception const& Prev, char const* Position)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Position);
}

ParserError
ParserError::AtRange(std::string const& Why, char const* Start, char const* End)
{
   return ParserError(std::list<std::string>(1, Why), Start, End);
}

ParserError
ParserError::AtRange(ParserError const& Prev, char const* Start, char const* End)
{
   return ParserError(Prev.CallStack, Start, End, Prev.hint());
}

ParserError
ParserError::AtRange(std::exception const& Prev, char const* Start, char const* End)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Start, End);
}

ParserError
ParserError::Finalize(ParserError const& Prev, std::string const& Why, char const* beg, char const* end)
{
   return ParserError(Prev.CallStack, Why, Prev.Pos, Prev.End, beg, end, Prev.hint());
}

ParserError
ParserError::Finalize(std::exception const& Prev, std::string const& Why, char const* beg, char const* end)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Why, NULL, NULL, beg, end);
}

void
ParserError::AssembleMessage()
{
   Msg = ColorError("Parser error: ");
   for (std::list<std::string>::const_iterator I = CallStack.begin(); I != CallStack.end(); ++I)
   {
      if (I == CallStack.begin())
      {
         Msg = Msg + (*I);
      }
      else
      {
         Msg = Msg + '\n' + ColorNote("note: ") + (*I);
      }
   }
   if (!Hint.empty())
   {
      Msg += '\n' + ColorHint("hint: ") + Hint;
   }
}

// LookupFileGrid
std::map<std::string, std::vector<std::vector<std::vector<std::string>>>> FileGridCache;

std::vector<std::string>
SplitLineIntoStrings(std::string line)
{
   std::vector<std::string> Result;

   std::string::const_iterator I = line.begin();
   while (I != line.end())
   {
      if (std::isspace(*I))
      {
         ++I;
         continue;
      }

      // we have a non-whitespace character.  This must be the start of a string.
      std::string Word;
      bool InQuote = false;
      while (I != line.end() && (InQuote || !std::isspace(*I)))
      {
         // A quoted string?
         if (*I == '"')
         {
            // skip over the quote
            ++I;
            // Toggle the InQuote flag
            InQuote = !InQuote;
            continue;
         }

         // An escape character? If so, treat the next character verbatim, even if it is whitespace or a quote.
         // Unless it is at the end of the string, in which case we ignore it.
         if (!InQuote && *I == '\\')
         {
            ++I;
            if (I == line.end())
               continue;
         }
         // add the character to the current word
         Word += *I;
         ++I;
      }
      Result.emplace_back(std::move(Word));
   }
   return Result;
}

void
LoadFileGrid(std::string Filename)
{
   std::vector<std::vector<std::vector<std::string>>> Lines;
   std::ifstream In(Filename);
   if (!In.good())
      throw ParserError("filegrid: file \"" + Filename + "\" does not exist!");
   std::string line;
   int z = 0;
   Lines.push_back(std::vector<std::vector<std::string>>());
   while (std::getline(In, line))
   {
      if (line == "")
      {
         // empty line = next block
         ++z;
         Lines.push_back(std::vector<std::vector<std::string>>());
         continue;
      }
      Lines[z].emplace_back(SplitLineIntoStrings(line));
   }
   FileGridCache.emplace(Filename, std::move(Lines));
}

std::string LookupFileGrid(std::string Filename, int x, int y, int z)
{
   auto I = FileGridCache.find(Filename);
   while (I == FileGridCache.end())
   {
      LoadFileGrid(Filename);
      I = FileGridCache.find(Filename);
   }

   if (z < 0 || z >= I->second.size())
   {
      throw ParserError("filegrid: z coordinate " + std::to_string(z) + " is out of range [0," + std::to_string(I->second.size()) + "]");
   }
   if (x < 0 || x >= I->second[z].size())
   {
      throw ParserError("filegrid: x coordinate " + std::to_string(x) + " is out of range [0," + std::to_string(I->second[z].size()) + "]");
   }
   if (y < 0 || y >= I->second[z][x].size())
   {
      throw ParserError("filegrid: y coordinate " + std::to_string(y) + " is out of range [0," + std::to_string(I->second[z][x].size()) + "]");
   }

   // Note the order of indices, to match the order used in the FileGridCache
   return I->second[z][x][y];
}

} // namespace Parser
