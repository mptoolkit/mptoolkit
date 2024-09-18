// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/doindent.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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

std::string indent_typeid(std::string In)
{
   std::ostringstream Out;
   std::string::const_iterator I = In.begin();
   int Indent = 0;
   while (I != In.end())
   {
      if (*I == '<')
      {
         ++Indent;
         Out << '\n';
         for (int i = 0; i < Indent*3; ++i) Out << ' ';
         Out << '<';
      }
      else if (*I == '>')
      {
         --Indent;
         Out << '\n';
         for (int i = 0; i < Indent*3; ++i) Out << ' ';
         out << ">";
      }
      else if (*I == ',')
      {
         Out << '\n';
         for (int i = 0; i < (Indent-1)*3; ++i) Out << ' ';
         Out << " , ";
      }
   }
   return Out.str();
}
