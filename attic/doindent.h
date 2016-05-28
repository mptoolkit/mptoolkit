// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/doindent.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
