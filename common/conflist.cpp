// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/conflist.cpp
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

#include "conflist.h"
#include "environment.h"
#include "common/stringutil.h"
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdlib.h>

std::string ExpandEnvironment(std::string const& s)
{
   std::string result;
   std::string::const_iterator I = std::find(s.begin(), s.end(), '$');
   std::string::const_iterator marker = s.begin();
   while (I != s.end())
   {
      result += std::string(marker, I);
      marker = I;

      // see if we have a valid environment string
      if (++I != s.end() && *I == '{')
      {
         ++I;
         std::string::const_iterator EndEnv = std::find(I, s.end(), '}');
         if (EndEnv != s.end())
         {
            std::string EnvString(I, EndEnv);
            char* Replacement = getenv(EnvString.c_str());
            if (Replacement)
            {
               ++EndEnv;
               marker = I = EndEnv;
               result += std::string(Replacement);
            }
         }
      }

      I = std::find(I, s.end(), '$');
   }
   result += std::string(marker, I);
   return result;
}

ConfList::ConfList(std::string const& file)
{
   std::ifstream InFile(file.c_str());
   if (!InFile.good()) throw std::runtime_error("Cannot open configuration file " + file + "!");
   std::string Str;
   std::getline(InFile, Str);
   while (InFile)
   {
      if (!Str.empty())
      {
         // remove the comment, if any
         Str = Str.substr(0, Str.find('!'));
         Str = Str.substr(0, Str.find('#'));

         std::size_t Pos = Str.find('=');
         if (Pos != std::string::npos && Pos > 0)
         {
            std::string Name = Str.substr(0, Pos);
            std::string Value = Str.substr(Pos+1);

            RemoveWhiteSpace(Name);
            RemoveWhiteSpace(Value);
            Value = ExpandEnvironment(Value);
            if (Name.length() > 0)
            {
               if (Data.find(Name) != Data.end())
               {
                  std::cerr << "warning: configuration file " << file << ": option "
                            << Name << " is redefined (was " << Data[Name]
                            << ", now " << Value << ")\n";
               }
               Data[Name] = Value;
            }
         }
      }

      std::getline(InFile, Str);
   }
}

std::string ConfList::Get(const std::string& s, const std::string& Default) const
{
   std::map<std::string, std::string>::const_iterator I = Data.find(s);
   if (I == Data.end()) return Default;

   return I->second;
}

std::string ConfList::Get(const std::string& s, char const* Default) const
{
   std::map<std::string, std::string>::const_iterator I = Data.find(s);
   if (I == Data.end()) return Default;

   return I->second;
}

bool ConfList::Get(const std::string& s, bool Default) const
{
   std::map<std::string, std::string>::const_iterator I = Data.find(s);
   if (I == Data.end()) return Default;

   std::string Str = I->second;

   std::transform(Str.begin(), Str.end(), Str.begin(), toupper);

   if (Str == "1" || Str == "YES" || Str == "TRUE") return true;
   if (Str == "0" || Str == "NO" || Str == "FALSE") return false;

   return Default;
}

unsigned long
ConfList::GetBytes(std::string const& s, unsigned long Default) const
{
   std::map<std::string, std::string>::const_iterator I = Data.find(s);
   if (I == Data.end()) return Default;

   std::string Str = I->second;
   std::istringstream In(Str);
   unsigned long Value = Default;
   In >> Value;
   std::string Suffix;
   In >> Suffix;
   std::transform(Suffix.begin(), Suffix.end(), Suffix.begin(), toupper);

   if (Suffix == "KB" || Suffix == "K") Value *= 1024;
   else if (Suffix == "MB" || Suffix == "M") Value *= 1024UL * 1024UL;
   else if (Suffix == "GB" || Suffix == "G") Value *= 1024UL * 1024UL * 1024UL;
   else if (Suffix == "TB" || Suffix == "T") Value *= 1024UL * 1024UL * 1024UL * 1024UL;

   return Value;
}
