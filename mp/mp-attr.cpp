// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-attr.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "quantumnumbers/all_symmetries.h"
#include "wavefunction/mpwavefunction.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include <boost/algorithm/string.hpp>

int main(int argc, char** argv)
{
   if (argc <= 1)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-attr wavefunction [attribute [= value]] ...\n";
      return 1;
   }

   typedef std::map<std::string, std::string> AttribToSetType;
   typedef AttribToSetType::const_iterator AttribSetIter;
   std::map<std::string, std::string> AttributesToSet;

   typedef std::vector<std::string>::const_iterator AttribPrintIter;
   std::vector<std::string> AttributesToPrint;

   int arg = 2;
   while (arg < argc)
   {
      std::string a = argv[arg++];
      std::string Attrib;
      std::string::const_iterator Delim = std::find(a.begin(), a.end(), '=');
      if (Delim == a.end() && (arg >= argc || argv[arg][0] != '='))
      {
         Attrib = a;
         boost::trim(Attrib);
      }
      else
      {
         // attribute = value
         Attrib = std::string(static_cast<std::string const&>(a).begin(), Delim);
         boost::trim(Attrib);
         if (Delim == a.end())
         {
            // in this case, the first character of the next argument must be '='
            a = argv[arg++];
            Delim = a.begin();
         }
         ++Delim;
         if (Delim == a.end() && arg < argc)
         {
            a = argv[arg++];
            Delim = a.begin();
         }
         std::string Value(Delim, static_cast<std::string const&>(a).end());
         boost::trim(Value);
         AttributesToSet[Attrib] = Value;
      }
      AttributesToPrint.push_back(Attrib);
   }

   bool Readonly = AttributesToSet.empty();
   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], mp_pheap::CacheSize(),
                                                          Readonly);

   if (!AttributesToSet.empty())
   {
      pvalue_lock<MPWavefunction> PsiLock(Psi.lock());
      for (AttribSetIter I = AttributesToSet.begin(); I != AttributesToSet.end(); ++I)
         PsiLock->Attributes()[I->first] = I->second;
   }

   if (AttributesToPrint.size() == 0)
   {
      std::cout << Psi->Attributes();
   }

   // if there is only one attribute to print, show it without the Attrib=
   if (AttributesToSet.empty() && AttributesToPrint.size() == 1)
   {
      std::cout << Psi->Attributes()[AttributesToPrint[0]] << '\n';
   }
   else
   {
      for (AttribPrintIter I = AttributesToPrint.begin(); I != AttributesToPrint.end(); ++I)
      {
         std::cout << (*I) << "=" << Psi->Attributes()[*I] << '\n';
      }
   }

   if (!Readonly)
   {
      pheap::ShutdownPersistent(Psi);
   }
}
