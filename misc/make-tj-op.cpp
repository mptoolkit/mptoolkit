// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/make-tj-op.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include <boost/lexical_cast.hpp>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      std::cerr << "usage: make-tj-op <L> <HoppingRange> <SpinRange> <CoulombRange>\n";
      exit(1);
   }

   int L = boost::lexical_cast<int>(argv[1]);
   int HopR = boost::lexical_cast<int>(argv[2]);
   int SpinR = boost::lexical_cast<int>(argv[3]);
   int CoulombR = boost::lexical_cast<int>(argv[4]);

   for (int i = 1; i <= L; ++i)
   {
      for (int j = i+1; j <= L; ++j)
      {
         if (j-i <= HopR)
            std::cout <<  " PairHopg("<<i<<","<<j<<")";
         if (j-i <= SpinR)
            std::cout <<  " SpinHop("<<i<<","<<j<<")";
         if (j-i <= CoulombR)
            std::cout <<  " VgHop("<<i<<","<<j<<")";
      }
   }
}
