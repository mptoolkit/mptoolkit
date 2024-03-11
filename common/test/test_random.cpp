// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/test/test_random.cpp
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
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

#include "common/randutil.h"
#include <iostream>

int main()
{
   // the 10000th consecutive invocation of a default-contructed std::mt19937 is required to produce the value 4123659995.

   for (int i = 0; i <  9999; ++i)
   {
      randutil::u_rand();
   }

   unsigned x = randutil::u_rand();
   std::cout << x << '\n';
}
