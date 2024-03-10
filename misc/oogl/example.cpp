// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/oogl/example.cpp
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

#include "oogl.h"
#include "arrow.h"
#include <iostream>

using namespace OOGL;

double frand()
{
   return double(rand()) / RAND_MAX;
}

double const Pi = 3.1415927;

int main()
{

   // a List represents a list of objects
   List MyList;

   for (int i = 0; i < 1000; ++i)
   {
      // choose a color
      Color c1 = Color::HSV(frand() * 360.0, 1, 1);
      Color c2 = Color::HSV(frand() * 360.0, 1, 1);
      Color c3 = Color::HSV(frand() * 360.0, 1, 1);

      // Create an arrow
      ColorPolyhedra MyArrow = CreateColorArrow(0.9, 0.9, 0.2, 0.4, 8, c1, c2, c3);

      // Choose a location
      Vertex Location(frand() * 100, frand() * 100, 0);

      // Rotate the arrow, shift and translate
      double Phi = frand() * Pi;
      double Theta =  frand() * Pi * 2;

      MyArrow *= Transform::rot_y(Phi) *
         Transform::rot_z(Theta) * Transform::scale(0.8) * Transform::translate(Location);

      // Add the arrow to the list
      MyList.push_back(MyArrow);

   }

   // Write the OOGL data to standard output
   std::cout << MyList << '\n';

   return 0;
}
