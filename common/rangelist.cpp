// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/rangelist.cpp
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "rangelist.h"
#include <cmath>
#include <iostream>

RangeList::RangeList(std::string Str)
{
   try
   {
      int Sep1 = Str.find_first_of(":");
      Start = std::stod(Str.substr(0, Sep1));

      // If only a single number is specified.
      if (Sep1 < 0)
      {
         End = Start;
         Step = 0.0;
         Num = 1;
         List = std::vector<double>(1, Start);
         return;
      }
      // If a range of numbers is specfied.
      else
      {
         int Sep2 = Str.find_first_of(":,", Sep1+1);
         if (Sep2 < 0)
            throw std::runtime_error("Missing second separator in RangeList " + Str);

         End = std::stod(Str.substr(Sep1+1, Sep2-Sep1-1));
         
         // If the step is specified.
         if (Str.at(Sep2) == ':')
         {
            Step = std::stod(Str.substr(Sep2+1));
            // Add an epsilon to ensure the end point is included if intended.
            Num = floor((End-Start)/Step + 1e-8) + 1;

            // If the step does not evenly divide the distance, we update the endpoint.
            if (std::abs(End-Start - (Num-1)*Step) / std::abs(Step) > 1e-8)
               End = Start + (Num-1)*Step;
         }
         // If the number of points is specified.
         else
         {
            Num = std::stoi(Str.substr(Sep2+1));
            Step = (End-Start)/(Num-1);
         }
      }
   }
   catch (std::invalid_argument& e)
   {
      throw std::invalid_argument("Invalid number in RangeList " + Str);
   }

   List = std::vector<double>();
   for (int i = 0; i < Num; ++i)
   {
      List.push_back(Start * ((double) (Num-1-i) / (Num-1)) + End * ((double) i / (Num-1)));
   }
}
