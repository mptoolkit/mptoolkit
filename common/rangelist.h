// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/rangelist.h
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

/*

RangeList is a class which parses a string specifying a linearly spaced range of numbers.

The allowed input strings are:

1. A single number (in this case Start == End == the input, Num == 1 and Step == 0.0).

2. A range specified as Start:End:Step.
   (If Step does not divide (End-Start) exactly, the final number will be the
   greatest value (Start + i*Step) for integer i which is less than End.)

3. A range specified as Start:End,Num.
*/

#if !defined(MPTOOLKIT_COMMON_RANGELIST_H)
#define MPTOOLKIT_COMMON_RANGELIST_H

#include <string>
#include <vector>

class RangeList
{
   public:
      typedef std::vector<double>::const_iterator const_iterator;

      RangeList() {}

      explicit RangeList(std::string Str);
      
      //std::vector<double> get_list() const { return List; }
      double get_start() const { return Start; }
      double get_end() const { return End; }
      double get_step() const { return Step; }
      int get_num() const { return Num; }

      const_iterator begin() const { return List.cbegin(); }
      const_iterator end() const { return List.cend(); }

   private:
      std::vector<double> List;
      double Start;
      double End;
      double Step;
      int Num;
};

#endif
