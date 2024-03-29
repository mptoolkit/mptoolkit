// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/boosttest.cpp
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

#include <boost/utility/result_of.hpp>

#include <iostream>
#include <iomanip>

#include "common/trace.h"

using tracer::typeid_name;

int F(int, int) { return 0; };

struct f
{
   int operator()(int, int) const { return 1; };

   template <typename T> struct result;

   template <typename T, typename U>
   struct result<f(T const&,U const&)> { typedef float type; };
};

int main()
{
   TRACE(typeid_name(boost::result_of<f(int,int)>::type()));
}
