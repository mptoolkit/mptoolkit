// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testslice.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/vector.h"
#include "linearalgebra/slice.h"

int main()
{
   LinearAlgebra::Slice s(0,5,4);
   LinearAlgebra::Vector<double> v1(s);
   std::cout << s << std::endl;
   std::cout << v1 << std::endl;
}
