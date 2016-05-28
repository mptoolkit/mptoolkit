// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/denseoperator.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "denseoperator.h"

std::ostream& operator<<(std::ostream& out, DenseOperator const& x)
{
   out << "DenseOperator:\nBasis1() = " << x.Basis1()
       << "\nBasis2() = " << x.Basis2()
       << "data() =\n" << x.data();
   return out;
}

   
