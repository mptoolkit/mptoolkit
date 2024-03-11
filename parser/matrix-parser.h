// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// parser/matrix-parser.h
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

// Parser for matrix expressions

// The element type is a MatrixOperator

// we also accept usual arguments, and also take a UnitCell so that
// we can use local functions (do we need this?)

#if !defined(MPTOOLKIT_MATRIX_PARSER_H)
#define MPTOOLKIT_MATRIX_PARSER_H

#include "mps/state_component.h"
#include "lattice/function.h"
#include <string>
#include <map>

MatrixOperator ParseMatrixOperator(std::string const& Str,
                                   Function::ArgumentList const& Args,
                                   std::map<std::string, MatrixOperator> const& Matrices);

#endif
