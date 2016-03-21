// -*- C++ -*-

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
