// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/matrixproductoperations.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "matrixproduct.h"
#include "mpoperator.h"
#include "matrixproductoperator.h"

MPStateComponent
reduced_matrix_elements(MPWavefunction const& A, MPOperator const& M, MPWavefunction const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.GetSymmetryList(), B.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(A.GetSymmetryList(), M.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(A.size(), B.size());
   DEBUG_PRECONDITION_EQUAL(A.size(), M.size());

   //CHECK_EQUAL(A.Basis2(), B.Basis2());

   VectorBasis Vacuum = VectorBasis(make_vacuum_basis(A.GetSymmetryList()));
   MPStateComponent Elem(make_vacuum_basis(A.GetSymmetryList()), Vacuum, Vacuum);
   Elem[0](0,0) = LinearAlgebra::identity_matrix<double>(1);

   for (int Loc = A.size()-1; Loc >= 0; --Loc)
   {
      Elem = operator_prod(M[Loc], B.LookupLinear(Loc), Elem, herm(A.LookupLinear(Loc)));
   }

   return Elem;
}

MPStateComponent
reduced_matrix_elements(MPWavefunction const& A, MPOperator const& M, MPWavefunction const& B,
                        MPStateComponent Elem)
{
   DEBUG_PRECONDITION_EQUAL(A.GetSymmetryList(), B.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(A.GetSymmetryList(), M.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(A.size(), B.size());
   DEBUG_PRECONDITION_EQUAL(A.size(), M.size());

   for (int Loc = A.size()-1; Loc >= 0; --Loc)
   {
      Elem = operator_prod(M[Loc], B.LookupLinear(Loc), Elem, herm(A.LookupLinear(Loc)));
   }

   return Elem;
}

std::complex<double>
expectation(MPWavefunction const& A, MPOperator const& M, MPWavefunction const& B)
{
   CHECK_EQUAL(A.TransformsAs(), B.TransformsAs());
   return trace(reduced_matrix_elements(A,M,B)[0]);
}
