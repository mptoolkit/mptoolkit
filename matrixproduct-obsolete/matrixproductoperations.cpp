// -*- C++ -*- $Id$

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

