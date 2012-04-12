// -*- C++ -*- $Id$

SimpleOperator
reduced_matrix_elements(MPWavefunction const& A, MPOperator const& M, MPWavefunction const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.GetSymmetryList(), B.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(A.GetSymmetryList(), M.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(A.size(), B.size());
   DEBUG_PRECONDITION_EQUAL(A.size(), M.size());

   CHECK_EQUAL(A.Basis2(), B.Basis2());

   SimpleOperator Elem = SimpleOperator::construct_identity(A.Basis2());

   for (int Loc = A.size()-1; Loc >= 0; --Loc)
   {
      Elem = operator_prod(M[Loc], B.LookupLinear(Loc), Elem, herm(A.LookupLinear(Loc)));
   }

   return Elem;
}

SimpleOperator
reduced_matrix_elements(MPWavefunction const& A, MPOperator const& M, MPWavefunction const& B,
                        SimpleOperator Elem)
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
   reutrn reduced_matrix_elements(A,M,B)(0,0);
}

