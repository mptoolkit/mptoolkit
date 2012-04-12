// -*- C++ -*- $Id$

#include "matrixproduct.h"

template <typename Component>
MatrixProduct<Component>
fancy_sum(MatrixProduct<Component> A, MatrixProduct<Component> B, int NumStates)
{
   typedef typename Component::OperatorType OperatorType;

   if (A.is_null()) return B;
   if (B.is_null()) return A;

   QuantumNumber Ident(A.GetSymmetryList());  // the scalar quantum number

   // We start from the right hand side here...
   while (A.RightSize() > 1) A.RotateRight();
   while (B.RightSize() > 1) B.RotateRight();
   
