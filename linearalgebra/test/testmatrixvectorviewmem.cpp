// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/matrix_utility.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M = random_matrix<double>(10, 15);

   TRACE(M);
   TRACE(vector_view(M));
   TRACE(vector_view(M)[0]);

   TRACE(vector_view(transpose(M))[0]);

   TRACE(vector_view(swap_sort_order(M))[0]);

   //typedef VectorView<SwapSortOrder<Matrix<double>&>::result_type>::result_type X;
   //TRACE(transform(vector_view(swap_sort_order(M)), NormFrob<X>()));

   Vector<Vector<double> > x = vector_view(swap_sort_order(M));

   TRACE(x);
   TRACE(transform(x, NormFrob<Vector<double> >()));
}

