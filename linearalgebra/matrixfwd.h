/* -*- C++ -*- $Id$

  forward declarations for some matrix types that are used as temporaries.

  Created 2005-03-10 Ian McCulloch
*/

#if !defined(MATRIXFWD_H_FHURH34Y89Y98Y9OY89PQYW9)
#define MATRIXFWD_H_FHURH34Y89Y98Y9OY89PQYW9

#include "matrixinterface.h"

namespace LinearAlgebra
{

template <typename Scalar, typename Orientation = RowMajor>
class Matrix;

// the 'default' matrix type for matrix_abstract_dense

template <typename T>
struct make_matrix_from_abstract<T, matrix_abstract_dense>
{
   typedef Matrix<T> type;
};

template <typename T, typename Orientation = RowMajor, 
          typename InnerType = MapVector<T>, typename OuterType = Vector<InnerType> >
class SparseMatrix;

// the default matrix type for matrix_abstract_sparse

template <typename T>
struct make_matrix_from_abstract<T, matrix_abstract_sparse>
{
   typedef SparseMatrix<T> type;
};

} // namespace LinearAlgebra

#endif
