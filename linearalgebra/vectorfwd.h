/* -*- C++ -*- $Id$

  vectorfwd.h

  Forward declarations so we can define the default value_types of expressions
  in vectorinterface.h
  
  Created 2005-01-12 Ian McCulloch
*/

#if !defined(VECTORFWD_H_CJKHUTYIUTYT897YT9874O)
#define VECTORFWD_H_CJKHUTYIUTYT897YT9874O

namespace LinearAlgebra
{

// dense vectors default to Vector
template <typename Scalar>
class Vector;

// sparse vectors default to Vector
template <typename T>
class MapVector;

} // namespace LinearAlgebra

#endif
