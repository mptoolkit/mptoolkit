// -*- C++ -*- $Id$

#include "dataops.h"

namespace LinearAlgebra
{

// VectorExpression

// no need to explicitly call cow() here, the begin() function should do it for us.

template <class Scalar, class Derived>
template <class S2, class D2>
inline
void 
VectorExpression<Scalar, Derived>::assign(VectorConstExpression<S2, D2> const& V)
{
  PRECONDITION(V.size() == size());
  ops::fast_copy(V.begin(), V.end(), this->begin());
}

template <class Scalar, class Derived>
template <class S2, class D2>
inline
void 
VectorExpression<Scalar, Derived>::add(VectorConstExpression<S2, D2> const& V)
{
  PRECONDITION(V.size() == size());
  ops::fast_add(V.begin(), V.end(), this->begin());
}

template <class Scalar, class Derived>
template <class S2, class D2>
inline
void 
VectorExpression<Scalar, Derived>::subtract(VectorConstExpression<S2, D2> const& V)
{
  PRECONDITION(V.size() == size());
  ops::fast_subtract(V.begin(), V.end(), this->begin());
}

template <class Scalar, class Derived>
template <typename S>
inline
void 
VectorExpression<Scalar, Derived>::multiply(S const& s)
{
  ops::fast_multiply(this->begin(), this->end(), s);
}

template <class Scalar, class Derived>
template <typename S>
inline
void 
VectorExpression<Scalar, Derived>::fill(S const& s)
{
  ops::fast_fill(this->begin(), this->end(), s);
}

} // namespace LinearAlgebra
