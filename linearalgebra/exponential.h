// -*- C++ -*- $id$
//
// Exponential function for matrices.

#if !defined(EXPONENTIAL_H_SDJKCH48935Y78Y7RFEH89P)
#define EXPONENTIAL_H_SDJKCH48935Y78Y7RFEH89P

#include "matrix.h"

namespace LinearAlgebra
{


//
// Exponentiate
//
// Calculate exp(M) of a matrix M, using Expokit.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementExponentiate {};

template <typename M>
inline
typename ImplementExponentiate<M>::result_type
Exponentiate(double t, M const& m)
{
   return ImplementExponentiate<M>()(t, m);
}

} // namespace LinearAlgebra

#include "exponential.cc"

#endif
