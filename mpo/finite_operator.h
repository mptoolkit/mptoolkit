// -*- C++ -*- $Id$
//
// FiniteOperator: an MPO that acts on a finite section of lattice.
// 

#if !defined(FINITE_OPERATOR_H_JDCHJKEHY589758YUER89H489)
#define FINITE_OPERATOR_H_JDCHJKEHY589758YUER89H489

#include "linear_operator.h"

class FiniteOperator : public LinearOperator
{
      FiniteOperator() {}

      // returns true if the operator transforms irreducibly.  This is true
      // iff the left basis contains only a single state.
      bool is_irreducible() const;

      // precondition: is_irreducible
      QuantumNumbers::QuantumNumber TransformsAs() const;
};

PStream::opstream& operator<<(PStream::opstream& out, FiniteOperator const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, FiniteOperator& op);

FiniteOperator& operator*=(FiniteOperator& x, double a);
FiniteOperator& operator*=(FiniteOperator& x, std::complex<double> a);

FiniteOperator& operator+=(FiniteOperator& x, FiniteOperator const& y);
FiniteOperator& operator-=(FiniteOperator& x, FiniteOperator const& y);

FiniteOperator operator+(FiniteOperator const& x, FiniteOperator const& y);
FiniteOperator operator-(FiniteOperator const& x, FiniteOperator const& y);

FiniteOperator operator-(FiniteOperator const& x);

FiniteOperator operator*(double a, FiniteOperator const& x);
FiniteOperator operator*(FiniteOperator const& x, double a);
FiniteOperator operator*(std::complex<double> a, FiniteOperator const& x);
FiniteOperator operator*(FiniteOperator const& x, std::complex<double> a);

FiniteOperator prod(FiniteOperator const& x, FiniteOperator const& y, QuantumNumbers::QuantumNumber const& q);
FiniteOperator prod(FiniteOperator const& x, FiniteOperator const& y);
FiniteOperator operator*(FiniteOperator const& x, FiniteOperator const& y);

// dot product - takes into account the multiplicity to rescale the result
FiniteOperator dot(FiniteOperator const& x, FiniteOperator const& y);

// project a (reducible) quantum number onto an irreducible component
FiniteOperator project(FiniteOperator const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.  Only useful for n small!
FiniteOperator pow(FiniteOperator const& x, int n);

// Conjugate
FiniteOperator conj(FiniteOperator const& x);

// Adjoint
FiniteOperator adjoint(FiniteOperator const& x);

// output to a stream
inline
std::ostream& operator<<(std::ostream& out, FiniteOperator const& x)
{
   return out << x.data();
}

#endif
