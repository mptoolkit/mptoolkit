// -*- C++ -*- $Id$

#include "linear_operator.h"
#include "pstream/variant.h"
#include <algorithm>

using QuantumNumbers::QuantumNumber;

bool LinearOperator::is_irreducible() const
{
   //   return (!this->is_null() && this->Basis1().size() == 1);
   return (this->is_null() || this->Basis1().size() == 1);
}

QuantumNumbers::QuantumNumber 
LinearOperator::TransformsAs() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   CHECK(this->is_irreducible());
   return this->Basis1()[0];
}

PStream::opstream& operator<<(PStream::opstream& out, LinearOperator const& op)
{
   return out << op.data();
}

PStream::ipstream& operator>>(PStream::ipstream& in, LinearOperator& op)
{
   return in >> op.data();
}

