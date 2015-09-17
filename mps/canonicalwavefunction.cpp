// -*- C++

#include "canonicalwavefunction.h"
#include "mps/linearwavefunction.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "mp-algorithms/arnoldi.h"


PStream::VersionTag
CanonicalWavefunctionBase::VersionT(1);

QuantumNumbers::QuantumNumber
CanonicalWavefunctionBase::TransformsAs() const
{
   PRECONDITION(this->is_irreducible());
   return transform_targets(this->Basis1()[0], adjoint(this->Basis2()[0]))[0];
}

void
CanonicalWavefunctionBase::ReadStream(PStream::ipstream& in)
{
   PStream::VersionSentry Sentry(in, VersionT, in.read<int>());

   in >> Data;
   in >> Lambda;
}

void
CanonicalWavefunctionBase::WriteStream(PStream::opstream& out) const
{
   out << VersionT.default_version();

   out << Data;
   out << Lambda;
}

