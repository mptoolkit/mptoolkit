// -*- C++ -*-

#include "canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "mp-algorithms/arnoldi.h"


PStream::VersionTag
CanonicalWavefunctionBase::VersionT(2);

QuantumNumbers::QuantumNumber
CanonicalWavefunctionBase::TransformsAs() const
{
   PRECONDITION(this->is_irreducible());
   return transform_targets(this->Basis1()[0], adjoint(this->Basis2()[0]))[0];
}

std::vector<BasisList>
ExtractLocalBasis(CanonicalWavefunctionBase const& Psi)
{
   std::vector<BasisList> Result;
   for (CanonicalWavefunctionBase::const_mps_iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
      Result.push_back(I->LocalBasis());
   }
   return Result;
}

void
CanonicalWavefunctionBase::ReadStream(PStream::ipstream& in)
{
   PStream::VersionSentry Sentry(in, VersionT, in.read<int>());

   in >> Data;
   in >> Lambda;

   if (Sentry.version() == 1)
   {
      if (Data.empty())
      {
	 Basis1_ = VectorBasis();
	 Basis2_ = VectorBasis();
      }
      else
      {
	 Basis1_ = Data.front().lock()->Basis1();
	 Basis2_ = Data.back().lock()->Basis2();
      }
   }
   else
   {
      in >> Basis1_;
      in >> Basis2_;
   }
}

void
CanonicalWavefunctionBase::WriteStream(PStream::opstream& out) const
{
   out << VersionT.default_version();

   out << Data;
   out << Lambda;

   out << Basis1_;
   out << Basis2_;
}

