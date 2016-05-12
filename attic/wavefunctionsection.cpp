// -*- C++ -*-

#include "wavefunctionsection.h"
#include "wavefunction/linearwavefunction.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "mp-algorithms/arnoldi.h"

// Serialization versions:
// Version 1:
// std::vector<pvalue_handle<StateComponent> >       Data
// std::vector<pvalue_handle<RealDiagonalOperator> > Lambda
// VectorBasis                                       Basis1
// VectorBasis                                       Basis2

PStream::VersionTag
WavefunctionSection::VersionT(1);

WavefunctionSection::WavefunctionSection(WavefunctionSection const& Psi) 
   : Data(Psi.Data), Lambda(Psi.Lambda), 
     Basis1_(Psi.Basis1_), Basis2_(Psi.Basis2_) 
{
}

WavefunctionSection::WavefunctionSection(CanonicalWavefunctionBase const& Psi)
   : Data(Psi.base_begin(), Psi.base_end()),
     Basis1_(Psi.Basis1()), Basis2_(Psi.Basis2())
{
   if (Psi.empty())
      return;

   auto i = Psi.lambda_base_begin();
   ++i;

   auto j = Psi.lambda_base_end();
   --i;

   Lambda = lambda_container_type(i, j);
}

void
WavefunctionSection::check_structure() const
{
   if (this->empty())
      return;

   CHECK_EQUAL(Lambda.size(), Data.size()-1);

   for (int i = 1; i < this->size()-1; ++i)
   {
      CHECK_EQUAL(this->lambda(i).Basis1(), this->operator[](i-1).Basis2())(i);
      CHECK_EQUAL(this->lambda(i).Basis2(), this->operator[](i).Basis1())(i);
      this->lambda(i).check_structure();
      this->operator[](i).check_structure();
   }
   this->operator[](0).check_structure();
   CHECK_EQUAL(Basis1_, this->operator[](0).Basis1());
   CHECK_EQUAL(Basis2_, this->operator[](this->size()-1).Basis1());
}

WavefunctionSection&
WavefunctionSection::operator=(WavefunctionSection const& Psi)
{ 
   Data = Psi.Data;
   Lambda = Psi.Lambda;
   Basis1_ = Psi.Basis1_;
   Basis2_ = Psi.Basis2_;
   return *this; 
}

std::vector<BasisList>
ExtractLocalBasis(WavefunctionSection const& Psi)
{
   std::vector<BasisList> Result;
   for (WavefunctionSection::const_mps_iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
      Result.push_back(I->LocalBasis());
   }
   return Result;
}

int
WavefunctionSection::ReadStream(PStream::ipstream& in)
{
   PStream::VersionSentry Sentry(in, VersionT, in.read<int>());

   CHECK(Sentry.version() >= 1);
   CHECK(Sentry.version() <= 1)
      ("Unrecognised version number in WavefunctionSection - this software is too old to read this file!")
      (Sentry.version());

   in >> Data;
   in >> Lambda;
   in >> Basis1_;
   in >> Basis2_;

   this->check_structure();

   return Sentry.version();
}

void
WavefunctionSection::WriteStream(PStream::opstream& out) const
{
   out << VersionT.default_version();

   out << Data;
   out << Lambda;

   out << Basis1_;
   out << Basis2_;
}
