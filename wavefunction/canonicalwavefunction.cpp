// -*- C++ -*-

#include "canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "mp-algorithms/arnoldi.h"

// Serialization versions:
// Version 1:
// std::vector<pvalue_handle<StateComponent> >           Data
// std::vector<pvalue_handle<RealSemiDiagonalOperator> > Lambda
//
// Vesion 2:
// std::vector<pvalue_handle<StateComponent> >           Data
// std::vector<pvalue_handle<RealSemiDiagonalOperator> > Lambda
// VectorBasis                                           Basis1
// VectorBasis                                           Basis2
//
// Version 3:
// Same as version 2, but the Lambda array has size()+1 components, with
// Lambda(0) being the same as Lambda(size()) [up to a quantum number shift].
// Note that the CanonicalWavefunctionBase stream operator cannot reconstruct
// Lambda(0) on its own, it sets Lambda(0) to be a dummy.  The derived class
// serialization must handle this initialization itself, based on the returned
// version number from CanonicalWavefunctionBase::ReadStream().
//
// Version 4:
// Same as version 3, but with RealDiagonalOperator
// std::vector<pvalue_handle<StateComponent> >       Data
// std::vector<pvalue_handle<RealDiagonalOperator> > Lambda
// VectorBasis                                       Basis1
// VectorBasis                                       Basis2


PStream::VersionTag
CanonicalWavefunctionBase::VersionT(4);

CanonicalWavefunctionBase::CanonicalWavefunctionBase(CanonicalWavefunctionBase const& Psi) 
   : Data(Psi.Data), Lambda(Psi.Lambda), 
     Basis1_(Psi.Basis1_), Basis2_(Psi.Basis2_) 
{
}

void
CanonicalWavefunctionBase::check_structure() const
{
   if (this->empty())
      return;

   CHECK_EQUAL(Lambda.size(), Data.size()+1);

   for (int i = 0; i < this->size(); ++i)
   {
      CHECK_EQUAL(this->lambda(i).Basis2(), this->operator[](i).Basis1())(i);
      CHECK_EQUAL(this->operator[](i).Basis2(), this->lambda(i+1).Basis1())(i);
      this->lambda(i).check_structure();
      this->operator[](i).check_structure();
   }
   this->lambda(this->size()).check_structure();
   CHECK_EQUAL(Basis1_, this->lambda(0).Basis1());
   CHECK_EQUAL(Basis2_, this->lambda(this->size()).Basis2());
}

CanonicalWavefunctionBase&
CanonicalWavefunctionBase::operator=(CanonicalWavefunctionBase const& Psi)
{ 
   Data = Psi.Data;
   Lambda = Psi.Lambda;
   Basis1_ = Psi.Basis1_;
   Basis2_ = Psi.Basis2_;
   return *this; 
}

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

int
CanonicalWavefunctionBase::ReadStream(PStream::ipstream& in)
{
   PStream::VersionSentry Sentry(in, VersionT, in.read<int>());

   CHECK(Sentry.version() >= 1);
   CHECK(Sentry.version() <= 4)
      ("Unrecognised version number in CanonicalWavefunctionBase - this software is too old to read this file!")
      (Sentry.version());

   in >> Data;

   if (Sentry.version() < 4)
   {
      using old_lambda_type           = RealSemiDiagonalOperator;
      using old_lambda_handle_type    = pvalue_handle<old_lambda_type>;
      using old_lambda_container_type = std::vector<old_lambda_handle_type>;

      // if we're reading version 3, then Lambda0 isn't yet part of the serialization format.
      // Add a dummy value for now.
      if (Sentry.version() < 3)
      {
	 Lambda.clear();
	 old_lambda_container_type LambdaTemp;
	 in >> LambdaTemp;

	 // argh - except for a bug where we forgot to increment the version number on saving.  So
	 // hack around it!
	 if (LambdaTemp.size() == Data.size()+1)
	 {
	    // this means we forgot to increment the version number
	    Sentry.change_version(3);
	 }
	 else
	 {
	    this->push_back_lambda(RealDiagonalOperator());
	 }
	 for (unsigned i = 0; i < LambdaTemp.size(); ++i)
	 {
	    Lambda.push_back(new RealDiagonalOperator(*LambdaTemp[i].lock()));
	 }
      }
      else
      {
	 Lambda.clear();
	 old_lambda_container_type LambdaTemp;
	 in >> LambdaTemp;
	 for (unsigned i = 0; i < LambdaTemp.size(); ++i)
	 {
	    Lambda.push_back(new RealDiagonalOperator(*LambdaTemp[i].lock()));
	 }
      }
   }
   else
   {
      // version 4
      in >> Lambda;
   }

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

   CHECK_EQUAL(Lambda.size(), Data.size()+1);

   return Sentry.version();
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

