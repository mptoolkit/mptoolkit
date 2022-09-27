// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/infinitewavefunctionleft.cpp
//
// Copyright (C) 2015-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "infinitewavefunctionleft.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "common/environment.h"

#include "common/statistics.h"

#include <fstream>

#include "wavefunction/operator_actions.h"
#include "pheap/pheapstream.h"
#include "interface/attributes.h"

// Streaming versions:
// Note: the base class CanonicalWavefunctionBase has a separate version number.
//
// Version 1:
// No verioning information, this version needs to be deduced from the stream metadata version.
// Old style InfiniteWavefunction:
//      MatrixOperator        C_old
//      QuantumNumber         QShift
//      LinearWavefunctionOld PsiLinear
//      MatrixOperator        C_right
//      AttributeList         Attr
//
// Version 2:
//      CanonicalWavefunctionBase (base class)
//      QuantumNumber              QShift
//
// Version 3:
//      double                     Amplitude
// Version 4:
//      remove 'Amplitude', replace by
//      double                     LogAmplitude

PStream::VersionTag
InfiniteWavefunctionLeft::VersionT(4);

extern double const ArnoldiTol = getenv_or_default("MP_ARNOLDI_TOL", 1E-15);

extern double const InverseTol = getenv_or_default("MP_INVERSE_TOL", 1E-7);

// the tol used in the orthogonalization can apparently be a bit smaller
extern double const OrthoTol = getenv_or_default("MP_ORTHO_TOL", 1E-8);

std::string InfiniteWavefunctionLeft::Type = "InfiniteWavefunctionLeft";

namespace
{

struct LeftMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiply(LinearWavefunction const& L_, QuantumNumber const& QShift_, int Verbose_ = 0)
      : L(L_), QShift(QShift_), Verbose(Verbose_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      int s = 0;
      for (LinearWavefunction::const_iterator I = L.begin(); I != L.end(); ++I)
      {
         if (Verbose > 0)
            std::cout << "site " << s << std::endl;
         r = operator_prod(herm(*I), r, *I);
         ++s;
      }
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   int Verbose;
};

struct RightMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiply(LinearWavefunction const& R_, QuantumNumber const& QShift_, int Verbose_ = 0)
      : R(R_), QShift(QShift_), Verbose(Verbose_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = x;
      int s = R.size();
      LinearWavefunction::const_iterator I = R.end();
      while (I != R.begin())
      {
         --s;
         if (Verbose > 0)
            std::cout << "site " << s << std::endl;
         --I;
         r = operator_prod(*I, r, herm(*I));
      }
      return delta_shift(r, adjoint(QShift));
   }

   LinearWavefunction const& R;
   QuantumNumber QShift;
   int Verbose;
};

} // namespace

InfiniteWavefunctionLeft::InfiniteWavefunctionLeft(QuantumNumbers::QuantumNumber const& QShift_, double LogAmplitude_)
   : QShift(QShift_), LogAmplitude(LogAmplitude_)
{
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::ConstructFromOrthogonal(LinearWavefunction const& Psi, MatrixOperator const& Lambda,
                                                  QuantumNumbers::QuantumNumber const& QShift_,
                                                  double LogAmplitude,
                                                  int Verbose)
{
   std::cerr << "warning: ConstructFromOrthogonal might not be reliable!\n";
   InfiniteWavefunctionLeft Result(QShift_, LogAmplitude);
   Result.Initialize(Psi, Lambda, Verbose-1);
   return Result;
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::ConstructFromOrthogonal(LinearWavefunction const& Psi, MatrixOperator const& Lambda,
                                                  QuantumNumbers::QuantumNumber const& QShift_,
                                                  MatrixOperator& Vh, double LogAmplitude,
                                                  int Verbose)
{
   InfiniteWavefunctionLeft Result(QShift_, LogAmplitude);
   Vh = Result.Initialize(Psi, Lambda, Verbose-1);
   return Result;
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::Construct(LinearWavefunction const& Psi,
                                    QuantumNumbers::QuantumNumber const& QShift_,
                                    double LogAmplitude,
                                    int Verbose)
{
   return InfiniteWavefunctionLeft::Construct(Psi, MatrixOperator::make_identity(Psi.Basis2()), QShift_, LogAmplitude, Verbose);
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::Construct(LinearWavefunction const& Psi, MatrixOperator const& GuessRho,
                                    QuantumNumbers::QuantumNumber const& QShift,
                                    double LogAmplitude,
                                    int Verbose)
{
   LinearWavefunction PsiL = Psi;

   MatrixOperator Guess = MatrixOperator::make_identity(PsiL.Basis2());

   // initialize LeftEigen to a guess eigenvector.  Since L satisfies the left orthogonality
   // constraint (except for the final matrix), we can do one iteration beyond the identity
   // and intialize it to herm(Xu) * Xu
   MatrixOperator LeftEigen = Guess;

   if (Verbose > 0)
      std::cout << "Obtaining left orthogonality eigenvector..." << std::endl;
   // get the eigenmatrix.  Do some dodgy explict restarts.
   int Iterations = 20;
   double Tol = ArnoldiTol;
   LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, QShift, Verbose-2),
                                                      Iterations, Tol,
                                                      LinearSolvers::LargestAlgebraicReal, false, Verbose-1);
   while (Tol < 0)
   {
      if (Verbose > 0)
         std::cout << "LeftEigen: Arnoldi not converged, restarting.  EValue="
                   << EtaL << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, QShift, Verbose-2), Iterations,
                                    Tol, LinearSolvers::LargestAlgebraicReal, false, Verbose-1);
   }

   CHECK(EtaL.real() > 0)("Eigenvalue must be positive")(EtaL);

   DEBUG_TRACE(EtaL);

   // Adjust the phase of the eigenvector to maximize the trace (ie, make it positive)
   std::complex<double> Alpha = std::conj(trace(LeftEigen));
   Alpha *= 1.0 / norm_frob(Alpha);
   LeftEigen *= Alpha;

   //TRACE(LeftEigen)(EigenvaluesHermitian(LeftEigen));

   DEBUG_CHECK(norm_frob(LeftEigen - adjoint(LeftEigen)) < 1e-12)
      (norm_frob(LeftEigen - adjoint(LeftEigen)));

   //   DEBUG_TRACE(EigenvaluesHermitian(RightEigen));
   //   DEBUG_TRACE(EigenvaluesHermitian(LeftEigen));

   // Do the Cholesky factorization of the eigenmatrices.
   // Actually a SVD is much more stable here
   //   MatrixOperator A = CholeskyFactorizeUpper(LeftEigen);
   //   MatrixOperator B = CholeskyFactorizeUpper(RightEigen);

   MatrixOperator D = LeftEigen;
   MatrixOperator U = DiagonalizeHermitian(D);
   D = SqrtDiagonal(D, OrthoTol);
   MatrixOperator DInv = InvertDiagonal(D, OrthoTol);

   // At this point, any matrix elements in D that are smaller than OrthoTol can be removed
   // from the basis, because they have negligible weight in the wavefunction.
   // However, we should be able to cope with states in the basis that have zero weight;
   // we can still orthogonalize the basis, and the matrix elements of any unitary that
   // act in the direction of small elements of D are arbitary, so we can set them to whatever is
   // needed to canonicalize the state.
   // If we are orthogonalizing the state as an intermediate step and expecting to continue calculations,
   // then we should keep all states even if they have small (or zero) weight.  But if this is the final calculation
   // then we should remove small elements as they have no effect on the final wavefunction.

   // LeftEigen = triple_prod(U, D*D, herm(U))
   DEBUG_CHECK(norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)))(D)(DInv);

   MatrixOperator A = delta_shift(D * U, QShift);
   MatrixOperator AInv = adjoint(U) * DInv;

#if 1
   // Explicitly re-othogonalize.  In principle this is not necessary,
   // but it gives some more robustness
   A = left_orthogonalize(A, PsiL);
   PsiL.set_back(prod(PsiL.get_back(), A*AInv));
#else
   PsiL.set_front(prod(A, PsiL.get_front()));
   PsiL.set_back(prod(PsiL.get_back(), AInv));
#endif

   // At this point, the left eigenvector is the identity matrix.
#if !defined(NDEBUG)
   MatrixOperator I = MatrixOperator::make_identity(PsiL.Basis1());
   MatrixOperator R = delta_shift(D*D, QShift);
   A = delta_shift(LeftMultiply(PsiL, QShift)(delta_shift(I, adjoint(QShift))), QShift);
   CHECK(norm_frob(inner_prod(A-EtaL*I, R)) < 10*A.Basis1().total_dimension() * ArnoldiTol)(norm_frob(A-EtaL*I))(A)(I)(D);
   TRACE(norm_frob(inner_prod(A-EtaL*I, R)));
#endif

   // same for the right eigenvector, which will be the density matrix

   // initialize the guess eigenvector
   MatrixOperator RightEigen = GuessRho;

   if (Verbose > 0)
      std::cout << "Obtaining right orthogonality eigenvector..." << std::endl;
   // get the eigenmatrix
   Iterations = 20; Tol = ArnoldiTol;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, QShift, Verbose-2),
                                                      Iterations, Tol,
                                                      LinearSolvers::LargestAlgebraicReal, false, Verbose-1);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Tol < 0)
   {
      if (Verbose > 0)
         std::cout << "RightEigen: Arnoldi not converged, restarting.  EValue="
                   << EtaR << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, QShift, Verbose-2),
                                    Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false, Verbose-1);
   }
   DEBUG_TRACE(EtaR);

   // Adjust the phase of the eigenvector to maximize the trace (ie, make it positive)
   Alpha = std::conj(trace(RightEigen));
   Alpha *= 1.0 / norm_frob(Alpha);
   RightEigen *= Alpha;

   //TRACE(trace(RightEigen));

   DEBUG_CHECK(norm_frob(RightEigen - adjoint(RightEigen)) < 1e-12)
      (norm_frob(RightEigen - adjoint(RightEigen)))(RightEigen);

   D = RightEigen;
   U = DiagonalizeHermitian(D);
   D = SqrtDiagonal(D, OrthoTol);

   // normalize
   D *= 1.0 / norm_frob(D);

#if 1
   // incorporate U into the MPS

   PsiL.set_back(prod(PsiL.get_back(), adjoint(U)));
   PsiL.set_front(prod(delta_shift(U, QShift), PsiL.get_front()));

   // D is now the Lambda matrix.  Check that
#if !defined(NDEBUG)
   MatrixOperator Rho = D*D;
   A = RightMultiply(PsiL, QShift)(Rho);
   CHECK(norm_frob(A-EtaR*Rho) < 10*A.Basis2().total_dimension() * ArnoldiTol);
#endif


#else


   DEBUG_CHECK(norm_frob(RightEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(RightEigen - triple_prod(herm(U), D*D, U)));

   // RightEigen = triple_prod(U, D*D, herm*U)

   A = delta_shift(U, QShift);
   AInv = adjoint(U);

   //   PsiL = inject_left_old_interface(A, PsiL);
   //   PsiL.set_back(prod(PsiL.get_back(), A * AInv));

   AInv = right_orthogonalize(PsiL, AInv);
   PsiL.set_front(prod(A*AInv, PsiL.get_front()));

   //   PsiL.set_back(prod(PsiL.get_back(), adjoint(U)));
   //   PsiL.set_front(prod(delta_shift(U, QShift), PsiL.get_front()));



   // orthonormalize each component of PsiL
   MatrixOperator I = MatrixOperator::make_identity(PsiL.Basis1());
   I = left_orthogonalize(I, PsiL);
   // now I should be unitary (or even an identity operator).  Normalize it and fold it back into the last component
   I *= 1.0 / std::sqrt(EtaR.real());
   PsiL.set_back(prod(PsiL.get_back(), I));
#endif

   return InfiniteWavefunctionLeft::ConstructFromOrthogonal(PsiL, D, QShift, LogAmplitude, Verbose-1);
}

MatrixOperator
InfiniteWavefunctionLeft::Initialize(LinearWavefunction const& Psi_, MatrixOperator const& Lambda, int Verbose)
{
   if (Verbose > 0)
   {
      std::cout << "Constructing canonical wavefunction..." << std::endl;
      std::cout << "Constructing right ortho matrices..." << std::endl;
   }

   LinearWavefunction Psi = Psi_;
   MatrixOperator M = right_orthogonalize(Psi, Lambda, Verbose-1);
   // normalize
   M *= 1.0 / norm_frob(M);

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   // we can't initialize lambda_l yet, as we don't have it in the correct (diagonal) basis.
   // We will only get this at the end, when we obtain lambda_r.  So just set it to a dummy
   this->push_back_lambda(D);

   if (Verbose > 0)
      std::cout << "Constructing left ortho matrices..." << std::endl;

   int n = 0;
   for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I, ++n)
   {
      if (Verbose > 1)
         std::cout << "orthogonalizing site " << n << std::endl;
      StateComponent A = prod(M, *I);
      M = ExpandBasis2(A);
      SingularValueDecomposition(M, U, D, Vh);
      this->push_back(prod(A, U));
      this->push_back_lambda(D);
      M = D*Vh;
   }

   Vh = delta_shift(Vh, QShift);
   this->set(0, prod(Vh, this->operator[](0)));
   this->setBasis1(Vh.Basis1());
   this->set_lambda(0, delta_shift(D, QShift));
   this->setBasis2(D.Basis1());

   if (Verbose > 0)
      std::cout << "Finished constructing canonical wavefunction." << std::endl;

   return Vh;
}

void read_version(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi, int Version)
{
   if (Version == 1)
   {
      MatrixOperator C_old;
      QuantumNumbers::QuantumNumber QShift;
      LinearWavefunction PsiLinear;
      MatrixOperator C_right;
      AttributeList Attr;
      AttributeList AttrDummy;

      in >> C_old;
      in >> QShift;
      in >> PsiLinear;
      in >> AttrDummy;
      in >> C_right;
      in >> Attr;

      Psi = InfiniteWavefunctionLeft::ConstructFromOrthogonal(PsiLinear, C_right, QShift);
   }
   else
   {
      int BaseVersion = Psi.CanonicalWavefunctionBase::ReadStream(in);
      in >> Psi.QShift;
      if (BaseVersion < 3)
         Psi.set_lambda(0, delta_shift(Psi.lambda_r(), Psi.qshift()));

      if (Version == 3)
      {
         // new in version 3, but replaced by LogAmplitude in version 4
         double Amplitude;
         in >> Amplitude;
         Psi.LogAmplitude = std::log(Amplitude);
      }
      else if (Version >= 4)
      {
         in >> Psi.LogAmplitude;
      }
      else
      {
         // prior to version 3, the amplitude is assumed to be always 1
         Psi.LogAmplitude = 0.0;
      }

      if (Version > 4)
      {
         PANIC("This program is too old to read this wavefunction, expected Version <= 4")(Version);
      }
   }

   Psi.debug_check_structure();
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi)
{
   // We formerly didn't have a version number for InfiniteWavefunction, so we have a hack,
   // and read the metadata version of the filesystem.  This is 1 if the file is old and
   // doesn't contain an InfiniteWavefunction version number.
   int Version;

   PHeapFileSystem::ipheapstream* S = dynamic_cast<PHeapFileSystem::ipheapstream*>(&in);
   if (S && S->version() == 1)
   {
      // old file, need to set the version number manually
      Version = 1;
   }
   else
   {
      // new file, read the version number
      Version = in.read<int>();
   }

   PStream::VersionSentry Sentry(in, InfiniteWavefunctionLeft::VersionT, Version);
   read_version(in, Psi, Version);

   return in;
}

PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunctionLeft const& Psi)
{
   out << InfiniteWavefunctionLeft::VersionT.default_version();

   Psi.CanonicalWavefunctionBase::WriteStream(out);

   out << Psi.QShift;
   out << Psi.LogAmplitude;

   return out;
}

std::pair<LinearWavefunction, RealDiagonalOperator>
get_left_canonical(InfiniteWavefunctionLeft const& Psi)
{
   return std::make_pair(LinearWavefunction(Psi.Basis1().GetSymmetryList(),
                                            Psi.base_begin(), Psi.base_end()), Psi.lambda_r());
}

std::tuple<MatrixOperator, RealDiagonalOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionLeft const& Psi)
{
   LinearWavefunction Result;
   RealDiagonalOperator D = Psi.lambda_r();
   MatrixOperator U = MatrixOperator::make_identity(D.Basis1());
   MatrixOperator Vh;
   InfiniteWavefunctionLeft::const_mps_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;

      StateComponent A = prod(*I, U*D);
      MatrixOperator M = ExpandBasis1(A);
      SingularValueDecomposition(M, U, D, Vh);
      Result.push_front(prod(Vh, A));
   }

   return std::make_tuple(U, D, Result);
}

void InfiniteWavefunctionLeft::scale(std::complex<double> x)
{
   double Amplitude = std::abs(x);
   x = x / Amplitude;  // x is now the phase
   LogAmplitude += std::log(Amplitude);
   // Because the components have value semantics we can't directly load() the pvalue_ptr
   // and then mutate it, because pvalue_handle's are immutable - we would end up referring
   // do a different object.
   pvalue_ptr<StateComponent> s = this->base_begin_()->load();
   *s.mutate() *= x;
   *this->base_begin_() = s;
}

InfiniteWavefunctionLeft& operator*=(InfiniteWavefunctionLeft& psi, std::complex<double> x)
{
   psi.scale(x);
   return psi;
}

void
InfiniteWavefunctionLeft::rotate_left(int Count)
{
   // Rotation is fairly straightforward, we just rotate the vectors around
   if (Count < 0)
   {
      this->rotate_right(-Count);
      return;
   }

   Count = Count % this->size();
   if (Count == 0)
      return;

   // the first Count elements are going to get shifted to the right hand side, so we need to
   // delta_shift them
   for (mps_iterator I = this->begin_(); I != this->begin_()+Count; ++I)
   {
      I->delta_shift(adjoint(this->qshift()));
   }
   // now do the actual rotation
   std::rotate(this->base_begin_(), this->base_begin_()+Count, this->base_end_());

   // for the Lambda matrices, start by removing the double-counted boundary lambda
   this->pop_back_lambda();
   // and delta-shift
   for (lambda_iterator I = this->lambda_begin_(); I != this->lambda_begin_()+Count; ++I)
   {
      I->delta_shift(adjoint(this->qshift()));
   }
   // and rotate
   std::rotate(this->lambda_base_begin_(), this->lambda_base_begin_()+Count, this->lambda_base_end_());
   // and put back the boundary lambda
   this->push_back_lambda(delta_shift(this->lambda_l(), adjoint(this->qshift())));

   // set the left and right basis
   this->setBasis1(lambda_l().Basis1());
   this->setBasis2(lambda_r().Basis2());

   this->debug_check_structure();
}

void
InfiniteWavefunctionLeft::rotate_right(int Count)
{
   if (Count < 0)
   {
      this->rotate_left(-Count);
      return;
   }

   Count = Count % this->size();
   if (Count == 0)
      return;

   this->rotate_left(this->size() - Count);
}

void
InfiniteWavefunctionLeft::check_structure() const
{
   this->CanonicalWavefunctionBase::check_structure();

   CHECK_EQUAL(this->Basis1(), delta_shift(this->Basis2(), this->qshift()));
}

void inplace_reflect(InfiniteWavefunctionLeft& Psi)
{
   // Although we construct the reflection into a new wavefunction,
   // we do gain an advantage for the inplace operation because we remove
   // elements from Psi as we go, so memory and disk use is ~ constant
   InfiniteWavefunctionLeft Result;
   Result.QShift = Psi.qshift();

   int Size = Psi.size();

   RealDiagonalOperator D = Psi.lambda_r();
   MatrixOperator U = MatrixOperator::make_identity(D.Basis1());

   // left-most lambda matrix, this is a place-holder
   Result.push_back_lambda(flip_conj(D));
   RealDiagonalOperator DSave = D;

   InfiniteWavefunctionLeft::const_base_mps_iterator I = Psi.base_end();
   while (I != Psi.base_begin())
   {
      --I;

      StateComponent A = prod(*I->lock(), U*D);
      MatrixOperator M = ExpandBasis1(A);

      MatrixOperator Vh;
      SingularValueDecomposition(M, U, D, Vh);
      A = prod(Vh, A);

      Result.push_back(reflect(A));
      Result.push_back_lambda(flip_conj(D));

      Psi.pop_back();
      Psi.pop_back_lambda();
   }

   // We have the final U matrix to deal with.  U.Basis1() is our original basis,
   // U.Basis2() has some gauge transform.  We want to write the final wavefunction in the
   // (flip conjugate) of the original basis.  If the wavefunction is perfectly orthogonal then
   // we should satisfy DSave = U*D*herm(U) here.  We want to use the DSave basis, not the D basis.
   // So we need to form
   // U*D (from the final SVD) = U*D*herm(U)*U = DSave * U
   // and merge the U into the final A-matrix.

#if 1
   //   Result.set_lambda(0, delta_shift(flip_conj(DSave), Psi.qshift()));
   //   Result.set_lambda(Size, flip_conj(DSave));
   Result.set_lambda(0, flip_conj(DSave));
   Result.set_lambda(Size, delta_shift(flip_conj(DSave), adjoint(Psi.qshift())));

   Result.set(Size-1, prod(Result[Size-1], herm(flip_conj(U))));

   Result.setBasis1(Result[0].Basis1());
   Result.setBasis2(Result[Size-1].Basis2());
#else
   // old code that used the D basis (and hence introduces a gauge transformation)
   Result.set_lambda(0, delta_shift(flip_conj(D), Psi.qshift()));
   Result.set(0, prod(herm(delta_shift(flip_conj(U), Psi.qshift())), Result[0]));

   Result.setBasis1(Result[0].Basis1());
   Result.setBasis2(adjoint(D.Basis2()));
#endif

   Psi = Result;

   Psi.check_structure();
}

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionLeft& Psi)
{
   for (InfiniteWavefunctionLeft::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = conj(*I);
   }
}

void inplace_qshift(InfiniteWavefunctionLeft& Psi, QuantumNumbers::QuantumNumber const& Shift)
{
   Psi.setBasis1(delta_shift(Psi.Basis1(), Shift));
   Psi.setBasis2(delta_shift(Psi.Basis2(), Shift));

   for (InfiniteWavefunctionLeft::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = delta_shift(*I, Shift);
   }

   for (InfiniteWavefunctionLeft::lambda_iterator I = Psi.lambda_begin_(); I != Psi.lambda_end_(); ++I)
   {
      *I = delta_shift(*I, Shift);
   }

   Psi.check_structure();
}

InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count)
{
   CHECK(Count >= 0);

   if (Count == 0 || Psi.empty())
      return InfiniteWavefunctionLeft();

   if (Count == 1)
      return Psi;

   // Claculate the final QShift, we need this for later
   QuantumNumber FinalQShift(Psi.GetSymmetryList());
   for (int i = 0; i < Count; ++i)
      FinalQShift = delta_shift(FinalQShift, Psi.qshift());

   // For the lambda matrices, we don't want to add the boundaries twice
   InfiniteWavefunctionLeft::const_lambda_iterator LambdaE = Psi.lambda_end();
   --LambdaE;

   // we need to handle the delta shift
   InfiniteWavefunctionLeft Result;
   bool First = true;
   for ( ; Count > 1; --Count)
   {
      QuantumNumber q(Psi.GetSymmetryList());
      for (int i = 1; i < Count; ++i)
         q = delta_shift(q, Psi.qshift());

      for (InfiniteWavefunctionLeft::const_mps_iterator I = Psi.begin(); I != Psi.end(); ++I)
      {
         // The very first time through the loop, set the Basis1.
         // We are guaranteed to go through this loop at least once, because
         // we return early if Count < 2
         if (First)
         {
            StateComponent A = delta_shift(*I, q);
            Result.setBasis1(A.Basis1());
            Result.push_back(A);
            Result.QShift = q;
            First = false;
         }
         else
         {
            Result.push_back(delta_shift(*I, q));
         }
      }

      for (InfiniteWavefunctionLeft::const_lambda_iterator I = Psi.lambda_begin(); I != LambdaE; ++I)
      {
         Result.push_back_lambda(delta_shift(*I, q));
      }

   }

   // The last repetition can make a simple copy
   for (InfiniteWavefunctionLeft::const_base_mps_iterator I = Psi.base_begin(); I != Psi.base_end(); ++I)
   {
      Result.push_back(*I);
   }

   // Here, we do want to go right to the end
   for (InfiniteWavefunctionLeft::const_base_lambda_iterator I = Psi.lambda_base_begin();
        I != Psi.lambda_base_end(); ++I)
   {
      Result.push_back_lambda(*I);
   }


   Result.setBasis2(Psi.Basis2());
   Result.QShift = FinalQShift;
   Result.LogAmplitude = Count * Psi.log_amplitude();

   return Result;
}

std::tuple<std::complex<double>, int, StateComponent>
overlap(InfiniteWavefunctionLeft const& x, ProductMPO const& StringOp,
        InfiniteWavefunctionLeft const& y,
        QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, int Verbose)
{
   int Length = statistics::lcm(x.size(), y.size(), StringOp.size());

   LinearWavefunction xPsi = get_left_canonical(x).first;
   LinearWavefunction yPsi = get_left_canonical(y).first;

   ProductMPO Str = StringOp * ProductMPO::make_identity(StringOp.LocalBasis2List(), Sector);

   StateComponent Init = MakeRandomStateComponent(Str.Basis1(), x.Basis1(), y.Basis1());

   int Iterations = Iter;
   int TotalIterations = 0;
   double MyTol = Tol;
   if (Verbose > 1)
   {
      std::cerr << "Starting Arnoldi, Tol=" << MyTol << ", Iterations=" << Iter << '\n';
   }
   std::complex<double> Eta = LinearSolvers::Arnoldi(Init,
                                                     LeftMultiplyOperator(xPsi, x.qshift(), Str,
                                                                          yPsi, y.qshift(), Length),
                                                     Iterations,
                                                     MyTol,
                                                     LinearSolvers::LargestMagnitude, false, Verbose);
   TotalIterations += Iterations;
   DEBUG_TRACE(Eta)(Iterations);

   while (MyTol < 0)
   {
      if (Verbose > 0)
         std::cerr << "Restarting Arnoldi, eta=" << Eta << ", Tol=" << -MyTol << '\n';
      Iterations = Iter;
      MyTol = Tol;
      Eta = LinearSolvers::Arnoldi(Init, LeftMultiplyOperator(xPsi, x.qshift(), Str,
                                                              yPsi, y.qshift(), Length),
                                   Iterations, MyTol, LinearSolvers::LargestMagnitude, false, Verbose);
      TotalIterations += Iterations;
      DEBUG_TRACE(Eta)(Iterations);
   }
   if (Verbose > 0)
      std::cerr << "Converged.  TotalIterations=" << TotalIterations
                << ", Tol=" << MyTol << '\n';

   double LogAmplitude = x.log_amplitude() * (Length / x.size()) + y.log_amplitude() * (Length / y.size());
   return std::make_tuple(Eta * std::exp(LogAmplitude), Length, Init);
}

#if 0
std::complex<double> overlap(InfiniteWavefunctionLeft const& x, BasicFiniteMPO const& StringOp,
                             InfiniteWavefunctionLeft const& y,
                             QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, int Verbose)
{
   CHECK_EQUAL(x.qshift(), y.qshift())("The wavefunctions must have the same quantum number per unit cell");

   LinearWavefunction xPsi = get_left_canonical(x).first;
   LinearWavefunction yPsi = get_left_canonical(y).first;

   MatrixOperator Init = MakeRandomMatrixOperator(x.Basis1(), y.Basis1(), Sector);

   int Iterations = Iter;
   int TotalIterations = 0;
   double MyTol = Tol;
   std::complex<double> Eta = LinearSolvers::Arnoldi(Init,
                LeftMultiplyString(xPsi, StringOp * MakeIdentityFrom(StringOp, Sector),
                                   yPsi, x.qshift()),
                                                     Iterations,
                                                     MyTol,
                                                     LinearSolvers::LargestMagnitude, false, Verbose);
   TotalIterations += Iterations;
   DEBUG_TRACE(Eta)(Iterations);

   while (MyTol < 0)
   {
      if (Verbose > 0)
         std::cerr << "Restarting Arnoldi, eta=" << Eta << ", Tol=" << -MyTol << '\n';
      Iterations = Iter;
      MyTol = Tol;
      Eta = LinearSolvers::Arnoldi(Init, LeftMultiplyString(xPsi, StringOp * MakeIdentityFrom(StringOp, Sector),
                                                            yPsi, x.qshift()),
                                   Iterations, MyTol, LinearSolvers::LargestMagnitude, false, Verbose);
      TotalIterations += Iterations;
      DEBUG_TRACE(Eta)(Iterations);
   }
   if (Verbose > 0)
      std::cerr << "Converged.  TotalIterations=" << TotalIterations
                << ", Tol=" << MyTol << '\n';

   return Eta;
}
#endif

std::pair<std::complex<double>, StateComponent>
overlap(InfiniteWavefunctionLeft const& x,  InfiniteWavefunctionLeft const& y,
        QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, int Verbose)
{
   CHECK_EQUAL(x.size(), y.size());
   std::tuple<std::complex<double>, int, StateComponent> Result =
      overlap(x, ProductMPO::make_identity(ExtractLocalBasis(y)), y, Sector, Iter, Tol, Verbose);
   return std::make_pair(std::get<0>(Result), std::get<2>(Result));
}

InfiniteWavefunctionLeft
reflect(InfiniteWavefunctionRight const& Psi)
{
   PANIC("not implemented");
}

void
InfiniteWavefunctionLeft::SetDefaultAttributes(AttributeList& A) const
{
   A["WavefunctionType"] = "Infinite";
   A["UnitCellSize"] = this->size();
   A["QShift"] = this->qshift();
}

// inject_left for a BasicFiniteMPO.  This can have support on multiple wavefunction unit cells
MatrixOperator
inject_left(MatrixOperator const& m,
            InfiniteWavefunctionLeft const& Psi1,
            BasicFiniteMPO const& Op,
            InfiniteWavefunctionLeft const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.qshift(), Psi2.qshift());
   DEBUG_CHECK_EQUAL(m.Basis1(), Psi1.Basis1());
   DEBUG_CHECK_EQUAL(m.Basis2(), Psi2.Basis1());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   MatrixOperator Result = m;
   StateComponent E(Op.Basis1(), m.Basis1(), m.Basis2());
   E[0] = m;
   E.debug_check_structure();
   InfiniteWavefunctionLeft::const_mps_iterator I1 = Psi1.begin();
   InfiniteWavefunctionLeft::const_mps_iterator I2 = Psi2.begin();
   BasicFiniteMPO::const_iterator OpIter = Op.begin();
   while (OpIter != Op.end())
   {
      if (I1 == Psi1.end())
      {
         I1 = Psi1.begin();
         I2 = Psi2.begin();
         E = delta_shift(E, Psi1.qshift());
      }
      E = contract_from_left(*OpIter, herm(*I1), E, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return delta_shift(E[0], Psi1.qshift());
}

std::complex<double>
expectation(InfiniteWavefunctionLeft const& Psi, BasicFiniteMPO const& Op)
{
#if 0
   MatrixOperator X = MatrixOperator::make_identity(Psi.Basis1());
   X = inject_left(X, Psi, Op, Psi);

   MatrixOperator Rho = Psi.lambda_r();
   Rho = Rho*Rho;

   return inner_prod(delta_shift(Rho, Psi.qshift()), X);
#else
   // find the first non-trivial site of the MPO
   BasicFiniteMPO::const_iterator i = Op.begin();
   int n = 0;
   OperatorClassification c;
   c.Identity_ = true;
   if (i != Op.end())
      c = classify(*i);
   while (i != Op.end() && c.is_identity())
   {
      ++i; ++n;
      if (i != Op.end())
         c = classify(*i);
   }
   // find the last non-trivial site of the MPO
   BasicFiniteMPO::const_iterator j = Op.end();
   int m = Op.size();
   c = OperatorClassification();
   c.Identity_ = true;
   while (j != i && c.is_identity())
   {
      --j; --m;
      c = classify(*j);
   }
   if (i == j && c.is_identity())
      return 1.0;
   // we could optimize for the case i==j and c.is_identity()
   ++j; ++m;
   // get an iterator into the wavefunction
   InfiniteWavefunctionLeft::const_mps_iterator I = Psi.begin();
   for (int a = 0; a < n; ++a)
   {
      ++I;
      if (I == Psi.end())
         I == Psi.begin();
   }
   //TRACE(m-n);
   MatrixOperator X = MatrixOperator::make_identity(I->Basis1());
   CHECK_EQUAL(i->Basis1().size(), 1);
   CHECK(is_scalar(i->Basis1()[0]));
   StateComponent E(i->Basis1(), X.Basis1(), X.Basis2());
   E[0] = X;
   while (n < m)
   {
      if (I == Psi.end())
      {
         I = Psi.begin();
         E = delta_shift(E, Psi.qshift());
      }
      E = contract_from_left(*i, herm(*I), E, *I);
      ++I;
      ++i;
      ++n;
   }
   m = m % Psi.size();
   if (m == 0)
      m = Psi.size();    // if we're at the edge of the unit cell, use lambda_r
   MatrixOperator Rho = Psi.lambda(m);
   Rho = Rho*Rho;

   return inner_prod(Rho, E[0]);
#endif
}

#if 0
InfiniteWavefunctionLeft
fine_grain(InfiniteWavefunctionLeft const& x, int N, std::vector<BasisList> const& FineGrainBasis,
	   SimpleOperator const& U, StatesInfo const& SInfo)
{
   // better approach: construct a right-canonical wavefunction and then do a left-right sweep

   CHECK_EQUAL(x.size()*N, FineGrainBasis.size());
   LinearWavefunction Result(x.GetSymmetryList());
   for (int n = 0; n < x.size(); ++n)
   {
      RealDiagonalOperator Lambda = x.lambda(n+1);
      StateComponent A = x[n];

      std::stack<BasisList> B1;
      B1.push(FineGrainBasis[n*N]);
      for (int i = 0; i < N-1; ++i)
      {
	 B1.push(make_product_basis(B1.top(), FineGrainBasis[n*N+i]));
      }
      std::stack<StateComponent> R;
      while (!B1.empty())
      {
	 AMatSVD SL(A*Lambda,
		    Tensor::ProductBasis<BasisList, BasisList>(B1.top, FineGrainBasis[n*N+B1.size()]));
	 AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
								     SInfo, Info);
	 StateComponent B;
	 SL.ConstructMatrices(SL.begin(), Cutoff, A, Lambda, B);
	 // normalize
	 Lambda *= 1.0 / norm_frob(Lambda);
	 R.push(B);
	 B1.pop();
      }

)
#endif
