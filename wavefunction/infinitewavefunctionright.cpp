// -*- C++ -*- $Id$

#include "infinitewavefunctionright.h"
#include "infinitewavefunctionleft.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "common/environment.h"

#include "mps/packunpack.h"
#include <fstream>

#include "wavefunction/operator_actions.h"
#include "pheap/pheapstream.h"

// Streaming versions:
// Note: the base class CanonicalWavefunctionBase has a separate version number.
//
// Version 1:
//      CanonicalWavefunctionBase (base class)
//      QuantumNumber              QShift

extern double const ArnoldiTol;
extern double const InverseTol;
extern double const OrthoTol;

PStream::VersionTag
InfiniteWavefunctionRight::VersionT(1);

struct LeftMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiply(LinearWavefunction const& L_, QuantumNumber const& QShift_) 
      : L(L_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      for (LinearWavefunction::const_iterator I = L.begin(); I != L.end(); ++I)
      {
	 r = operator_prod(herm(*I), r, *I);
      }
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
};
   
struct RightMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiply(LinearWavefunction const& R_, QuantumNumber const& QShift_) 
      : R(R_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = x;
      LinearWavefunction::const_iterator I = R.end();
      while (I != R.begin())
      {
	 --I;
	 r = operator_prod(*I, r, herm(*I));
      }
      return delta_shift(r, adjoint(QShift));
   }

   LinearWavefunction const& R;
   QuantumNumber QShift;
};

InfiniteWavefunctionRight::InfiniteWavefunctionRight(MatrixOperator const& Lambda,
						     LinearWavefunction const& Psi, 
						     QuantumNumbers::QuantumNumber const& QShift_)
   : QShift(QShift_)
{
   this->Initialize(Lambda, Psi);
}

void
InfiniteWavefunctionRight::Initialize(MatrixOperator const& Lambda, 
				      LinearWavefunction const& Psi_)
{
   LinearWavefunction Psi = Psi_;
   MatrixOperator M = left_orthogonalize(Lambda, Psi);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   // we need to assemble the MPS in reverse order, so make a temporary container
   std::list<StateComponent> AMat;
   std::list<RealDiagonalOperator> Lam;

   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      StateComponent A = prod(*I, M);
      M = ExpandBasis1(A);
      SingularValueDecomposition(M, U, D, Vh);
      AMat.push_front(prod(Vh, A));
      Lam.push_front(D);
      M = U*D;
   }
   U = delta_shift(U, adjoint(QShift));
   AMat.back() = prod(AMat.back(), U);
   Lam.push_back(delta_shift(D, adjoint(QShift)));

   for (auto const& A : AMat)
   {
      this->push_back(A);
   }

   for (auto const& L : Lam)
   {
      this->push_back_lambda(L);
   }

   this->setBasis1(D.Basis1());
   this->setBasis2(U.Basis2());
}

InfiniteWavefunctionRight::InfiniteWavefunctionRight(LinearWavefunction const& Psi, 
						     QuantumNumbers::QuantumNumber const& QShift_)
   : QShift(QShift_)
{
   LinearWavefunction PsiR = Psi;

   MatrixOperator Guess = MatrixOperator::make_identity(PsiR.Basis2());
   
   MatrixOperator RightEigen = Guess;

   // get the eigenmatrix.  Do some dodgy explict restarts.
   int Iterations = 20;
   double Tol = ArnoldiTol;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen)); // make the eigenvector symmetric
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, QShift), 
                                                      Iterations, Tol, 
						      LinearSolvers::LargestAlgebraicReal, false);
   while (Tol < 0)
   {
      std::cerr << "RightEigen: Arnoldi not converged, restarting.  EValue=" 
                << EtaR << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen)); // make the eigenvector symmetric
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, QShift), 
				    Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   CHECK(EtaR.real() > 0)("Eigenvalue must be positive");

   DEBUG_TRACE(EtaR);

   if (trace(RightEigen).real() < 0)
      RightEigen = -RightEigen;

   DEBUG_CHECK(norm_frob(RightEigen - adjoint(RightEigen)) < 1e-12)
      (norm_frob(RightEigen - adjoint(RightEigen)));

   MatrixOperator D = RightEigen;
   MatrixOperator U = DiagonalizeHermitian(D);
   D = SqrtDiagonal(D, OrthoTol);
   MatrixOperator DInv = InvertDiagonal(D, OrthoTol);

   // RightEigen = triple_prod(U, D*D, herm(U))
   DEBUG_CHECK(norm_frob(RightEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(RightEigen - triple_prod(herm(U), D*D, U)));

   MatrixOperator R = adjoint(U)*D;
   MatrixOperator RInv = delta_shift(DInv * U, QShift);

   // Incorporate into the MPS: PsiR -> R^{-1} * PsiR * R
   // We don't need to actually right-orthogonalize everything, that is done in Initialize() anyway
   PsiR.set_back(prod(PsiR.get_back(), R));
   PsiR.set_front(prod(RInv, PsiR.get_front()));

   // Get the left eigenvector, which is the density matrix
   MatrixOperator LeftEigen = Guess;

   // get the eigenmatrix
   Iterations = 20; Tol = ArnoldiTol;
   LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen));
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiR, QShift), 
                                                      Iterations, Tol, 
						      LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(LeftEigen - adjoint(LeftEigen)));
   while (Tol < 0)
   {
      std::cerr << "LeftEigen: Arnoldi not converged, restarting.  EValue=" 
                << EtaL << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen));
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiR, QShift), 
				    Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   D = LeftEigen;
   U = DiagonalizeHermitian(D);
   D = SqrtDiagonal(D, OrthoTol);

   DEBUG_CHECK(norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)));

   // incorporate U into the MPS

   PsiR.set_front(prod(U, PsiR.get_front()));
   PsiR.set_back(prod(PsiR.get_back(), adjoint(U)));

   // And now we have the right-orthogonalized form and the left-most lambda matrix   
   this->Initialize(D, PsiR);
}

InfiniteWavefunctionRight::InfiniteWavefunctionRight(InfiniteWavefunctionLeft const& Psi)
   : InfiniteWavefunctionRight(LinearWavefunction(Psi.base_begin(), Psi.base_end()), Psi.qshift())
{
}

void read_version(PStream::ipstream& in, InfiniteWavefunctionRight& Psi, int Version)
{
   if (Version == 1)
   {
      int BaseVersion = Psi.CanonicalWavefunctionBase::ReadStream(in);
      in >> Psi.QShift;
      CHECK(BaseVersion >= 3);
   }
   else
   {
      PANIC("This program is too old to read this wavefunction, expected Version <= 2")(Version);
   }

   Psi.debug_check_structure();
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteWavefunctionRight& Psi)
{
   int Version = in.read<int>();
   PStream::VersionSentry Sentry(in, InfiniteWavefunctionRight::VersionT, Version);
   read_version(in, Psi, Version);
   return in;
}

PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunctionRight const& Psi)
{
   out << InfiniteWavefunctionRight::VersionT.default_version();

   Psi.CanonicalWavefunctionBase::WriteStream(out);

   out << Psi.QShift;

   return out;
}

std::pair<LinearWavefunction, RealDiagonalOperator>
get_right_canonical(InfiniteWavefunctionRight const& Psi)
{
   return std::make_pair(LinearWavefunction(Psi.base_begin(), Psi.base_end()), 
			 Psi.lambda(Psi.size()));
}

boost::tuple<MatrixOperator, RealDiagonalOperator, LinearWavefunction>
get_left_canonical(InfiniteWavefunctionRight const& Psi)
{
   LinearWavefunction Result;
   RealDiagonalOperator D = Psi.lambda_r();
   MatrixOperator U = MatrixOperator::make_identity(D.Basis1());
   MatrixOperator Vh;
   InfiniteWavefunctionRight::const_mps_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;

      StateComponent A = prod(*I, U*D);
      MatrixOperator M = ExpandBasis1(A);
      SingularValueDecomposition(M, U, D, Vh);
      Result.push_front(prod(Vh, A));
   }

   return boost::make_tuple(U, D, Result);
}

void
InfiniteWavefunctionRight::rotate_left(int Count)
{
   // Rotation is fairly straightforward, we just rotate the vectors around
   if (Count == 0)
      return;
   else if (Count < 0)
      this->rotate_right(-Count);

   Count = Count % this->size();

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
   std::rotate(this->lambda_base_begin_(), this->lambda_base_begin_()+Count, 
	       this->lambda_base_end_());
   // and put back the boundary lambda
   this->push_back_lambda(delta_shift(this->lambda_l(), adjoint(this->qshift())));

   // set the right and right basis
   this->setBasis1(lambda_l().Basis1());
   this->setBasis2(lambda_r().Basis2());

   this->debug_check_structure();
}

void
InfiniteWavefunctionRight::rotate_right(int Count)
{
   if (Count == 0)
      return;
   else if (Count < 0)
      this->rotate_left(-Count);

   Count = Count % this->size();

   this->rotate_left(this->size() - Count);
}

void
InfiniteWavefunctionRight::check_structure() const
{
   this->CanonicalWavefunctionBase::check_structure();

   CHECK_EQUAL(this->Basis1(), delta_shift(this->Basis2(), this->qshift()));
}

void inplace_conj(InfiniteWavefunctionRight& Psi)
{
   for (InfiniteWavefunctionRight::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = conj(*I);
   }
}

InfiniteWavefunctionRight
reflect(InfiniteWavefunctionLeft const& Psi)
{
   PANIC("not implemented");
}