// -*- C++ -*- $Id$

#include "infinitewavefunction.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
//#include "wavefunc-utils.h"
#include "common/environment.h"

#include "mps/packunpack.h"
#include <fstream>

#include "mps/spectrum_arpack.h"
#include "mps/operator_actions.h"

#include "pheap/pheapstream.h"


double const ArnoldiTol = getenv_or_default("MP_ARNOLDI_TOL", 1E-15);

double const InverseTol = getenv_or_default("MP_INVERSE_TOL", 1E-7);

// the tol used in the orthogonalization can apparently be a bit smaller
double const OrthoTol = getenv_or_default("MP_ORTHO_TOL", 1E-9);

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



PStream::VersionTag
InfiniteWavefunctionLeft::VersionT(2);

InfiniteWavefunctionLeft::InfiniteWavefunctionLeft(LinearWavefunction const& Psi, MatrixOperator const& Lambda,
						   QuantumNumbers::QuantumNumber const& QShift_)
   : QShift(QShift_)
{
   this->Initialize(Psi, Lambda);
}

void
InfiniteWavefunctionLeft::Initialize(LinearWavefunction const& Psi_, MatrixOperator const& Lambda)
{
   LinearWavefunction Psi = Psi_;
   MatrixOperator M = right_orthogonalize(Psi, Lambda);
   this->setBasis1(M.Basis1());
   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
      StateComponent A = prod(M, *I);
      M = ExpandBasis2(A);
      SingularValueDecomposition(M, U, D, Vh);
      this->push_back(prod(A, U));
      this->push_back_lambda(D);
      M = D*Vh;
   }

   this->set(0, prod(delta_shift(Vh, QShift), this->operator[](0)));

   this->setBasis2(M.Basis1());
}

InfiniteWavefunctionLeft::InfiniteWavefunctionLeft(LinearWavefunction const& Psi, QuantumNumbers::QuantumNumber const& QShift_)
   : QShift(QShift_)
{
   LinearWavefunction PsiL = Psi;

   MatrixOperator Guess = MatrixOperator::make_identity(PsiL.Basis2());
   
   // initialize LeftEigen to a guess eigenvector.  Since L satisfies the left orthogonality
   // constraint (except for the final matrix), we can do one iteration beyond the identity
   // and intialize it to herm(Xu) * Xu
   MatrixOperator LeftEigen = Guess;

      // get the eigenmatrix.  Do some dodgy explict restarts.
   int Iterations = 20;
   double Tol = ArnoldiTol;
   LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   while (Iterations == 20)
   {
      std::cerr << "LeftEigen: Arnoldi not converged, restarting.  EValue=" 
                << EtaL << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   CHECK(EtaL.real() > 0)("Eigenvalue must be positive");

   DEBUG_TRACE(EtaL);

   if (trace(LeftEigen).real() < 0)
      LeftEigen = -LeftEigen;

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

   // LeftEigen = triple_prod(U, D*D, herm(U))
   DEBUG_CHECK(norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)));

   MatrixOperator A = delta_shift(D * U, QShift);
   MatrixOperator AInv = adjoint(U) * DInv;

   A = left_orthogonalize(A, PsiL);
   PsiL.set_back(prod(PsiL.get_back(), A*AInv));

   // same for the right eigenvector

   // initialize the guess eigenvector
   MatrixOperator RightEigen = Guess;

   // get the eigenmatrix
   Iterations = 20; Tol = ArnoldiTol;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Iterations == 20)
   {
      std::cerr << "RightEigen: Arnoldi not converged, restarting.  EValue=" 
                << EtaR << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }
   DEBUG_TRACE(EtaR);

   if (trace(RightEigen).real() < 0)
      RightEigen = -RightEigen;

   DEBUG_CHECK(norm_frob(RightEigen - adjoint(RightEigen)) < 1e-12)
      (norm_frob(RightEigen - adjoint(RightEigen)));

   D = RightEigen;
   U = DiagonalizeHermitian(D);
   D = SqrtDiagonal(D, OrthoTol);

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

   // normalize
   D *= 1.0 / norm_frob(D);

   this->Initialize(PsiL, D);
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

   if (Version == 1)
   {
      MatrixOperator C_old;
      QuantumNumbers::QuantumNumber QShift;
      LinearWavefunction PsiLinear;
      MatrixOperator C_right;
      AttributeList Attr;

      in >> C_old;
      in >> QShift;
      in >> PsiLinear;
      in >> C_right;
      in >> Attr;
   
      Psi = InfiniteWavefunctionLeft(PsiLinear, C_right, QShift);
   }
   else if (Version == 2)
   {
      Psi.CanonicalWavefunctionBase::ReadStream(in);
      in >> Psi.QShift;
   }
   else
   {
      PANIC("This program is too old to read this wavefunction");
   }

   return in;
}

PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunctionLeft const& Psi)
{
   out << InfiniteWavefunctionLeft::VersionT.default_version();

   Psi.CanonicalWavefunctionBase::WriteStream(out);

   out << Psi.QShift;

   return out;
}

std::pair<LinearWavefunction, RealDiagonalOperator>
get_left_canonical(InfiniteWavefunctionLeft const& Psi)
{
   return std::make_pair(LinearWavefunction(Psi.begin_raw(), Psi.end_raw()), Psi.lambda(Psi.size()-1));
}

#if 0
std::pair<MatrixOperator, LinearWavefunction>
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

   U = U*D;
   return std::make_pair(U, Result);
}
#endif

std::pair<MatrixOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionLeft const& Psi)
{
   LinearWavefunction Result;
   MatrixOperator M = Psi.lambda_r();
   InfiniteWavefunctionLeft::const_mps_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;

      StateComponent A = prod(*I, M);
      M = TruncateBasis1(A);
      Result.push_front(A);
   }

   return std::make_pair(M, Result);
}



