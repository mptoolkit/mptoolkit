// -*- C++ -*- $Id$

#include "infinitewavefunctionleft.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "common/environment.h"
#include "common/statistics.h"

#include "mps/packunpack.h"
#include <fstream>

#include "wavefunction/operator_actions.h"
#include "pheap/pheapstream.h"

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

PStream::VersionTag
InfiniteWavefunctionLeft::VersionT(2);

extern double const ArnoldiTol = getenv_or_default("MP_ARNOLDI_TOL", 1E-15);

extern double const InverseTol = getenv_or_default("MP_INVERSE_TOL", 1E-7);

// the tol used in the orthogonalization can apparently be a bit smaller
extern double const OrthoTol = getenv_or_default("MP_ORTHO_TOL", 1E-9);

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
   // normalize
   M *= 1.0 / norm_frob(M);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   // we can't initialize lambda_l yet, as we don't have it in the correct (diagonal) basis.
   // We will only get this at the end, when we obtain lambda_r.  So just set it to a dummy
   this->push_back_lambda(D);
   for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
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
}

InfiniteWavefunctionLeft::InfiniteWavefunctionLeft(LinearWavefunction const& Psi, 
						   QuantumNumbers::QuantumNumber const& QShift_)
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
                                                      Iterations, Tol, 
						      LinearSolvers::LargestAlgebraicReal, false);
   while (Tol < 0)
   {
      std::cerr << "LeftEigen: Arnoldi not converged, restarting.  EValue=" 
                << EtaL << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, QShift), Iterations, 
				    Tol, LinearSolvers::LargestAlgebraicReal, false);
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
      (norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)));

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
   A = LeftMultiply(PsiL, QShift)(I);
   CHECK(norm_frob(A-EtaL*I) < 10*A.Basis1().total_dimension() * ArnoldiTol)(norm_frob(A-EtaL*I))(A)(I)(D);
#endif 

   // same for the right eigenvector, which will be the density matrix

   // initialize the guess eigenvector
   MatrixOperator RightEigen = Guess;

   // get the eigenmatrix
   Iterations = 20; Tol = ArnoldiTol;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, QShift), 
                                                      Iterations, Tol, 
						      LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Tol < 0)
   {
      std::cerr << "RightEigen: Arnoldi not converged, restarting.  EValue=" 
                << EtaR << ", Tol=" << Tol << "\n";
      Iterations = 20; Tol = ArnoldiTol;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, QShift), 
				    Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }
   DEBUG_TRACE(EtaR);

   if (trace(RightEigen).real() < 0)
      RightEigen = -RightEigen;

   DEBUG_CHECK(norm_frob(RightEigen - adjoint(RightEigen)) < 1e-12)
      (norm_frob(RightEigen - adjoint(RightEigen)));

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

   this->Initialize(PsiL, D);
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

      in >> C_old;
      in >> QShift;
      in >> PsiLinear;
      in >> C_right;
      in >> Attr;
   
      Psi = InfiniteWavefunctionLeft(PsiLinear, C_right, QShift);
   }
   else if (Version == 2)
   {
      int BaseVersion = Psi.CanonicalWavefunctionBase::ReadStream(in);
      in >> Psi.QShift;
      if (BaseVersion < 3)
	 Psi.set_lambda(0, delta_shift(Psi.lambda_r(), Psi.qshift()));
   }
   else
   {
      PANIC("This program is too old to read this wavefunction, expected Version <= 2")(Version);
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

   return out;
}

std::pair<LinearWavefunction, RealDiagonalOperator>
get_left_canonical(InfiniteWavefunctionLeft const& Psi)
{
   return std::make_pair(LinearWavefunction(Psi.base_begin(), Psi.base_end()), Psi.lambda(Psi.size()));
}

boost::tuple<MatrixOperator, RealDiagonalOperator, LinearWavefunction>
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

   return boost::make_tuple(U, D, Result);
}

void
InfiniteWavefunctionLeft::rotate_left(int Count)
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
   if (Count == 0)
      return;
   else if (Count < 0)
      this->rotate_left(-Count);

   Count = Count % this->size();

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

   RealDiagonalOperator D = Psi.lambda_r();
   MatrixOperator U = MatrixOperator::make_identity(D.Basis1());

   // left-most lambda matrix, this is a place-holder
   Result.push_back_lambda(flip_conj(D));
   MatrixOperator DSave = D;

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

   // In principle DSave is U*D*herm(U)

   Result.set_lambda(0, delta_shift(flip_conj(D), Psi.qshift()));
   Result.set(0, prod(herm(delta_shift(flip_conj(U), Psi.qshift())), Result[0]));
   
   Result.setBasis1(Result[0].Basis1());
   Result.setBasis2(adjoint(D.Basis2()));

   Psi = Result;

   Psi.debug_check_structure();
}

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionLeft& Psi)
{
   for (InfiniteWavefunctionLeft::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = conj(*I);
   }
}

InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count)
{
   CHECK(Count >= 0);

   if (Count == 0 || Psi.empty())
      return InfiniteWavefunctionLeft();

   if (Count == 1)
      return Psi;

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

   return Result;
}

std::pair<std::complex<double>, int>
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

   return std::make_pair(Eta, Length);
}

std::complex<double> overlap(InfiniteWavefunctionLeft const& x, FiniteMPO const& StringOp,
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
      if (Verbose)
         std::cerr << "Restarting Arnoldi, eta=" << Eta << ", Tol=" << -MyTol << '\n';
      Iterations = Iter;
      MyTol = Tol;
      Eta = LinearSolvers::Arnoldi(Init, LeftMultiplyString(xPsi, StringOp * MakeIdentityFrom(StringOp, Sector), 
							    yPsi, x.qshift()), 
				   Iterations, MyTol, LinearSolvers::LargestMagnitude, false, Verbose);
      TotalIterations += Iterations;
      DEBUG_TRACE(Eta)(Iterations);
   }
   if (Verbose)
      std::cerr << "Converged.  TotalIterations=" << TotalIterations
                << ", Tol=" << MyTol << '\n';

   return Eta;
}

std::complex<double> overlap(InfiniteWavefunctionLeft const& x,  InfiniteWavefunctionLeft const& y,
                             QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, int Verbose)
{
   return overlap(x, ProductMPO::make_identity(ExtractLocalBasis(y)), y, Sector, Iter, Tol, Verbose).first;
}

InfiniteWavefunctionLeft
reflect(InfiniteWavefunctionRight const& Psi)
{
   PANIC("not implemented");
}

