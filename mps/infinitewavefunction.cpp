// -*- C++ -*- $Id$

#include "infinitewavefunction.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
//#include "wavefunc-utils.h"
#include "common/environment.h"

#include "mps/packunpack.h"
#include <fstream>

#include "mps/spectrum_arpack.h"

double InverseTol = getenv_or_default("MP_INVERSE_TOL", InverseTolDefault);

// the tol used in the orthogonalization can apparently be a bit smaller
double OrthoTol = getenv_or_default("MP_ORTHO_TOL", 1E-9);

InfiniteWavefunction rotate_left(InfiniteWavefunction const& Psi, int Count)
{
   // this function assumes the wavefunction is orthogonal
   LinearWavefunction Result(Psi.GetSymmetryList());
   LinearWavefunction::const_iterator Pivot = Psi.Psi.begin();
   for (int i = 0; i < Count; ++i)
      ++Pivot;

   LinearWavefunction::const_iterator I = Pivot;
   while (I != Psi.Psi.end())
      Result.push_back(*I++);

   I = Psi.Psi.begin();
   MatrixOperator C = delta_shift(Psi.C_right, Psi.QShift) * InvertDiagonal(Psi.C_old, InverseTol);
   while (I != Pivot)
   {
      StateComponent x = prod(C, *I);
      C = TruncateBasis2(x);
      Result.push_back(delta_shift(x, adjoint(Psi.QShift)));
      ++I;
   }

   Result.set_back(prod(Result.get_back(), (delta_shift(C, adjoint(Psi.QShift)))));
   
   InfiniteWavefunction Ret;
   Ret.C_old = MatrixOperator::make_identity(Result.Basis1());
   Ret.QShift = Psi.QShift;
   Ret.Psi = Result;
   Ret.C_right = MatrixOperator::make_identity(Result.Basis2());
   Ret.Attr = Psi.Attr;

   DEBUG_CHECK_EQUAL(Ret.C_old.Basis1(), Ret.C_old.Basis2());
   DEBUG_CHECK_EQUAL(Ret.C_old.Basis2(), Ret.Psi.Basis1());
   DEBUG_CHECK_EQUAL(Ret.Psi.Basis2(), Ret.C_right.Basis1());
   DEBUG_CHECK_EQUAL(Ret.Psi.Basis1(), DeltaShift(Ret.C_right.Basis2(), Ret.QShift));

   orthogonalize(Ret);

   return Ret;
}

PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunction const& psi)
{
   out << psi.C_old;
   out << psi.QShift;
   out << psi.Psi;
   out << psi.C_right;
   out << psi.Attr;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteWavefunction& psi)
{
   in >> psi.C_old;
   in >> psi.QShift;
   in >> psi.Psi;
   in >> psi.C_right;
   in >> psi.Attr;
   return in;
}

LinearWavefunction get_orthogonal_wavefunction(InfiniteWavefunction const& Psi)
{
   return Psi.Psi;
#if 0
   MatrixOperator x_unit =  Psi.C_right * delta_shift(InvertDiagonal(Psi.C_old, InverseTol), adjoint(Psi.QShift));
   LinearWavefunction xPsi = Psi.Psi;
   //xPsi.set_back(prod(xPsi.get_back(), x_unit));   
   xPsi.set_front(prod(delta_shift(x_unit, Psi.QShift), xPsi.get_front()));
   MatrixOperator I = MatrixOperator::make_identity(xPsi.Basis1());
   I = left_orthogonalize(I, xPsi);
   xPsi.set_back(prod(xPsi.get_back(), I));
   return xPsi;
#endif
}

LinearWavefunction get_right_orthogonal_wavefunction(InfiniteWavefunction const& Psi)
{
   MatrixOperator r = Psi.C_right;
   LinearWavefunction PsiR = Psi.Psi;
   r = right_orthogonalize(PsiR, r);
   PsiR.set_front(prod(InvertDiagonal(Psi.C_old, InverseTol) * r, PsiR.get_front()));
   return PsiR;
}

double orthogonality_fidelity(InfiniteWavefunction const& x)
{
   // The fidelity of two density matrices \rho_1 and \rho_2 is
   // Tr \sqrt{ \sqrt{\rho_1} \rho_2 \sqrt{\rho_1}}
   // which corresponds to Tr \sqrt{ C_right * herm(C_old) * C_old * herm(C_right)}
   // which we can do as the trace of the singular values of C_right * herm(C_old)
   // We don't assume that the state is normalized, although in practice it probably
   // always is.
   MatrixOperator U, D, Vh;
   SingularValueDecomposition(scalar_prod(x.C_right, 
                                          herm(delta_shift(x.C_old, adjoint(x.QShift)))),
                              U, D, Vh);
   return trace(D).real() / (norm_frob(x.C_right) * norm_frob(x.C_old));
}

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
   
struct GeneralizedLeftMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   GeneralizedLeftMultiply(LinearWavefunction const& L1_, 
                           LinearWavefunction const& L2_,
                           QuantumNumber const& QShift_) 
      : L1(L1_), L2(L2_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = x; //delta_shift(x, QShift);
      LinearWavefunction::const_iterator I1 = L1.begin();
      LinearWavefunction::const_iterator I2 = L2.begin();
      for ( ; I1 != L1.end(); ++I1, ++I2)
      {
	 r = operator_prod(herm(*I1), r, *I2);
      }
      return delta_shift(r, QShift);
   }

   LinearWavefunction const& L1;
   LinearWavefunction const& L2;
   QuantumNumber QShift;
};

#if 0
void orthogonalize(InfiniteWavefunction& x)
{
   // firstly, get the eigenmatrices of the left and right orthogonality superoperators

   // Incorporate C_old and C_right into the A-matrices
   LinearWavefunction PsiL = x.Psi;
   StateComponent Last = PsiL.get_back();  // right-most A-matrix
   MatrixOperator C_old_inverse = InvertDiagonal(x.C_old, InverseTol);
   MatrixOperator Xu = x.C_right * delta_shift(C_old_inverse, adjoint(x.QShift));
   Last = prod(Last, Xu);
   PsiL.set_back(Last);

   // initialize LeftEigen to a guess eigenvector.  Since L satisfies the left orthogonality
   // constraint (except for the final matrix), we can do one iteration beyond the identity
   // and intialize it to herm(Xu) * Xu
   MatrixOperator LeftEigen = scalar_prod(herm(Xu), Xu);

   // get the eigenmatrix.  Do some dodgy explict restarts.
   int Iterations = 20;
   double Tol = 1E-12;
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, x.QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-12;
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, x.QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   if (trace(LeftEigen).real() < 0)
      LeftEigen = -LeftEigen;

   DEBUG_CHECK(norm_frob(LeftEigen - adjoint(LeftEigen)) < 1e-10)
      (norm_frob(LeftEigen - adjoint(LeftEigen)));

   // same for the right eigenvector
   MatrixOperator C_left = x.C_right;
   LinearWavefunction PsiR x.Psi;
   C_left = right_orthogonalize(x.Psi, C_left);
   MatrixOperator Ru = C_old_inverse * C_left;
   PsiR.set_front(prod(Ru, PsiR.get_front()));

   // initialize the guess eigenvector
   MatrixOperator RightEigen = delta_shift(scalar_prod(Ru, herm(Ru)), adjoint(x.QShift));

   //   DEBUG_TRACE(EigenvaluesHermitian(RightEigen));
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));

   // get the eigenmatrix
   Iterations = 20; Tol = 1E-12;
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, x.QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-12;
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, x.QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
      //      DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   }

   //   DEBUG_TRACE(EigenvaluesHermitian(RightEigen));
   //   DEBUG_TRACE(EigenvaluesHermitian(LeftEigen));

   if (trace(RightEigen).real() < 0)
      RightEigen = -RightEigen;

   DEBUG_CHECK(norm_frob(RightEigen - adjoint(RightEigen)) < 1e-10)
      (norm_frob(RightEigen - adjoint(RightEigen)));

   // Do the Cholesky factorization of the eigenmatrices.
   // Actually a SVD is much more stable here
   //   MatrixOperator A = CholeskyFactorizeUpper(LeftEigen);
   //   MatrixOperator B = CholeskyFactorizeUpper(RightEigen);
   MatrixOperator A = SingularFactorize(LeftEigen);
   MatrixOperator B = SingularFactorize(RightEigen);

   // LeftEigen == scalar_prod(herm(A), A)
   // RightEigen == scalar_prod(herm(B), B)

   DEBUG_CHECK(norm_frob(LeftEigen - scalar_prod(herm(A), A)) < 1e-10)
      (norm_frob(LeftEigen - scalar_prod(herm(A), A)));
   DEBUG_CHECK(norm_frob(RightEigen - scalar_prod(herm(B), B)) < 1e-10)
      (norm_frob(RightEigen - scalar_prod(herm(B), B)));

   // transform C_old to orthogonal form
   MatrixOperator new_C_old = delta_shift(A, x.QShift) * x.C_old * herm(delta_shift(B, x.QShift));
   // and do the singular value decomposition, updating x.C_old to be diagonal
   MatrixOperator U, Vh;
   SingularValueDecomposition(new_C_old, U, x.C_old, Vh);

   // shift to the new basis
   A = herm(delta_shift(U, adjoint(x.QShift))) * A;
   //B = B * herm(Vh);
   B = delta_shift(Vh, adjoint(x.QShift)) * B;

   // normalize
   x.C_old *= 1.0 / norm_frob(x.C_old);
   //   DEBUG_TRACE(EigenvaluesHermitian(x.C_old));

   // apply A to L   
   //x.C_right = x.C_right * herm(B);
   x.C_right = right_orthogonalize(x.Psi, x.C_right);
   x.C_right = delta_shift(A, x.QShift) * x.C_right;
   x.C_right = left_orthogonalize(x.C_right, x.Psi);
   x.C_right = x.C_right * herm(B);

   // shift C_right to the diagonal basis.  This is a bit of a hack. 
   //TRACE(SingularValues(x.C_old));
   //Vh = x.C_right * InvertDiagonal(x.C_old, 1e-9);
   //x.Psi.set_back(prod(x.Psi.get_back(), Vh));
   //x.C_right = x.C_old;

   //MatrixOperator Ctemp = x.C_right;
   //SingularValueDecomposition(Ctemp, U, x.C_right, Vh);
   //TRACE(Vh);
   //x.Psi.set_back(prod(x.Psi.get_back(), U));
   //x.C_right = x.C_right*Vh;
   
   // normalize
   x.C_right *= 1.0 / norm_frob(x.C_right);

   //   TRACE(x.C_right);
   //   DEBUG_TRACE(SingularValues(x.C_right));
}

#else
#if 1

void orthogonalize(InfiniteWavefunction& x)
{
   // firstly, get the eigenmatrices of the left and right orthogonality superoperators
   // We use the left-orthogonalized matrices, which has the fixed point A^\dagger E A = E
   // where E is the identity operator.  For the right-side, we use the same orthogonalization
   // which gives the fixed point A \rho A^\dagger = \rho, where \rho is the density matrix.
   // We then diagonalize rho, and write rho = U \Lambda \Lambda U^\dagger

   // Incorporate C_old and C_right into the A-matrices
   LinearWavefunction PsiL = x.Psi;
   StateComponent Last = PsiL.get_back();  // right-most A-matrix
   MatrixOperator C_old_inverse = InvertDiagonal(x.C_old, InverseTol);
   //   TRACE(SingularValues(C_old_inverse));
   MatrixOperator Xu = x.C_right * delta_shift(C_old_inverse, adjoint(x.QShift));
   Last = prod(Last, Xu);
   PsiL.set_back(Last);

   // initialize LeftEigen to a guess eigenvector.  Since L satisfies the left orthogonality
   // constraint (except for the final matrix), we can do one iteration beyond the identity
   // and intialize it to herm(Xu) * Xu
   MatrixOperator LeftEigen = scalar_prod(herm(Xu), Xu);

   // get the eigenmatrix.  Do some dodgy explict restarts.
   int Iterations = 20;
   double Tol = 1E-14;
   LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, x.QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-14;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, x.QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   CHECK(EtaL.real() > 0)("Eigenvalue must be positive");

   TRACE(EtaL);

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

   MatrixOperator A = delta_shift(D * U, x.QShift);
   MatrixOperator AInv = adjoint(U) * DInv;

   A = left_orthogonalize(A, PsiL);
   PsiL.set_back(prod(PsiL.get_back(), A*AInv));

   // same for the right eigenvector

   // initialize the guess eigenvector
   MatrixOperator RightEigen = scalar_prod(herm(Xu), Xu); //delta_shift(scalar_prod(Xu, herm(Xu)), adjoint(x.QShift));

   // get the eigenmatrix
   Iterations = 20; Tol = 1E-14;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, x.QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-14;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiL, x.QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }
   TRACE(EtaR);

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

   A = delta_shift(U, x.QShift);
   AInv = adjoint(U);

   //   PsiL = inject_left_old_interface(A, PsiL);
   //   PsiL.set_back(prod(PsiL.get_back(), A * AInv));

   AInv = right_orthogonalize(PsiL, AInv);
   PsiL.set_front(prod(A*AInv, PsiL.get_front()));

   //   PsiL.set_back(prod(PsiL.get_back(), adjoint(U)));
   //   PsiL.set_front(prod(delta_shift(U, x.QShift), PsiL.get_front()));

   // orthonormalize each component of PsiL
   MatrixOperator I = MatrixOperator::make_identity(PsiL.Basis1());
   I = left_orthogonalize(I, PsiL);
   // now I should be unitary (or even an identity operator).  Normalize it and fold it back into the last component
   I *= 1.0 / std::sqrt(EtaR.real());
   PsiL.set_back(prod(PsiL.get_back(), I));

   // normalize
   D *= 1.0 / norm_frob(D);

   x.Psi = PsiL;
   x.C_right = D;
   x.C_old = delta_shift(D, x.QShift);
}

#else

void orthogonalize(InfiniteWavefunction& x)
{
   // firstly, get the eigenmatrices of the left and right orthogonality superoperators

   // Incorporate C_old and C_right into the A-matrices
   LinearWavefunction PsiL = x.Psi;
   StateComponent Last = PsiL.get_back();  // right-most A-matrix
   //   TRACE(SingularValues(x.C_old));
   MatrixOperator C_old_inverse = InvertDiagonal(x.C_old, InverseTol);
   MatrixOperator Xu = x.C_right * delta_shift(C_old_inverse, adjoint(x.QShift));
   Last = prod(Last, Xu);
   PsiL.set_back(Last);

   // initialize LeftEigen to a guess eigenvector.  Since L satisfies the left orthogonality
   // constraint (except for the final matrix), we can do one iteration beyond the identity
   // and intialize it to herm(Xu) * Xu
   MatrixOperator LeftEigen = scalar_prod(herm(Xu), Xu);

   // get the eigenmatrix.  Do some dodgy explict restarts.
   int Iterations = 20;
   double Tol = 1E-14;
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, x.QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-14;
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMultiply(PsiL, x.QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   if (trace(LeftEigen).real() < 0)
      LeftEigen = -LeftEigen;

   DEBUG_CHECK(norm_frob(LeftEigen - adjoint(LeftEigen)) < 1e-10)
      (norm_frob(LeftEigen - adjoint(LeftEigen)));

   // same for the right eigenvector
   MatrixOperator C_left = x.C_right;
   LinearWavefunction PsiR = x.Psi;
   C_left = right_orthogonalize(x.Psi, C_left);
   MatrixOperator Ru = C_old_inverse * C_left;
   PsiR.set_front(prod(Ru, PsiR.get_front()));

   // initialize the guess eigenvector
   MatrixOperator RightEigen = delta_shift(scalar_prod(Ru, herm(Ru)), adjoint(x.QShift));

   //   DEBUG_TRACE(EigenvaluesHermitian(RightEigen));
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));

   // get the eigenmatrix
   Iterations = 20; Tol = 1E-15;
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, x.QShift), 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-15;
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMultiply(PsiR, x.QShift), Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
      //      DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
      //      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   }

   //   DEBUG_TRACE(EigenvaluesHermitian(RightEigen));
   //   DEBUG_TRACE(EigenvaluesHermitian(LeftEigen));

   if (trace(RightEigen).real() < 0)
      RightEigen = -RightEigen;

   DEBUG_CHECK(norm_frob(RightEigen - adjoint(RightEigen)) < 1e-10)
      (norm_frob(RightEigen - adjoint(RightEigen)));

   // Do the Cholesky factorization of the eigenmatrices.
   // Actually a SVD is much more stable here
   //   MatrixOperator A = CholeskyFactorizeUpper(LeftEigen);
   //   MatrixOperator B = CholeskyFactorizeUpper(RightEigen);

   MatrixOperator D = LeftEigen;
   MatrixOperator U = DiagonalizeHermitian(D);
   //   TRACE(SingularValues(D));
   D = SqrtDiagonal(D, OrthoTol);
   MatrixOperator DInv = InvertDiagonal(D, OrthoTol);

   // LeftEigen = triple_prod(U, D*D, herm(U))
   DEBUG_CHECK(norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)) < 1e-10)
      (norm_frob(LeftEigen - triple_prod(herm(U), D*D, U)));

   MatrixOperator A = delta_shift(D * U, x.QShift);
   MatrixOperator AInv = adjoint(U) * DInv;

   //   MatrixOperator A = SingularFactorize(LeftEigen);
   MatrixOperator B = SingularFactorize(RightEigen);

   // LeftEigen == scalar_prod(herm(A), A)
   // RightEigen == scalar_prod(herm(B), B)

   DEBUG_CHECK(norm_frob(LeftEigen - scalar_prod(herm(A), A)) < 1e-10)
      (norm_frob(LeftEigen - scalar_prod(herm(A), A)));
   DEBUG_CHECK(norm_frob(RightEigen - scalar_prod(herm(B), B)) < 1e-10)
      (norm_frob(RightEigen - scalar_prod(herm(B), B)));

   // transform C_old to orthogonal form
   MatrixOperator new_C_old = delta_shift(A, x.QShift) * x.C_old * herm(delta_shift(B, x.QShift));
   // and do the singular value decomposition, updating x.C_old to be diagonal
   MatrixOperator Vh;
   SingularValueDecomposition(new_C_old, U, x.C_old, Vh);
   
   // normalize
   x.C_old *= 1.0 / norm_frob(x.C_old);

   // shift to the new basis
   A = herm(delta_shift(U, adjoint(x.QShift))) * A;
   AInv = AInv * U;

   //B = B * herm(Vh);
   //   B = delta_shift(Vh, adjoint(x.QShift)) * B;

   // Transform the wavefunction
   x.C_right = delta_shift(x.C_old, adjoint(x.QShift));
   x.Psi = PsiL;
   x.Psi.set_front(prod(A, x.Psi.get_front()));
   x.Psi.set_back(prod(x.Psi.get_back(), AInv));

#if 1 // orthogonalize the internal matrices
   MatrixOperator I = MatrixOperator::make_identity(x.Psi.Basis1());
   I = left_orthogonalize(I, x.Psi);
   x.Psi.set_back(prod(x.Psi.get_back(), I));
#endif
}

#endif
#endif

std::complex<double> overlap(InfiniteWavefunction const& x, InfiniteWavefunction const& y,
                             QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, bool Verbose)
{
   //   TRACE(x.C_old);
   CHECK_EQUAL(x.QShift, y.QShift)("The wavefunctions must have the same quantum number per unit cell");
   MatrixOperator x_unit =  x.C_right * delta_shift(InvertDiagonal(x.C_old, InverseTol), adjoint(x.QShift));
   MatrixOperator y_unit =  y.C_right * delta_shift(InvertDiagonal(y.C_old, InverseTol), adjoint(y.QShift));
   
   LinearWavefunction xPsi = x.Psi;
   xPsi.set_back(prod(xPsi.get_back(), x_unit));

   LinearWavefunction yPsi = y.Psi;
   yPsi.set_back(prod(yPsi.get_back(), y_unit));

   MatrixOperator Init = MakeRandomMatrixOperator(x.Psi.Basis1(), y.Psi.Basis1(), Sector);

   int Iterations = Iter;
   int TotalIterations = 0;
   double MyTol = Tol;
   std::complex<double> Eta = LinearSolvers::Arnoldi(Init, 
                                                     GeneralizedLeftMultiply(xPsi, yPsi, x.QShift), 
                                                     Iterations, MyTol, LinearSolvers::LargestMagnitude, false, Verbose);
   TotalIterations += Iterations;
   DEBUG_TRACE(Eta)(Iterations);

   while (MyTol < 0)
   {
      if (Verbose)
         std::cerr << "Restarting Arnoldi, eta=" << Eta << ", Tol=" << -MyTol << '\n';
      Iterations = Iter;
      MyTol = Tol;
      Eta = LinearSolvers::Arnoldi(Init, GeneralizedLeftMultiply(xPsi, yPsi, x.QShift), 
				   Iterations, MyTol, LinearSolvers::LargestMagnitude, false, Verbose);
      TotalIterations += Iterations;
      DEBUG_TRACE(Eta)(Iterations);
   }
   if (Verbose)
      std::cerr << "Converged.  TotalIterations=" << TotalIterations
                << ", Tol=" << MyTol << '\n';

   return Eta;
}
