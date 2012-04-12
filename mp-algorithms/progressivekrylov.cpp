// -*- C++ -*- $Id$

#include "progressivekrylov.h"
#include "matrixproduct/wavefunc-utils.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/exponential.h"
#include "linearalgebra/matrix_utility.h"
#include "matrixproduct/matrixproduct-sum.h"

typedef std::complex<double> complex;
using LinearAlgebra::range;

PStream::opstream& operator<<(PStream::opstream& out, KrylovSolver const& s)
{
   return out << s.Krylov
	      << s.H
	      << s.kn1_H_kn
	      << s.ki_kj
	      << s.sub_H
	      << s.Ident
      ;
}

PStream::ipstream& operator>>(PStream::ipstream& in, KrylovSolver& s)
{
   return in >> s.Krylov
	     >> s.H
	     >> s.kn1_H_kn
	     >> s.ki_kj
	     >> s.sub_H
	     >> s.Ident
      ;
}

KrylovSolver::KrylovSolver(SplitOperator const& H_, SplitOperator const& H2_, CenterWavefunction const& Psi_)
   : Krylov(1, Psi_), NumConvergedKrylov(1), H(H_), H2(H2_), kn1_H_kn(), ki_kj(1,1), 
     sub_I(1, 1, 1.0),
     sub_H(1, 1, expectation(Psi_, H_, Psi_).real()),
     sub_H2(1, expectation(Psi_, H2_, Psi_).real()),
     Ident(H_.GetSymmetryList())
{
   //Krylov[0].normalize();
   sub_L = sub_I;
   CholeskyFactorizeLower(sub_L);
   zero_upper_triangular(sub_L);
   sub_Linv = sub_L;
   InvertLowerTriangular(sub_Linv);
   ortho_H = sub_Linv * sub_H * herm(sub_Linv);
}

KrylovSolver::KrylovSolver(SplitOperator const& H_, SplitOperator const& H2_, CenterWavefunction const& Psi_, double Energy)
   : Krylov(1, Psi_),  NumConvergedKrylov(1), H(H_), kn1_H_kn(), ki_kj(1,1), 
     sub_I(1, 1, 1.0),
     sub_H(1, 1, Energy),
     sub_H2(1, expectation(Psi_, H2_, Psi_).real()),
     Ident(H_.GetSymmetryList())
{
   //   Krylov[0].normalize();
   sub_L = sub_I;
   CholeskyFactorizeLower(sub_L);
   zero_upper_triangular(sub_L);
   sub_Linv = sub_L;
   InvertLowerTriangular(sub_Linv);
   ortho_H = sub_Linv * sub_H * herm(sub_Linv);
}

KrylovSolver::KrylovSolver(SplitOperator const& H_, SplitOperator const& H2_, CenterWavefunction const& Psi_, 
			   double Energy, double ExpectationH2)
   : Krylov(1, Psi_),  NumConvergedKrylov(1), H(H_), kn1_H_kn(), ki_kj(1,1), 
     sub_I(1, 1, 1.0),
     sub_H(1, 1, Energy),
     sub_H2(1, ExpectationH2),
     Ident(H_.GetSymmetryList())
{
   Krylov[0].RotateToNormalForm();
   H.RotateToNormalForm();
   H2.RotateToNormalForm();

   //Krylov[0].normalize();
   sub_L = sub_I;
   CholeskyFactorizeLower(sub_L);
   zero_upper_triangular(sub_L);
   sub_Linv = sub_L;
   InvertLowerTriangular(sub_Linv);
   ortho_H = sub_Linv * sub_H * herm(sub_Linv);
}

void
KrylovSolver::AddKrylovVector(CenterWavefunction const& K)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size());

   Krylov.push_back(K);
   while (Krylov.back().LeftSize() > Krylov[0].LeftSize())
      Krylov.back().RotateLeft();
   while (Krylov.back().RightSize() > Krylov[0].RightSize())
      Krylov.back().RotateRight();
   while (H.LeftSize() > Krylov[0].LeftSize())
      H.RotateLeft();
   while (H.RightSize() > Krylov[0].RightSize())
      H.RotateRight();

   InitializeSuperblockStack(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   Matrix<TransformOperator> New_ki_kj(Krylov.size(), Krylov.size());
   New_ki_kj(range(0,Krylov.size()-1), range(0, Krylov.size()-1)) = ki_kj;
   ki_kj = New_ki_kj;
   for (unsigned i = 0; i < Krylov.size()-1; ++i)
   {
      InitializeTransformStack(ki_kj(i, NumConvergedKrylov), Krylov[i], Krylov[NumConvergedKrylov]);
   }
}

void KrylovSolver::ShiftRightAndExpand(bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   // Hamiltonian operator
   H.RotateRight();
   // Wavefunctions
   for (unsigned i = jStart; i < Krylov.size()-1; ++i) // this is OK, Krylov.size() is always >= 1
   {
      Krylov[i].RotateRight();
   }
   Krylov.back().RotateRightExpand();
   // matrix elements of < Krylov[last] | H | Krylov[last-1] >
   SuperblockStackRotateRight(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
	 TransformStackRotateRight(ki_kj(i,j), Krylov[i], Krylov[j]);
      }
   }
   this->DebugCheckBasis();
}

void KrylovSolver::ShiftLeftAndExpand(bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   // Hamiltonian operator
   H.RotateLeft();
   // Wavefunctions
   for (unsigned i = jStart; i < Krylov.size()-1; ++i) // this is OK, Krylov.size() is always >= 1
   {
      Krylov[i].RotateLeft();
   }
   Krylov.back().RotateLeftExpand();
   // matrix elements of < Krylov[last] | H | Krylov[last-1] >
   SuperblockStackRotateLeft(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
	 TransformStackRotateLeft(ki_kj(i,j), Krylov[i], Krylov[j]);
      }
   }
   this->DebugCheckBasis();
}

void KrylovSolver::ExpandLeft(bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   Krylov.back().Left() = prod(Krylov.back().Left(), Krylov.back().Center());
   Krylov.back().Center() = ExpandBasis2(Krylov.back().Left());
   SuperblockStackUpdateLeft(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      TransformStackUpdateLeft(ki_kj(i, NumConvergedKrylov), Krylov[i], Krylov[NumConvergedKrylov]);
   }
   this->DebugCheckBasis();
}

void KrylovSolver::ExpandRight(bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   Krylov.back().Right() = prod(Krylov.back().Center(), Krylov.back().Right());
   Krylov.back().Center() = ExpandBasis1(Krylov.back().Right());
   SuperblockStackUpdateRight(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      TransformStackUpdateRight(ki_kj(i, Krylov.size()-1), Krylov[i], Krylov.back());
   }
   this->DebugCheckBasis();
}

#if 0
Matrix<complex> KrylovSolver::OrthogonalityMatrix() const
{
   Matrix<complex> Result(Krylov.size(), Krylov.size());
   for (unsigned i = 0; i < Krylov.size(); ++i)
   {
      Result(i,i) = norm_frob_sq(Krylov[i].Center());
      for (unsigned j = i+1; j < Krylov.size(); ++j)
      {
	 Result(i,j) = inner_prod(triple_prod(herm(ki_kj(j,i).Left()), 
					      Krylov[j].Center(), 
					      ki_kj(j,i).Right()),
				  Krylov[i].Center());
	 Result(j,i) = Result(i,j);
      }
   }
   return Result;
}
#endif

void
KrylovSolver::Solve(bool FullOrthogonalization)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrthogonalization || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   // Psi[k+1] = H|Psi[k]> - sum_{j<=k} <Psi[j] | H | Psi[k] > | Psi[j] >
   Krylov.back().Center() = operator_prod(conj(H.Center()),
					  kn1_H_kn.Left(), 
					  Krylov[Krylov.size()-2].Center(), 
					  herm(kn1_H_kn.Right()));
   for (unsigned j = jStart; j < NumConvergedKrylov; ++j)
   {
      Krylov.back().Center() -= sub_H(j, NumConvergedKrylov-1) * 
	 triple_prod(herm(ki_kj(j, NumConvergedKrylov).Left()),
		     Krylov[j].Center(),
		     ki_kj(j, NumConvergedKrylov).Right());
   }
#if 0
   // if the local H is zero, add a bit of Psi[k] itself...
   if (norm_frob_sq(Krylov.back().Center()) < 1E-14)
   {
      Krylov.back().Center() += 1E-10 * 
         triple_prod(herm(ki_kj(NumConvergedKrylov-1, NumConvergedKrylov).Left()),
		     Krylov[NumConvergedKrylov-1].Center(),
		     ki_kj(NumConvergedKrylov-1, NumConvergedKrylov).Right());
   }
#endif
}

double
KrylovSolver::Variance(bool FullOrtho) const
{
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   // <k+1|k+1>
   double NormPart = norm_frob_sq(Krylov.back().Center());
   // <k|H^2|k>
   double H2Part = sub_H2[NumConvergedKrylov-1];

   // sum_{j<=k} sum_{l<=k} <j|H|k><k|H|l><l|j>
   complex SumPart = 0.0;
   for (unsigned j = jStart; j < NumConvergedKrylov; ++j)
   {
      for (unsigned l = jStart; l < NumConvergedKrylov; ++l)
      {
	 SumPart += sub_H(j, NumConvergedKrylov-1) * sub_H(NumConvergedKrylov-1, l) * sub_I(l,j);
      }
   }

   // <k+1|H|k> + h.c.
   double k1HkPart = 2.0 * inner_prod(Krylov.back().Center(), 
				      operator_prod(conj(H.Center()),
						    kn1_H_kn.Left(), 
						    Krylov[Krylov.size()-2].Center(), 
						    herm(kn1_H_kn.Right()))).real();

   // sum_{j<=k} <j|H|k><k+1|j> + h.c.
   double HPart = 0;
   for (unsigned j = jStart; j < NumConvergedKrylov; ++j)
   {
      HPart += 2.0 * (sub_H(j, NumConvergedKrylov-1) * 
		      inner_prod(Krylov[Krylov.size()-1].Center(),
				 triple_prod(herm(ki_kj(j, Krylov.size()-1).Left()),
					     Krylov[j].Center(),
					     ki_kj(j, Krylov.size()-1).Right()))).real();
   }

   // sum_{j<=k} <j|H|k> <k|H|j> + h.c.
   double HjkPart = 0;
   for (unsigned j = jStart; j < NumConvergedKrylov; ++j)
   {
      HjkPart += 2.0 * LinearAlgebra::norm_2_sq(sub_H(j,NumConvergedKrylov-1));
   }

   return (NormPart + H2Part + SumPart.real() + HPart - k1HkPart - HjkPart)
      / (0.5 * (NormPart + H2Part + SumPart.real() - HPart + k1HkPart - HjkPart));
}

/*
double KrylovSolver::SubMatrixElement() const
{
   // This function is incorrect, do not use it.

   // <k|H^2|k>
   double H2Part = sub_H2[NumConvergedKrylov-1] 
      / LinearAlgebra::norm_2_sq(sub_L(NumConvergedKrylov-1, NumConvergedKrylov-1));

   // sum_{j<=k} sum_{l<=k} <j|H|k><k|H|l><l|j>
   complex SumPart = 0.0;
   for (unsigned j = 0; j < NumConvergedKrylov; ++j)
   {
      SumPart += ortho_H(j, NumConvergedKrylov-1) * ortho_H(NumConvergedKrylov-1, j);
   }

   // sum_{j<=k} <j|H|k> <k|H|j> + h.c.
   double HjkPart = 0;
   for (unsigned j = 0; j < NumConvergedKrylov; ++j)
   {
      HjkPart += 2.0 * LinearAlgebra::norm_2_sq(ortho_H(j,NumConvergedKrylov-1));
   }

   TRACE(sub_H2[0])(ortho_H)(sub_H)(sub_I);

   TRACE(H2Part)(SumPart.real())(HjkPart);

   TRACE(H2Part + SumPart.real() - HjkPart);

   return std::sqrt(H2Part + SumPart.real() - HjkPart);
}
*/

Vector<complex> KrylovSolver::Exponentiate(complex Tau, double SubElement) const
{
   Matrix<complex> M(size1(ortho_H)+1, size2(ortho_H)+1, 0.0);
   M(range(0,size1(ortho_H)), range(0, size2(ortho_H))) = ortho_H;
   // Set the (m+1,m) component.  Remember e_m is invariant under the LU transformation.
   M(size1(ortho_H), size2(ortho_H)-1) = SubElement;
   M = Tau * M;
   TRACE(M)(Tau)(sub_L);
   Matrix<complex> Res = LinearAlgebra::Exponentiate(1.0, M);
   //TRACE(Res);
   //TRACE(M);
   Vector<complex> e1 = direct_sum(conj(sub_Linv(0, LinearAlgebra::all)), 0.0);
   //TRACE(e1);
   Vector<complex> Result(e1.size(), 0.0);
   for (unsigned i = 0; i < size1(Res); ++i)
   {
      for (unsigned j = 0; j < size2(Res); ++j)
      {
         Result[i] += Res(i,j) * e1[j];
      }
   }
   return Result;
}

Vector<complex> KrylovSolver::Exponentiate(complex Tau) const
{
   Matrix<complex> M = ortho_H;
   M = Tau * M;
   Matrix<complex> Res = LinearAlgebra::Exponentiate(1.0, M);
   //TRACE(Res);
   //TRACE(M);
   Vector<complex> e1 = conj(sub_Linv(0, LinearAlgebra::all));
   //TRACE(e1);
   Vector<complex> Result(e1.size(), 0.0);
   for (unsigned i = 0; i < size1(Res); ++i)
   {
      for (unsigned j = 0; j < size2(Res); ++j)
      {
         Result[i] += Res(i,j) * e1[j];
      }
   }
   return Result;
}

void KrylovSolver::FixKrylovVector(bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;

   // Normalize the new Krylov vector
   Krylov.back().Normalize();

   // New matrix elements for sub_H
   Matrix<std::complex<double> > new_sub_H(Krylov.size(), Krylov.size());
   new_sub_H(range(0, Krylov.size()-1), range(0, Krylov.size()-1)) = sub_H;
   sub_H = new_sub_H;

   sub_H(Krylov.size()-1, Krylov.size()-1) = expectation(Krylov.back(), H, Krylov.back()).real();
   sub_H(Krylov.size()-1, Krylov.size()-2) = inner_prod(Krylov.back().Center(), 
							operator_prod(conj(H.Center()),
								      kn1_H_kn.Left(), 
								      Krylov[Krylov.size()-2].Center(), 
								      herm(kn1_H_kn.Right())));
   DEBUG_CHECK(LinearAlgebra::equal(sub_H(Krylov.size()-1, Krylov.size()-2), 
		     expectation(Krylov.back(), H, Krylov[Krylov.size()-2]),
		     1E-5))
      (sub_H(Krylov.size()-1, Krylov.size()-2))
      (expectation(Krylov.back(), H, Krylov[Krylov.size()-2]));
   sub_H(Krylov.size()-2, Krylov.size()-1) = conj(sub_H(Krylov.size()-1, Krylov.size()-2));

   for (unsigned j = 0; j < Krylov.size()-2; ++j)
   {
      sub_H(Krylov.size()-1, j) = expectation(Krylov[Krylov.size()-1], H, Krylov[j]);
      sub_H(j, Krylov.size()-1) = conj(sub_H(Krylov.size()-1, j));
   }

   // New matrix elements for sub_I
   Matrix<std::complex<double> > new_sub_I(Krylov.size(), Krylov.size());
   new_sub_I(range(0, Krylov.size()-1), range(0, Krylov.size()-1)) = sub_I;
   sub_I = new_sub_I;

   sub_I(Krylov.size()-1, Krylov.size()-1) = 1.0;
   for (unsigned j = jStart; j < Krylov.size()-1; ++j)
   {
      sub_I(Krylov.size()-1, j) = inner_prod(Krylov[Krylov.size()-1].Center(),
					     triple_prod(herm(ki_kj(j, Krylov.size()-1).Left()),
							 Krylov[j].Center(),
							 ki_kj(j, Krylov.size()-1).Right()));
      sub_I(j, Krylov.size()-1) = conj(sub_I(Krylov.size()-1, j));

      DEBUG_CHECK(LinearAlgebra::equal(sub_I(Krylov.size()-1, j), 
			overlap(Krylov[Krylov.size()-1], Krylov[j]),
			1E-5))
		  (sub_I(Krylov.size()-1, j))
		  (overlap(Krylov[Krylov.size()-1], Krylov[j]));
   }
   for (unsigned j = 0; j < jStart; ++j)
   {
      sub_I(Krylov.size()-1, j) = overlap(Krylov[Krylov.size()-1], Krylov[j]);
      sub_I(j, Krylov.size()-1) = conj(sub_I(Krylov.size()-1, j));
   }

   // New expectation value < K | H^2 | K >
   sub_H2.push_back(expectation(Krylov.back(), H2, Krylov.back()).real());

   // Update the Cholesky factorization and the orthogonalized H
   sub_L = sub_I;
   CholeskyFactorizeLower(sub_L);
   zero_upper_triangular(sub_L);
   sub_Linv = sub_L;
   InvertLowerTriangular(sub_Linv);
   ortho_H = sub_Linv * sub_H * herm(sub_Linv);
   TRACE(sub_Linv);

   // Finally, our new Krylov vector is converged
   ++NumConvergedKrylov;
}

void
KrylovSolver::TruncateLeft(int MaxStates, double MixFactor, bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);

   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   MatrixOperator Rho = scalar_prod(Krylov.back().Center(), herm(Krylov.back().Center()));

   // Get the trace and normalize Rho.  But watch out if Rho is very small.  In that case,
   // we rely on the mix factor.
   double Tr = LinearAlgebra::norm_2(trace(Rho));
   if (Tr > std::numeric_limits<double>::epsilon()*10)
   {
      Rho *= (1.0 / Tr);
   }
      
   if (MixFactor != 0)
   {
      MatrixOperator RhoLast = Krylov[Krylov.size()-2].Center();
      RhoLast = scalar_prod(RhoLast, herm(RhoLast));
      MatrixOperator Correction =
         operator_prod(kn1_H_kn.Left(), 
                       RhoLast, 
                       herm(kn1_H_kn.Left()));
      Correction *= (MixFactor / trace(Correction));
      Rho += Correction;
   }      

   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);
   Krylov.back().Left() = prod(Krylov.back().Left(), herm(U));
   Krylov.back().Center() = prod(U, Krylov.back().Center(), Ident);
   kn1_H_kn.Left() = prod(U, kn1_H_kn.Left());
   //   SuperblockStackUpdateLeft(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      unsigned const j = Krylov.size()-1;
      ki_kj(i,j).Left() = ki_kj(i,j).Left() * herm(U);
      //TransformStackUpdateLeft(ki_kj(i,j), Krylov[i], Krylov[j]);
   }
   this->DebugCheckBasis();
}

void
KrylovSolver::TruncateRight(int MaxStates, double MixFactor, bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   MatrixOperator Rho = scalar_prod(herm(Krylov.back().Center()), Krylov.back().Center());

   // Get the trace and normalize Rho.  But watch out if Rho is very small.  In that case,
   // we rely on the mix factor.
   double Tr = LinearAlgebra::norm_2(trace(Rho));
   if (Tr > std::numeric_limits<double>::epsilon()*10)
   {
      Rho *= (1.0 / Tr);
   }

   if (MixFactor != 0)
   {
      MatrixOperator RhoLast = Krylov[Krylov.size()-2].Center();
      RhoLast = scalar_prod(herm(RhoLast), RhoLast);
         MatrixOperator Correction =
            operator_prod(kn1_H_kn.Right(),
                          RhoLast,
                          herm(kn1_H_kn.Right()));
      Correction *= (MixFactor / trace(Correction));
      Rho += Correction;
   }

   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);
   Krylov.back().Right() = prod(U, Krylov.back().Right());
   Krylov.back().Center() = Krylov.back().Center() * herm(U);
   SuperblockStackUpdateRight(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      unsigned const j = Krylov.size()-1;
      TransformStackUpdateRight(ki_kj(i,j), Krylov[i], Krylov[j]);
   }
   this->DebugCheckBasis();
}

TruncationInfo
KrylovSolver::TruncateLeft(int MinStates, int MaxStates, double MinTruncation, double MixFactor, 
                           bool FullMixing, bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   MatrixOperator Rho = scalar_prod(Krylov.back().Center(), herm(Krylov.back().Center()));
   MatrixOperator RhoPsi = Rho;

   // Get the trace and normalize Rho.  But watch out if Rho is very small.  In that case,
   // we rely on the mix factor.
   double Tr = LinearAlgebra::norm_2(trace(Rho));
   if (Tr > std::numeric_limits<double>::epsilon()*10)
   {
      Rho *= (1.0 / Tr);
   }

   if (MixFactor != 0)
   {
      MatrixOperator RhoLast = Krylov[Krylov.size()-2].Center();
      RhoLast = scalar_prod(RhoLast, herm(RhoLast));
      MatrixOperator Correction =
         operator_prod(kn1_H_kn.Left(), 
                       RhoLast, 
                       herm(kn1_H_kn.Left()));
      Correction *= (MixFactor / trace(Correction));
      Rho += Correction;
   }      

   DensityMatrix<MatrixOperator> DM = FullMixing ?  DensityMatrix<MatrixOperator>(Rho)
      : DensityMatrix<MatrixOperator>(Rho, RhoPsi);

   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                                                                           MinStates,
                                                                                           MaxStates,
                                                                                           MinTruncation,
                                                                                           Info));
   Krylov.back().Left() = prod(Krylov.back().Left(), herm(U));
   Krylov.back().Center() = U * Krylov.back().Center();
   kn1_H_kn.Left() = prod(U, kn1_H_kn.Left());
   //   SuperblockStackUpdateLeft(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      unsigned const j = Krylov.size()-1;
      ki_kj(i,j).Left() = ki_kj(i,j).Left() * herm(U);
      //TransformStackUpdateLeft(ki_kj(i,j), Krylov[i], Krylov[j]);
   }
   this->DebugCheckBasis();
   return Info;
}

TruncationInfo
KrylovSolver::TruncateRight(int MinStates, int MaxStates, double MinTruncation, double MixFactor,
                            bool FullMixing, bool FullOrtho)
{
   CHECK_EQUAL(NumConvergedKrylov, Krylov.size()-1);
   unsigned const jStart = (FullOrtho || NumConvergedKrylov < 2) ? 0 : NumConvergedKrylov-2;
   MatrixOperator Rho = scalar_prod(herm(Krylov.back().Center()), Krylov.back().Center());
   MatrixOperator RhoPsi = Rho;

   // Get the trace and normalize Rho.  But watch out if Rho is very small.  In that case,
   // we rely on the mix factor.
   double Tr = LinearAlgebra::norm_2(trace(Rho));
   if (Tr > std::numeric_limits<double>::epsilon()*10)
   {
      Rho *= (1.0 / Tr);
   }

   if (MixFactor != 0)
   {
     MatrixOperator RhoLast = Krylov[Krylov.size()-2].Center();
      RhoLast = scalar_prod(herm(RhoLast), RhoLast);
         MatrixOperator Correction =
            operator_prod(kn1_H_kn.Right(),
                          RhoLast,
                          herm(kn1_H_kn.Right()));
      Correction *= (MixFactor / trace(Correction));
      Rho += Correction;
   }

   DensityMatrix<MatrixOperator> DM = FullMixing ?  DensityMatrix<MatrixOperator>(Rho)
      : DensityMatrix<MatrixOperator>(Rho, RhoPsi);

   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                                                                           MinStates,
                                                                                           MaxStates,
                                                                                           MinTruncation,
                                                                                           Info));
   Krylov.back().Right() = prod(U, Krylov.back().Right());
   Krylov.back().Center() = Krylov.back().Center() * herm(U);
   SuperblockStackUpdateRight(kn1_H_kn, Krylov[Krylov.size()-1], H, Krylov[Krylov.size()-2]);
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = jStart; i < Krylov.size()-1; ++i)
   {
      unsigned const j = Krylov.size()-1;
      TransformStackUpdateRight(ki_kj(i,j), Krylov[i], Krylov[j]);
   }
   this->DebugCheckBasis();
   return Info;
}

CenterWavefunction 
KrylovSolver::ConstructExpansion(LinearAlgebra::Vector<std::complex<double> > const& n,
                                 StatesInfo const& SInfo) const
{
   std::vector<CenterWavefunction> X(Krylov);
   if (n.size() < X.size()) X.resize(n.size());
   for (unsigned i = 0; i < n.size(); ++i)
   {
      X[i] *= n[i];
   }
   return mpwavefunction_sum(X, SInfo, true);
}

void KrylovSolver::DebugCheckBasis() const
{
   DEBUG_CHECK_EQUAL(kn1_H_kn.Left().Basis2(), Krylov[Krylov.size()-2].Center().Basis1());
   DEBUG_CHECK_EQUAL(kn1_H_kn.Right().Basis2(), Krylov[Krylov.size()-2].Center().Basis2());

   DEBUG_CHECK_EQUAL(kn1_H_kn.Left().Basis1(), Krylov[Krylov.size()-1].Center().Basis1());
   DEBUG_CHECK_EQUAL(kn1_H_kn.Right().Basis1(), Krylov[Krylov.size()-1].Center().Basis2());

#if 0
   // this is broken for the non-full-orthogonalization case
   for (std::size_t i = 0; i < Krylov.size(); ++i)
   {
      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
	 DEBUG_CHECK_EQUAL(ki_kj(i,j).Left().Basis2(), Krylov[j].Center().Basis1());
	 DEBUG_CHECK_EQUAL(ki_kj(i,j).Left().Basis1(), Krylov[i].Center().Basis1());
	 DEBUG_CHECK_EQUAL(ki_kj(i,j).Right().Basis2(), Krylov[j].Center().Basis2());
	 DEBUG_CHECK_EQUAL(ki_kj(i,j).Right().Basis1(), Krylov[i].Center().Basis2());
      }
   }
#endif
}
