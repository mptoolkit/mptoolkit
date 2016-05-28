// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-pwfrg-old.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/wavefunc-utils.h"
#include "matrixproduct/triangularoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"

#include "models/spin-su2.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "tensor/match_basis.h"

namespace prog_opt = boost::program_options;

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(MPStateComponent const& Left_,
		      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }

   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(MPStateComponent const& Left_,
				       MPStateComponent const& Right_)
   : Left(Left_), Right(Right_)
{
}

bool ExpandL = true, ExpandR = true;

LinearAlgebra::Vector<double>
ExtractDiagonal(MatrixOperator const& m)
{
   std::list<double> x;
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      if (iterate_at(m.data(), i, i))
      {
	 for (int j = 0; j < m.Basis1().dim(i); ++j)
	 {
	    x.push_back(m(i,i)(j,j).real());
	 }
      }
      else
      {
	 for (int j = 0; j < m.Basis1().dim(i); ++j)
	 {
	    x.push_back(0.0);
	 }
      }
   }
   return LinearAlgebra::Vector<double>(x.begin(), x.end());
}

void DoIteration(MPOpComponent const& Ham, 
		 MPStateComponent& A1, MPStateComponent& B2,
		 MatrixOperator& OldCenter,
		 MatrixOperator& Center,
		 MPStateComponent& E,
		 MPStateComponent& F)
{
   DEBUG_CHECK_EQUAL(A1.Basis2(), Center.Basis1());
   DEBUG_CHECK_EQUAL(E.Basis2(), Center.Basis1());
   DEBUG_CHECK_EQUAL(E.Basis1(), Center.Basis1());

   DEBUG_CHECK_EQUAL(B2.Basis1(), Center.Basis2());
   DEBUG_CHECK_EQUAL(F.Basis1(), Center.Basis2());
   DEBUG_CHECK_EQUAL(F.Basis2(), Center.Basis2());

   // invert the singular values
   MatrixOperator OldCenterInverse = OldCenter;
   OldCenterInverse *= 0.0;
   for (unsigned i = 0; i < OldCenter.Basis1().size(); ++i)
   {
      if (iterate_at(OldCenter.data(), i, i))
      {
	 for (int j = 0; j < std::min(OldCenter.Basis1().dim(i), OldCenter.Basis2().dim(i)); ++j)
	 {
	    double x = OldCenter(i,i)(j,j).real();
	    if (fabs(x) < 1e-5) 
	       x = 0;
	    else 
	       x = 1.0 / x;
	    OldCenterInverse(i,i)(j,j) = x;
	 }
      }
   }

   // rotate to the right
   MPStateComponent A2 = prod(Center, B2);
   MatrixOperator Cr = TruncateBasis2(A2);

   // rotate to the left
   MPStateComponent B1 = prod(A1, Center);
   MatrixOperator Cl = TruncateBasis1(B1);
   
   // insert the sites
   OldCenter = Center;
   // wavefunction transform

#if 0
   Center = Cr * herm(OldCenterInverse) * Cl;  // wavefunction transformation
#else

   // alternative approach
   MatrixOperator Ur, Dr, Vhr;
   SingularValueDecomposition(Cr, Ur, Dr, Vhr);
   MatrixOperator Ul, Dl, Vhl;
   SingularValueDecomposition(Cl, Ul, Dl, Vhl);

   // take the square root of the singular values
   for (unsigned i = 0; i < Dr.Basis1().size(); ++i)
   {
      if (iterate_at(Dr.data(), i, i))
      {
	 for (int j = 0; j < std::min(Dr.Basis1().dim(i), Dr.Basis2().dim(i)); ++j)
	 {
	    Dr(i,i)(j,j) = std::sqrt(Dr(i,i)(j,j));
	 }
      }
   }

   for (unsigned i = 0; i < Dl.Basis1().size(); ++i)
   {
      if (iterate_at(Dr.data(), i, i))
      {
	 for (int j = 0; j < std::min(Dl.Basis1().dim(i), Dl.Basis2().dim(i)); ++j)
	 {
	    Dl(i,i)(j,j) = std::sqrt(Dl(i,i)(j,j));
	 }
      }
   }

   Center = Ur * Dr * Vhr * Ul * Dl * Vhl;
#endif


   TRACE("transformed wavefunction")(norm_frob(Center));

   // if the wavefunction has ~zero weight, things might get screwy
   if (norm_frob(Center) < 1e-10)
   {
      Center = MakeRandomMatrixOperator(Center.Basis1(), Center.Basis2());
   }

   Center *= 1.0 / norm_frob(Center);    // normalize

   if (ExpandL)
   {
      // expand the left basis
      MatrixOperator U = ExpandBasis2(A2);
      Center = U * Center;
   }
   if (ExpandR)
   {
      // expand the right basis
      MatrixOperator U = ExpandBasis1(B1);
      Center = Center * U;
   }

   E = operator_prod(herm(Ham), herm(A2), E, A2);
   F = operator_prod(Ham, B1, F, herm(B1));
   A1 = A2;
   B2 = B1;
}

int main(int argc, char** argv)
{
   int NumIter = 40;

   SiteBlock Site = CreateSU2SpinSite(0.5);

   int UnitCellSize = 2;

   // initial UnitCellSize-site block
   int MinStates = 1;
   int MaxStates = 20;
   double TruncCutoff = 0;
   double EigenCutoff = -1;

   StatesInfo SInfo;
   SInfo.MinStates = MinStates;
   SInfo.MaxStates = MaxStates;
   SInfo.TruncationCutoff = TruncCutoff;
   SInfo.EigenvalueCutoff = EigenCutoff;
   std::cout << SInfo << '\n';
   
   MpOpTriangular Ham = TriangularTwoSite(-sqrt(3.0)*Site["S"], 
					  Site["S"], 
					  Site["I"].TransformsAs());

   MPStateComponent E = Initial_E(Ham);
   MPStateComponent F = Initial_F(Ham);
   MatrixOperator Center = MatrixOperator::make_identity(E.Basis1());

   MatrixOperator OldCenter = Center;
   typedef std::deque<MPStateComponent> MatListType;

   MatListType A1, A2, B1, B2;

   // make 2 unit cells of left blocks
   A1.push_back(ConstructFromLeftBasis(Site.Basis2().Basis(), Center.Basis1()));
   for (int j = 0; j < UnitCellSize-1; ++j)
   {
      A1.push_back(ConstructFromLeftBasis(Site.Basis2().Basis(), A1.back().Basis2()));
   }

   A2.push_back(ConstructFromLeftBasis(Site.Basis2().Basis(), A1.back().Basis2()));
   for (int j = 0; j < UnitCellSize-1; ++j)
   {
      A2.push_back(ConstructFromLeftBasis(Site.Basis2().Basis(), A2.back().Basis2()));
   }

   // make 2 unit cells of right blocks
   B2.push_front(ConstructFromRightBasis(Site.Basis1().Basis(), Center.Basis2()));
   for (int j = 0; j < UnitCellSize-1; ++j)
   {
      B2.push_front(ConstructFromRightBasis(Site.Basis1().Basis(), B2.front().Basis1()));
   }

   B1.push_front(ConstructFromRightBasis(Site.Basis1().Basis(), B2.front().Basis1()));
   for (int j = 0; j < UnitCellSize-1; ++j)
   {
      B1.push_front(ConstructFromRightBasis(Site.Basis1().Basis(), B1.front().Basis1()));
   }

   for (int j = 0; j < UnitCellSize; ++j)
   {
      E = operator_prod(herm(Ham.data()), herm(A1[j]), E, A1[j]);
   }
   for (int j = 0; j < UnitCellSize; ++j)
   {
      E = operator_prod(herm(Ham.data()), herm(A2[j]), E, A2[j]);
   }

   for (int j = UnitCellSize-1; j >= 0; --j)
   {
      F = operator_prod(Ham.data(), B2[j], F, herm(B2[j]));
   }
   for (int j = UnitCellSize-1; j >= 0; --j)
   {
      F = operator_prod(Ham.data(), B1[j], F, herm(B1[j]));
   }

   // initial set of singular values
   std::deque<MatrixOperator> CenterQ;
   for (int i = 0; i < UnitCellSize; ++i)
   {
      CenterQ.push_back(MakeRandomMatrixOperator(A2[i].Basis2(), B1[UnitCellSize-i-1].Basis1()));
   }

      CHECK_EQUAL(A2.front().Basis2(), CenterQ.front().Basis1());
      CHECK_EQUAL(CenterQ.front().Basis2(), B1.back().Basis1());
      CHECK_EQUAL(A2.back().Basis2(), CenterQ.back().Basis1());
      CHECK_EQUAL(CenterQ.back().Basis2(), B1.front().Basis1());

   Center = MakeRandomMatrixOperator(A2.back().Basis2(), B1.front().Basis1());
   CHECK_EQUAL(Center.Basis1(), Center.Basis2());

   int Iterations = NumIter;
   Center = 0.5 * (Center + adjoint(Center));
   double Energy = Lanczos(Center, 
			   SuperblockMultiply(E, F),
			   Iterations);
   Center = 0.5 * (Center + adjoint(Center));
   TRACE(Energy);

   MPStateComponent A = A2.back(), B = B1.front();

   std::vector<MatrixOperator> SaveTrunc;
   for (int i = 0; i < 400; ++i)
   {
      // truncate
      MatrixOperator RhoL = scalar_prod(Center, herm(Center));
      DensityMatrix<MatrixOperator> DML(RhoL);
      TruncationInfo Info;
      MatrixOperator TruncL = DML.ConstructTruncator(DML.begin(), 
						     TruncateFixTruncationErrorAbsolute(DML.begin(),
											DML.end(),
											SInfo,
											Info));

      // randomize the signs of TruncL
      for (int j = 0; j < TruncL.Basis1().size(); ++j)
      {
	 for (int k = 0; k < TruncL.Basis1().dim(j); ++k)
	 {
	    int sign = ((rand() / 1000) % 2) * 2 - 1;
	    for (int l = 0; l < TruncL.Basis2().size(); ++l)
	    {
	       if (iterate_at(TruncL.data(), j,l))
	       {
		  for (int m = 0; m < TruncL.Basis2().dim(l); ++m)
		  {
		     TruncL(j,l)(k,m) *= sign;
		  }
	       }
	    }
	 }
      }

      TRACE(ExtractDiagonal(scalar_prod(TruncL, herm(TruncL))));

      Center = TruncL * Center;
      A2.back() = prod(A2.back(), herm(TruncL));
      E = triple_prod(TruncL, E, herm(TruncL));

      SaveTrunc.push_back(TruncL);

#if 0      
      MatrixOperator RhoR = scalar_prod(herm(Center), Center);
      DensityMatrix<MatrixOperator> DMR(RhoR);
      TruncationInfo InfoR;
      MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(), 
						     TruncateFixTruncationErrorAbsolute(DMR.begin(),
											DMR.end(),
											SInfo,
											InfoR));
#else
      MatrixOperator TruncR = TruncL;
#endif

      Center = Center * herm(TruncR);
      TRACE(ExtractDiagonal(Center));
      B1.front() = prod(TruncR, B1.front());
      F = triple_prod(TruncR, F, herm(TruncR));

#if 0
      // Make the singular values positive
      MatrixOperator U, D, Vh;
      MatrixOperator CenterCheck = operator_prod(E, Center, herm(F));
      Tensor::SingularValueDecomposition(Center, U, D, Vh);
      Center = D;
      A2.back() = prod(A2.back(), U);
      B1.front() = prod(Vh, B1.front());
      E = triple_prod(herm(U), E, U);
      F = triple_prod(Vh, F, herm(Vh));
      MatrixOperator Cx = U * operator_prod(E, Center, herm(F)) * herm(Vh);
      TRACE(norm_frob(CenterCheck - Cx));
#endif

      double BondE = inner_prod(Center, triple_prod(E[1], Center, herm(F[1]))).real();
      TRACE(BondE);

      //      if (i%2)
      //	 TRACE(ExtractDiagonal(Center));

      CenterQ.back() = Center;

      // The PWFRG transformation
      MatrixOperator pE = MatchBasisReverse(A2.front().Basis1(), A1.front().Basis1());
      for (int i = 0; i < UnitCellSize; ++i)
      {
	 TRACE(SingularValues(pE));
	 pE = operator_prod(herm(A2[i]), pE, A1[i]);
      }
      MPStateComponent A = prod(pE, A2.front());

#if 0
      MatrixOperator pF = MatchBasis(B1.back().Basis2(), B2.back().Basis2());
      for (int i = UnitCellSize-1; i >= 0; --i)
      {
	 pF = operator_prod(B1[i], pF, herm(B2[i]));
      }
#else
      MatrixOperator pF = pE;
#endif
      MPStateComponent B = prod(B1.back(), herm(pF));

      //      TRACE(pE.Basis1())(pE.Basis2())(pF.Basis1())(pF.Basis2());
      TRACE(SingularValues(pE));

      if (min(SingularValues(pE)) < 0.9)
      {
	 TRACE("Singular value failure")(min(SingularValues(pE)));
      }

      // new singular values
      CenterQ.push_back(Center);
      Center = CenterQ.front();
      CenterQ.pop_front();

      CHECK_EQUAL(Center.Basis1(), Center.Basis2());
      TRACE(norm_frob(Center - adjoint(Center)));

      //      TRACE(Center);

      CHECK_EQUAL(A.Basis2(), Center.Basis1());
      CHECK_EQUAL(Center.Basis2(), B.Basis1());

#if 0
      {
	 MatrixOperator U = TruncateBasis2(A);
	 Center = U * Center * herm(U);
      }
      //      {
      //	 MatrixOperator U = TruncateBasis1(B);
      //	 Center = Center * U;
      //      }
      if (norm_frob(Center) > 1E-10)
	 TRACE(ExtractDiagonal(Center));
#endif

      //TRACE(ExtractDiagonal(Center));

      // save R_new which becomes the updated R_old
      // note that we are not necessarily orthogonal at this point
      A1.pop_front();
      A1.push_back(A2.front());
      A2.pop_front();
      A2.push_back(A);

      B2.pop_back();
      B2.push_front(B1.back());
      B1.pop_back();
      B1.push_front(B);

      // expand basis
      CHECK_EQUAL(A.Basis2(), B.Basis1());
      {
	 MatrixOperator U = ExpandBasis2(A);
	 Center = U * Center;
      }
      {
      	 MatrixOperator U = ExpandBasis1(B);
      	 Center = Center * U;
      }

      CHECK_EQUAL(Center.Basis1(), Center.Basis2());

      // hamiltonian
      E = operator_prod(herm(Ham.data()), herm(A), E, A);
      F = operator_prod(Ham.data(), B, F, herm(B));

      //for (int i = 0; i < UnitCellSize; ++i)
      //{
      //	 TRACE(i)(A2[i].Basis2())(CenterQ[i].Basis1());
      //}

      //CHECK_EQUAL(A2.front().Basis2(), CenterQ.front().Basis1());
      //CHECK_EQUAL(CenterQ.front().Basis2(), B1.back().Basis1());
      //CHECK_EQUAL(A2.back().Basis2(), CenterQ.back().Basis1());
      //CHECK_EQUAL(CenterQ.back().Basis2(), B1.front().Basis1());

      if (norm_frob(Center) < 1E-10)
      {
	 TRACE(norm_frob(Center));
	 Center = MakeRandomMatrixOperator(Center.Basis1(), Center.Basis2());
      }

      Center = 0.5 * (Center + adjoint(Center));
      MatrixOperator CenterSave = Center;
      Iterations = NumIter;
      Energy = Lanczos(Center, 
		       SuperblockMultiply(E, F),
		       Iterations);
      Center = 0.5 * (Center + adjoint(Center));
      TRACE(Energy)(inner_prod(Center, CenterSave));

      //      TRACE(norm_frob(Center))(norm_frob(CenterSave));
      //      TRACE(Center)(CenterSave);
      //      TRACE(EigenvaluesHermitian(scalar_prod(Center, herm(Center))))
      //	 (EigenvaluesHermitian(scalar_prod(CenterSave, herm(CenterSave))));
      //      TRACE(Center.Basis1());
      //      TRACE(CenterSave.Basis1());

      // Adjustment for the energy per site
      E[0] -= (Energy * 0.5) * MatrixOperator::make_identity(Center.Basis1());
      F[2] -= (Energy * 0.5) * MatrixOperator::make_identity(Center.Basis2());

      // bond energy
      double BondEnergy = inner_prod(Center, triple_prod(E[1], Center, herm(F[1]))).real();
      TRACE(BondEnergy);

      // over-relaxation step (not effective)
      //if (inner_prod(Center, CenterSave).real() < 0)
      //   Center *= -1.0;
      //MatrixOperator Delta = Center - CenterSave;
      //Center = CenterSave + 0.9 * Delta;

   }
}
