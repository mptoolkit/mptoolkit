// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-pwfrg.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

int main(int argc, char** argv)
{
   int NumIter = 40;

   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(1);

   double J1 = 1;
   double J2 = 0.5;
   int UnitCellSize = 1;

   // initial UnitCellSize-site block
   int MinStates = 1;
   int MaxStates = 10;
   double TruncCutoff = 0;
   double EigenCutoff = -1;

   StatesInfo SInfo;
   SInfo.MinStates = MinStates;
   SInfo.MaxStates = MaxStates;
   SInfo.TruncationCutoff = TruncCutoff;
   SInfo.EigenvalueCutoff = EigenCutoff;
   std::cout << SInfo << '\n';

   MpOpTriangular Ham = ZigZagChain(-sqrt(3.0)*Site["S"], Site["S"], J1, J2);

   MpOpTriangular BoundaryHam = ZigZagChain(-sqrt(3.0)*Boundary["S"], Boundary["S"], J1, J2);

   //   MpOpTriangular Ham = TriangularTwoSite(-sqrt(3.0)*Site["S"],
   //                                     Site["S"],
   //                                     Site["I"].TransformsAs());

   MPStateComponent E = Initial_E(BoundaryHam);
   MPStateComponent F = Initial_F(BoundaryHam);

   //   MatrixOperator OldCenter = Center;
   typedef std::deque<MPStateComponent> MatListType;

   MatrixOperator Center = MatrixOperator::make_identity(E.Basis1());
   MPStateComponent A1 = ConstructFromLeftBasis(Boundary.Basis2().Basis(), Center.Basis1());
   MPStateComponent B2 = ConstructFromRightBasis(Boundary.Basis1().Basis(), Center.Basis2());
   E = operator_prod(herm(BoundaryHam.data()), herm(A1), E, A1);
   F = operator_prod(BoundaryHam.data(), B2, F, herm(B2));

   Center = MakeRandomMatrixOperator(A1.Basis2(), B2.Basis1());

   MatListType A, B;

   A.push_back(ConstructFromLeftBasis(Site.Basis2().Basis(), Center.Basis1()));
   for (int j = 0; j < UnitCellSize-1; ++j)
   {
      A.push_back(ConstructFromLeftBasis(Site.Basis2().Basis(), A.back().Basis2()));
   }

   B.push_front(ConstructFromRightBasis(Site.Basis1().Basis(), Center.Basis2()));
   for (int j = 0; j < UnitCellSize-1; ++j)
   {
      B.push_front(ConstructFromRightBasis(Site.Basis1().Basis(), B.front().Basis1()));
   }

   for (int j = 0; j < UnitCellSize; ++j)
   {
      E = operator_prod(herm(Ham.data()), herm(A[j]), E, A[j]);
   }

   for (int j = UnitCellSize-1; j >= 0; --j)
   {
      F = operator_prod(Ham.data(), B[j], F, herm(B[j]));
   }

   // initial set of singular values
   std::deque<MatrixOperator> CenterQ;
   for (int i = 0; i < UnitCellSize; ++i)
   {
      CenterQ.push_back(MakeRandomMatrixOperator(A[i].Basis2(), B[UnitCellSize-i-1].Basis1()));
   }

   MPStateComponent R = A.back(), S = B.front();

   Center = MakeRandomMatrixOperator(R.Basis2(), S.Basis1());
   CHECK_EQUAL(Center.Basis1(), Center.Basis2());

   int Iterations = NumIter;
   //Center = 0.5 * (Center + adjoint(Center));
   double Energy = Lanczos(Center,
                           SuperblockMultiply(E, F),
                           Iterations);
   //Center = 0.5 * (Center + adjoint(Center));
   TRACE(Energy);

   std::vector<MatrixOperator> SaveTrunc;
   for (int i = 0; i < 10000; ++i)
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
      for (unsigned j = 0; j < TruncL.Basis1().size(); ++j)
      {
         for (int k = 0; k < TruncL.Basis1().dim(j); ++k)
         {
            int sign = ((rand() / 1000) % 2) * 2 - 1;
            for (unsigned l = 0; l < TruncL.Basis2().size(); ++l)
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

      Center = TruncL * Center;
      A.back() = prod(A.back(), herm(TruncL));
      E = triple_prod(TruncL, E, herm(TruncL));

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
#if 0
      // randomize the signs of TruncL
      for (unsigned j = 0; j < TruncL.Basis1().size(); ++j)
      {
         for (int k = 0; k < TruncL.Basis1().dim(j); ++k)
         {
            int sign = ((rand() / 1000) % 2) * 2 - 1;
            for (unsigned l = 0; l < TruncL.Basis2().size(); ++l)
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
#endif
      MatrixOperator TruncR = TruncL;
#endif

      Center = Center * herm(TruncR);
      TRACE(ExtractDiagonal(Center));
      TRACE(SingularValues(Center));
      B.front() = prod(TruncR, B.front());
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
      //         TRACE(ExtractDiagonal(Center));

      CenterQ.back() = Center;

      if (i%1000 == 0)
      {
         R = A.back();
         S = B.front();
      }

      // The PWFRG transformation
      TRACE(SingularValues(scalar_prod(herm(A.back()), R)));
      MatrixOperator UE = scalar_prod(herm(A.back()), R) * MatchBasisReverse(R.Basis2(), A.front().Basis1());
      //TRACE(UE.Basis1())(UE.Basis2());

      R = prod(UE, A.front());
      //      R *= R.Basis1().total_degree() / norm_frob_sq(R);
      A.push_back(R);
      A.pop_front();

      S = prod(B.back(), herm(UE));
      //S *= S.Basis2().total_degree() / norm_frob_sq(S);
      B.push_front(S);
      B.pop_back();

      // new singular values
      CenterQ.push_back(Center);
      Center = CenterQ.front();
      CenterQ.pop_front();

      TRACE(SingularValues(Center));

      //      TRACE(Center);

      CHECK_EQUAL(A.back().Basis2(), Center.Basis1());
      CHECK_EQUAL(Center.Basis2(), B.front().Basis1());
      CHECK_EQUAL(A.back().Basis2(), B.front().Basis1());
      {
         MatrixOperator U = ExpandBasis2(A.back());
         Center = U * Center;
      }
      {
         MatrixOperator U = ExpandBasis1(B.front());
         Center = Center * U;
      }

      TRACE(norm_frob(Center));
      CHECK_EQUAL(Center.Basis1(), Center.Basis2());

      TRACE(norm_frob(R))(R.Basis1().total_degree())(R.Basis2().total_degree());
      TRACE(norm_frob(S))(S.Basis1().total_degree())(S.Basis2().total_degree());

      // normalize
      Center *= 1.0 / norm_frob(Center);

      // hamiltonian
      E = operator_prod(herm(Ham.data()), herm(A.back()), E, A.back());
      F = operator_prod(Ham.data(), B.front(), F, herm(B.front()));

      //for (int i = 0; i < UnitCellSize; ++i)
      //{
      //         TRACE(i)(A2[i].Basis2())(CenterQ[i].Basis1());
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

      //Center = 0.5 * (Center + adjoint(Center));
      MatrixOperator CenterSave = Center;
      Iterations = NumIter;
      Energy = Lanczos(Center,
                       SuperblockMultiply(E, F),
                       Iterations);
      //Center = 0.5 * (Center + adjoint(Center));
      TRACE(Energy)(inner_prod(Center, CenterSave));

      //      TRACE(norm_frob(Center))(norm_frob(CenterSave));
      //      TRACE(Center)(CenterSave);
      //      TRACE(EigenvaluesHermitian(scalar_prod(Center, herm(Center))))
      //         (EigenvaluesHermitian(scalar_prod(CenterSave, herm(CenterSave))));
      //      TRACE(Center.Basis1());
      //      TRACE(CenterSave.Basis1());

      // Adjustment for the energy per site
      E[0] -= (Energy * 0.5) * MatrixOperator::make_identity(Center.Basis1());
      F[3] -= (Energy * 0.5) * MatrixOperator::make_identity(Center.Basis2());

      // bond energy
      double BondEnergy = inner_prod(Center, triple_prod(E[1], Center, herm(F[1]))).real()
         + inner_prod(Center, triple_prod(E[2], Center, herm(F[2]))).real();
      TRACE(BondEnergy);

      // over-relaxation step (not effective)
      //if (inner_prod(Center, CenterSave).real() < 0)
      //   Center *= -1.0;
      //MatrixOperator Delta = Center - CenterSave;
      //Center = CenterSave + 0.9 * Delta;

   }
}
