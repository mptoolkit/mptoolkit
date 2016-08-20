// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-idmrg3.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "matrixproduct/infinitewavefunction.h"
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
#include "models/spin-u1.h"
#include "models/spin.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"

#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

struct SuperblockMultiplyX
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiplyX(MPStateComponent const& Left_,
                      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }

   MPStateComponent Left, Right;
};

inline
SuperblockMultiplyX::SuperblockMultiplyX(MPStateComponent const& Left_,
                                       MPStateComponent const& Right_)
   : Left(Left_), Right(Right_)
{
}


struct SuperblockMultiply
{
   typedef MPStateComponent result_type;
   typedef MPStateComponent value_type;
   typedef MPStateComponent argument_type;

   SuperblockMultiply(MPOpComponent const& H_,
                      MPStateComponent const& Left_,
                      MPStateComponent const& Right_);

   MPStateComponent operator()(MPStateComponent const& Psi) const
   {
      return local_operator_prod(H, Left, Psi, herm(Right));
   }

   MPOpComponent H;
   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(MPOpComponent const& H_,
                                       MPStateComponent const& Left_,
                                       MPStateComponent const& Right_)
   : H(H_), Left(Left_), Right(Right_)
{
}

// does a DMRG sweep moving right.  In input, the IncomingHam is defined on Basis1() of
// the incoming wavefunction.  The RightBlockHam is defined on Psi.front().Basis2().
// On exit, RightBlockHam is replaced by a container of the left block hamiltonians.

void DoDMRGSweepRight(LinearWavefunction& Psi,
                 SimpleMPOperator const& Ham,
                 std::deque<MPStateComponent>& LeftBlockHam,
                 std::deque<MPStateComponent>& RightBlockHam,
                 StatesInfo const& SInfo)
{
   LinearWavefunction::iterator I = Psi.begin();
   LinearWavefunction::iterator J = I; ++J;

   SimpleMPOperator::const_iterator H1 = Ham.begin();
   SimpleMPOperator::const_iterator H2 = H1; ++H2;

   MPStateComponent R = *I;
   while (J != Psi.end())
   {
      MPStateComponent A = local_tensor_prod(R, *J);
      MPOpComponent H = local_tensor_prod(*H1, *H2);

      // apply the solver
      int Iterations = 10;
      double Energy = Lanczos(A, SuperblockMultiply(H, LeftBlockHam.back(), RightBlockHam.front()),
                              Iterations);

      // truncate
      AMatSVD SL(A, Tensor::ProductBasis<BasisList, BasisList>(R.LocalBasis(), J->LocalBasis()));
      TruncationInfo Info;
      AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                                  SInfo, Info);

      std::cout << "Energy=" << Energy
                << " States=" << Info.KeptStates()
                << " TruncError=" << Info.TruncationError()
                << " Entropy=" << Info.KeptEntropy() << '\n';

      MatrixOperator C;
      SL.ConstructMatrices(SL.begin(), Cutoff, *I, C, R);
      R = prod(C, R);

      // new hamiltonian operators
      LeftBlockHam.push_back(operator_prod(*H1, *I, LeftBlockHam.back(), herm(*I)));
      RightBlockHam.pop_front();

      I=J;
      ++J;
      H1 = H2;
      ++H2;
   }
}

// does a DMRG sweep moving right.  In input, the IncomingHam is defined on Basis1() of
// the incoming wavefunction.  The RightBlockHam is defined on Psi.front().Basis2().
// On exit, RightBlockHam is replaced by a container of the left block hamiltonians.

void DoDMRGSweepLeft(LinearWavefunction& Psi,
                     SimpleMPOperator const& Ham,
                     std::deque<MPStateComponent>& LeftBlockHam,
                     std::deque<MPStateComponent>& RightBlockHam,
                     StatesInfo const& SInfo)
{
   LinearWavefunction::iterator J = Psi.end();
   LinearWavefunction::iterator I = J; --I;

   SimpleMPOperator::const_iterator H2 = Ham.end();
   SimpleMPOperator::const_iterator H1 = H2; --H1;

   MPStateComponent R = *I;
   while (I != Psi.begin())
   {
      J = I; --I;
      H2 = H1; --H1;

      MPStateComponent A = local_tensor_prod(*I, R);
      MPOpComponent H = local_tensor_prod(*H1, *H2);

      // apply the solver
      int Iterations = 10;
      double Energy = Lanczos(A, SuperblockMultiply(H, LeftBlockHam.back(), RightBlockHam.front()),
                              Iterations);

      // truncate
      AMatSVD SL(A, Tensor::ProductBasis<BasisList, BasisList>(I->LocalBasis(), R.LocalBasis()));
      TruncationInfo Info;
      AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                                  SInfo, Info);

      std::cout << "Energy=" << Energy
                << " States=" << Info.KeptStates()
                << " TruncError=" << Info.TruncationError()
                << " Entropy=" << Info.KeptEntropy() << '\n';

      MatrixOperator C;
      SL.ConstructMatrices(SL.begin(), Cutoff, R, C, *J);
      R = prod(R, C);

      // new hamiltonian operators
      RightBlockHam.push_front(operator_prod(herm(*H2), herm(*J), RightBlockHam.front(), *J));
      LeftBlockHam.pop_back();
   }
}

void DoStep(LinearWavefunction& Left, MatrixOperator& C, LinearWavefunction& Right,
                 SimpleMPOperator const& LeftHam, SimpleMPOperator const& RightHam,
                 std::deque<MPStateComponent>& LeftBlockHam,
                 std::deque<MPStateComponent>& RightBlockHam,
                 StatesInfo SInfo)
{
   MPStateComponent A = local_tensor_prod(prod(Left.get_back(), C), Right.get_front());
   MPOpComponent H = local_tensor_prod(LeftHam.back(), RightHam.front());

   int Iterations = 10;

   double Energy = Lanczos(A, SuperblockMultiply(H, LeftBlockHam.back(), RightBlockHam.front()),
                           Iterations);

   // truncate
   AMatSVD SL(A, Tensor::ProductBasis<BasisList, BasisList>(Left.get_back().LocalBasis(),
                                                            Right.get_front().LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                               SInfo, Info);

   std::cout << "Energy=" << Energy
             << " States=" << Info.KeptStates()
             << " TruncError=" << Info.TruncationError()
             << " Entropy=" << Info.KeptEntropy() << '\n';

   MPStateComponent L, R;
   SL.ConstructMatrices(SL.begin(), Cutoff, L, C, R);

   Left.set_back(L);
   Right.set_front(R);
}

void DoIteration(LinearWavefunction& Left, MatrixOperator& C, LinearWavefunction& Right,
                 SimpleMPOperator const& LeftHam, SimpleMPOperator const& RightHam,
                 std::deque<MPStateComponent>& LeftBlockHam,
                 std::deque<MPStateComponent>& RightBlockHam,
                 MatrixOperator& C_old,
                 StatesInfo SInfo)
{
   // These two could be done in parallel
   std::deque<MPStateComponent> NewRightBlockHam;
   NewRightBlockHam.push_front(operator_prod(RightHam.front(),
                                             Right.get_front(),
                                             RightBlockHam.front(),
                                             herm(Right.get_front())));
   std::deque<MPStateComponent> NewLeftBlockHam;
   NewLeftBlockHam.push_back(operator_prod(herm(LeftHam.back()),
                                           herm(Left.get_back()),
                                           LeftBlockHam.back(),
                                           Left.get_back()));

   Left.set_back(prod(Left.get_back(), C));
   //   DoDMRGSweepLeft(Left, LeftHam, LeftBlockHam, NewRightBlockHam, SInfo);

   Right.set_front(prod(C, Right.get_front()));
   //   DoDMRGSweepRight(Right, RightHam, NewLeftBlockHam, RightBlockHam, SInfo);

   MatrixOperator C_new = InvertDiagonal(C_old, 1E-7);
   C_old = C;
   C = C_new;
   DoStep(Right, C, Left, RightHam, LeftHam, NewLeftBlockHam, NewRightBlockHam, SInfo);
   //   std::swap(Left, Right);
   LeftBlockHam = NewLeftBlockHam;
   RightBlockHam = NewRightBlockHam;
}

void DoSweep(LinearWavefunction& Left, MatrixOperator& C, LinearWavefunction& Right,
             SimpleMPOperator const& LeftHam, SimpleMPOperator const& RightHam,
             std::deque<MPStateComponent>& LeftBlockHam,
             std::deque<MPStateComponent>& RightBlockHam,
             MatrixOperator& C_old,
             StatesInfo SInfo)
{
   DoIteration(Left, C, Right, LeftHam, RightHam, LeftBlockHam, RightBlockHam, C_old, SInfo);
   DoIteration(Right, C, Left, RightHam, LeftHam, LeftBlockHam, RightBlockHam, C_old, SInfo);
}

int main(int argc, char** argv)
{
   std::string OutName = "test.out";

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(OutName, 1, PageSize, CacheSize);

   int MainIter = argc > 1 ? boost::lexical_cast<int>(argv[1]) : 100;

   //   mp_pheap::InitializeTempPHeap(false);

   int NumIter = 10;

   double J1 = 1;
   double J2 = 0.5;

   int MinStates = 1;
   int MaxStates = 40;
   double TruncCutoff = 0;
   double EigenCutoff = 1e-12;

   std::cout.precision(12);

   StatesInfo SInfo;
   SInfo.MinStates = MinStates;
   SInfo.MaxStates = MaxStates;
   SInfo.TruncationCutoff = TruncCutoff;
   SInfo.EigenvalueCutoff = EigenCutoff;
   std::cout << SInfo << '\n';

   // Ising model
#if 1
   double Lambda = 1.0;
   SiteBlock Boundary = CreateSpinSite(0.5);
   SiteBlock Site = CreateSpinSite(0.5);

   MpOpTriangular Ham = 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
      + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);

   MpOpTriangular BoundaryHam = Ham;
#endif

   // SU(2) Heisenberg model
#if 0
   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(0.5);

   MpOpTriangular Ham = TriangularTwoSite(-sqrt(3.0)*Site["S"],
                                          Site["S"],
                                          Site["I"].TransformsAs());

   MpOpTriangular BoundaryHam = TriangularTwoSite(-sqrt(3.0)*Boundary["S"],
                                                  Boundary["S"],
                                                  Boundary["I"].TransformsAs());
#endif

   // Heisenberg model with no symmetries
#if 0
   SiteBlock Boundary = CreateSpinSite(0.5);
   SiteBlock Site = CreateSpinSite(0.5);

   MpOpTriangular Ham = 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"]) +
                               TriangularTwoSite(Site["Sm"], Site["Sp"]))
      + TriangularTwoSite(Site["Sz"], Site["Sz"]);

   MpOpTriangular BoundaryHam = 0.5 * (TriangularTwoSite(Boundary["Sp"], Boundary["Sm"]) +
                                       TriangularTwoSite(Boundary["Sm"], Boundary["Sp"]))
      + TriangularTwoSite(Boundary["Sz"], Boundary["Sz"]);

   TRACE(Ham.data());
#endif

   //MpOpTriangular Ham = ZigZagChain(-sqrt(3.0)*Site["S"], Site["S"], J1, J2);

   //MpOpTriangular BoundaryHam = ZigZagChain(-sqrt(3.0)*Boundary["S"], Boundary["S"], J1, J2);

   MPStateComponent E = Initial_E(BoundaryHam);
   MPStateComponent F = Initial_F(BoundaryHam);

   // initial 2-site block
   MatrixOperator Center = MatrixOperator::make_identity(E.Basis1());
   MPStateComponent A1 = ConstructFromLeftBasis(Boundary.Basis2().Basis(), Center.Basis1());
   MPStateComponent B2 = ConstructFromRightBasis(Boundary.Basis1().Basis(), Center.Basis2());

   Center = MakeRandomMatrixOperator(A1.Basis2(), B2.Basis1());
   MatrixOperator OldCenter = MakeRandomMatrixOperator(B2.Basis2(), A1.Basis1());
   OldCenter *= 1.0 / norm_frob(OldCenter);

   // set up

   LinearWavefunction Left, Right;
   Left.push_back(A1);
   Right.push_front(B2);

   SimpleMPOperator LeftHam, RightHam;
   LeftHam.push_back(Ham.data());
   RightHam.push_front(Ham.data());

   MatrixOperator C_LR = Center;
   MatrixOperator C_RL = OldCenter;

   std::deque<MPStateComponent> LeftBlockHam, RightBlockHam;
   LeftBlockHam.push_back(E);
   RightBlockHam.push_front(F);

   DoStep(Left, C_LR, Right, LeftHam, RightHam, LeftBlockHam, RightBlockHam, SInfo);

   for (int i = 0; i < MainIter; ++i)
   {
      DoSweep(Left, C_LR, Right, LeftHam, RightHam, LeftBlockHam, RightBlockHam, C_RL, SInfo);

      if (i < 100 || (i%100 == 0))
      {
         std::string FileName = "dmrg.psi." + boost::lexical_cast<std::string>(i);
         // Convert the wavefunction into an InfiniteWavefunction
         //C_RL *= 1.0 / norm_frob(C_RL);  // normalize after the last truncation
         InfiniteWavefunction Psi;
         Psi.C_old = C_RL;
         for (LinearWavefunction::const_iterator I = Left.begin(); I != Left.end(); ++I)
         {
            Psi.Psi.push_back(*I);
         }

         MatrixOperator U = C_LR;
         for (LinearWavefunction::const_iterator I = Right.begin(); I != Right.end(); ++I)
         {
            MPStateComponent x = prod(U, *I);
            U = TruncateBasis2(x);
            Psi.Psi.push_back(x);
         }
         Psi.C_right = U;
         Psi.QShift = QuantumNumbers::QuantumNumber(U.GetSymmetryList());

         pvalue_ptr<InfiniteWavefunction> Out = new InfiniteWavefunction(Psi);
         pheap::ExportHeap(FileName, Out);
      }
   }
}
