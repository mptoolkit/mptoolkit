// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/prodoptimizer.cpp
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

#include "prodoptimizer.h"

template <typename FwdIter>
double Entropy(FwdIter first, FwdIter last)
{
   double E = 0;
   while (first != last)
   {
      double x = first->Eigenvalue;
      if (x > 0)
         E -= x * log(x);

      ++first;
   }
   return E;
}

template <typename FwdIter>
double TruncError(FwdIter first, FwdIter last)
{
   double E = 0;
   while (first != last)
   {
      double x = first->Eigenvalue;
      if (x > 0)
         E += x;

      ++first;
   }
   return 1.0 - E;
}

ProductOptimizer::ProductOptimizer(CenterWavefunction const& Psi_, SplitOperator const& A_,
                                   CenterWavefunction const& Rhs_)
   : Psi(Psi_), A(1, A_), Rhs(1, Rhs_), Psi_A_Rhs(1), Ident(Psi_.GetSymmetryList())
{
   // construct the OpMatrix elements for the right hand side

   this->Initialize();
}

ProductOptimizer::ProductOptimizer(CenterWavefunction const& Psi_, std::vector<SplitOperator> const& A_,
                                   std::vector<CenterWavefunction> const& Rhs_)
   : Psi(Psi_), A(A_), Rhs(Rhs_), Psi_A_Rhs(Rhs_.size()), Ident(Psi_.GetSymmetryList())
{
   CHECK_EQUAL(A.size(), Rhs.size());
   // construct the OpMatrix elements for the right hand side

   this->Initialize();
}

void
ProductOptimizer::Initialize()
{
   typedef MPStateComponent Component;

   // Shift everything to normal form
   while (Psi.LeftSize() > 1) Psi.RotateLeft();

   // The initial matrix elements represent the vacuum basis on the right hand side
   BasisList Vacuum = make_vacuum_basis(Psi.GetSymmetryList());

   MatrixOperator VacIdent = MatrixOperator(VectorBasis(Vacuum), VectorBasis(Vacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   for (unsigned i = 0; i < A.size(); ++i)
   {
      while (A[i].LeftSize() > 1) A[i].RotateLeft();
      while (Rhs[i].LeftSize() > 1) Rhs[i].RotateLeft();


      Psi_A_Rhs[i].PushRight(make_vacuum_state(Psi.GetSymmetryList()));

      // apply the Right matrices
      for (int Loc = 0; Loc < Psi.RightSize(); ++Loc)
      {
         Psi_A_Rhs[i].PushRight(operator_prod(A[i].LookupRight(Loc),
                                              Psi.LookupRight(Loc),
                                              Psi_A_Rhs[i].Right(),
                                              herm(Rhs[i].LookupRight(Loc))));
      }

      //Now start from the left hand side and work inwards.
      // Fix 2007-01-24: make sure we get the right quantum numbers for the left-handed vacuum
      BasisList LeftVacuum(Psi.GetSymmetryList());
      CHECK_EQUAL(Psi.LookupLeft(0).Basis1().size(), 1);
      LeftVacuum.push_back(Psi.LookupLeft(0).Basis1()[0]);

      BasisList LeftVacRhs(Psi.GetSymmetryList());
      CHECK_EQUAL(Rhs[i].LookupLeft(0).Basis1().size(), 1);
      LeftVacRhs.push_back(Rhs[i].LookupLeft(0).Basis1()[0]);

      VacIdent = MatrixOperator(VectorBasis(LeftVacuum), VectorBasis(LeftVacRhs), A[i].TransformsAs());
      VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

      MPStateComponent LVac(A[i].LookupLeft(0).Basis1(), VectorBasis(LeftVacuum), VectorBasis(LeftVacRhs));
      LVac[0] = VacIdent;

      //   Psi_A_Rhs.PushLeft(make_vacuum_state(Psi.LookupLeft(0).Basis1()[0]));
      Psi_A_Rhs[i].PushLeft(LVac);

      for (int Loc = 0; Loc < Psi.LeftSize(); ++Loc)
      {
         Psi_A_Rhs[i].PushLeft(operator_prod(herm(A[i].LookupLeft(Loc)),
                                             herm(Psi.LookupLeft(Loc)),
                                             Psi_A_Rhs[i].Left(),
                                             Rhs[i].LookupLeft(Loc)));
      }
   }
}

double ProductOptimizer::Energy()
{
   //   DEBUG_TRACE(OpMatrices.Left())(OpMatrices.Right());
   //   DEBUG_TRACE(scalar_prod(Psi.Center(), Psi.Center()));
   //   DEBUG_TRACE(scalar_prod(PsiP, PsiP));

   double Res = 0;
   for (unsigned i = 0; i < A.size(); ++i)
   {
      Res += inner_prod(operator_prod(A[i].Center(),
                                      Psi_A_Rhs[i].Left(),
                                      Rhs[i].Center(),
                                      herm(Psi_A_Rhs[i].Right())),
                        Psi.Center()).real();
   }
   return Res;
}

double ProductOptimizer::Solve()
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   Psi.Center() = operator_prod(A[0].Center(),
                                Psi_A_Rhs[0].Left(),
                                Rhs[0].Center(),
                                herm(Psi_A_Rhs[0].Right()));

   for (unsigned i = 1; i < A.size(); ++i)
   {
      Psi.Center() += operator_prod(A[i].Center(),
                                    Psi_A_Rhs[i].Left(),
                                    Rhs[i].Center(),
                                    herm(Psi_A_Rhs[i].Right()));
   }

   return norm_frob(Psi.Center());
}

void ProductOptimizer::ExpandLeft()
{
   ExpandBasis2(Psi.Left());

   // matrix elements
   for (unsigned i = 0; i < A.size(); ++i)
   {
      Psi_A_Rhs[i].PopLeft();
      Psi_A_Rhs[i].PushLeft(operator_prod(herm(A[i].Left()),
                                       herm(Psi.Left()),
                                       Psi_A_Rhs[i].Left(),
                                       Rhs[i].Left()));
   }
}

void ProductOptimizer::ExpandRight()
{
   ExpandBasis1(Psi.Right());

   // matrix elements
   for (unsigned i = 0; i < A.size(); ++i)
   {
      Psi_A_Rhs[i].PopRight();
      Psi_A_Rhs[i].PushRight(operator_prod(A[i].Right(),
                                        Psi.Right(),
                                        Psi_A_Rhs[i].Right(),
                                        herm(Rhs[i].Right())));
   }
}

void ProductOptimizer::ShiftRightAndExpand()
{
   Psi.PushLeft(prod(Psi.Center(), Psi.Right()));
   Psi.PopRight();
   ExpandBasis2(Psi.Left());  // dont update yprime.Center() - do later


   // Rotate the wavefunctions to the right.  A and x are fixed and do not
   // need the basis expanded
   for (unsigned i = 0; i < A.size(); ++i)
   {
      A[i].RotateRight();
      Rhs[i].RotateRight();
      Psi_A_Rhs[i].PushLeft(operator_prod(herm(A[i].Left()), herm(Psi.Left()), Psi_A_Rhs[i].Left(), Rhs[i].Left()));
      Psi_A_Rhs[i].PopRight();
   }
}

void ProductOptimizer::ShiftLeftAndExpand()
{
   Psi.PushRight(prod(Psi.Left(), Psi.Center()));
   Psi.PopLeft();
   ExpandBasis1(Psi.Right());  // dont update yprime.Center() - do later

   for (unsigned i = 0; i < A.size(); ++i)
   {
      // Rotate the wavefunctions to the left.  A and x are fixed and do not
      // need the basis expanded
      A[i].RotateLeft();
      Rhs[i].RotateLeft();

      // new matrix elements
      Psi_A_Rhs[i].PushRight(operator_prod(A[i].Right(),
                                           Psi.Right(),
                                           Psi_A_Rhs[i].Right(),
                                           herm(Rhs[i].Right())));
      Psi_A_Rhs[i].PopLeft();
   }
}

TruncationInfo ProductOptimizer::TruncateLeft(StatesInfo const& SInfo, double CFactor)
{
   MatrixOperator ydm = scalar_prod(Psi.Center(), herm(Psi.Center()));
   if (CFactor != 0)
   {
      int Count = 0;
      MatrixOperator Rhop;
      for (unsigned i = 0; i < A.size(); ++i)
      {
         MatrixOperator RhoRhs = scalar_prod(Rhs[i].Center(), herm(Rhs[i].Center()));
         for (unsigned j = 0; j < Psi_A_Rhs[i].Left().size(); ++j)
         {
            MatrixOperator Correction = triple_prod(Psi_A_Rhs[i].Left()[j],
                                                    RhoRhs,
                                                    herm(Psi_A_Rhs[i].Left()[j]));
            std::complex<double> Tr = trace(Correction);
            if (LinearAlgebra::norm_2(Tr) > std::numeric_limits<double>::epsilon() * 100)
            {
               Correction *= CFactor / Tr;
               Rhop += Correction;
               Count++;
            }
         }
      }
      ydm += (trace(ydm).real() / Count) * Rhop;
   }

   DensityMatrix<MatrixOperator> DM(ydm);
   //DM.DensityMatrixReport(std::cerr);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                                                                           SInfo,
                                                                                           Info));
   // truncate
   Psi.Center() = prod(U, Psi.Center(), Ident);
   Psi.Left() = prod(Psi.Left(), adjoint(U));
   for (unsigned i = 0; i < A.size(); ++i)
   {
      Psi_A_Rhs[i].Left() = prod(U, Psi_A_Rhs[i].Left());
   }

   return Info;
}

TruncationInfo ProductOptimizer::TruncateRight(StatesInfo const& SInfo, double CFactor)
{
   MatrixOperator ydm = scalar_prod(herm(Psi.Center()), Psi.Center());
   if (CFactor != 0)
   {
      int Count = 0;
      MatrixOperator Rhop;
      for (unsigned i = 0; i < A.size(); ++i)
      {
         MatrixOperator RhoRhs = scalar_prod(herm(Rhs[i].Center()), Rhs[i].Center());
         for (unsigned j = 0; j < Psi_A_Rhs[i].Right().size(); ++j)
         {
            MatrixOperator Correction = triple_prod(Psi_A_Rhs[i].Right()[j],
                                                    RhoRhs,
                                                    herm(Psi_A_Rhs[i].Right()[j]));
            std::complex<double> Tr = trace(Correction);
            if (LinearAlgebra::norm_2(Tr) > std::numeric_limits<double>::epsilon() * 100)
            {
               Correction *= CFactor / Tr;
               Rhop += Correction;
               Count++;
            }
         }
      }
      ydm += (trace(ydm).real() / Count) * Rhop;
   }

   DensityMatrix<MatrixOperator> DM(ydm);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                                                                            SInfo,
                                                                                            Info));
   // truncate
   Psi.Center() = prod(Psi.Center(), adjoint(U), Ident);
   Psi.Right() = prod(U, Psi.Right());
   for (unsigned i = 0; i < A.size(); ++i)
   {
      Psi_A_Rhs[i].Right() = prod(U, Psi_A_Rhs[i].Right());
   }

   return Info;
}
