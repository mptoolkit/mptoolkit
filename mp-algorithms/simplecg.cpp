// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/simplecg.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "simplecg.h"
#include "linearalgebra/eigen.h"
#include "tensor/tensor_eigen.h"
#include "common/messagelogger.h"
#include "tensor/regularize.h"

using MessageLogger::Logger;
using MessageLogger::msg_log;

typedef std::complex<double> complex;

Solver::Solver(CenterWavefunction const& Psi_, SplitOperator const& Op_, CenterWavefunction const& Rhs_,
               double Freq_, double Broad_)
   : x(Psi_), A(Op_), y(Rhs_), Ident(Psi_.GetSymmetryList()), Frequency(Freq_),
     Broadening(Broad_),
     TotalSweepNumber(0),
     TotalSweepRecNumber(0),
     TotalNumIterations(0),
     TotalNumMultiplies(0)
{
   // construct the OpMatrix elements for the right hand side

   typedef MPStateComponent Component;

   // Shift A so the Center matrix lines up with that of x
   while (A.LeftSize() > x.LeftSize()) A.RotateLeft();
   while (A.RightSize() > x.RightSize()) A.RotateRight();
   // Shift y so the Center matrix lines up with that of x
   while (y.LeftSize() > x.LeftSize()) y.RotateLeft();
   while (y.RightSize() > x.RightSize()) y.RotateRight();

   // The initial matrix elements represent the vacuum basis on the right hand side
   BasisList Vacuum = make_vacuum_basis(x.GetSymmetryList());

   MatrixOperator VacIdent = MatrixOperator(VectorBasis(Vacuum), VectorBasis(Vacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   x_A_x.PushRight(make_vacuum_state(x.GetSymmetryList()));
   x_y.PushRight(VacIdent);

   // apply the Right matrices
   for (int Loc = 0; Loc < x.RightSize(); ++Loc)
   {
      x_A_x.PushRight(operator_prod(A.LookupRight(Loc),
                                    x.LookupRight(Loc),
                                    x_A_x.Right(),
                                    herm(x.LookupRight(Loc))));

      x_y.PushRight(operator_prod(x.LookupRight(Loc),
                                  x_y.Right(),
                                  herm(y.LookupRight(Loc))));
   }

   //Now start from the left hand side and work inwards
   BasisList LeftVacuum(x.GetSymmetryList());
   CHECK_EQUAL(x.LookupLeft(0).Basis1().size(), 1);
   LeftVacuum.push_back(x.LookupLeft(0).Basis1()[0]);

   VacIdent = MatrixOperator(VectorBasis(LeftVacuum), VectorBasis(LeftVacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   x_A_x.PushLeft(make_vacuum_state(x.LookupLeft(0).Basis1()[0]));
   x_y.PushLeft(VacIdent);

   for (int Loc = 0; Loc < x.LeftSize(); ++Loc)
   {
      x_A_x.PushLeft(operator_prod(herm(A.LookupLeft(Loc)),
                                   herm(x.LookupLeft(Loc)),
                                   x_A_x.Left(),
                                   x.LookupLeft(Loc)));

      x_y.PushLeft(operator_prod(herm(x.LookupLeft(Loc)),
                                 x_y.Left(),
                                 y.LookupLeft(Loc)));
   }

   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));

   DEBUG_TRACE(A.LeftSize())(x.LeftSize());

   this->DebugCheckBasis();
}

Solver::~Solver()
{
}

CenterWavefunction& Solver::Wavefunction()
{
   x.AttributesMutable()["GreensFunction"] = boost::lexical_cast<std::string>(this->GreensFunction());
   return x;
}

CenterWavefunction const& Solver::Wavefunction() const
{
   x.AttributesMutable()["GreensFunction"] = boost::lexical_cast<std::string>(this->GreensFunction());
   return x;
}

void Solver::CreateLogFiles(std::string const& BasePath, ConfList const& Conf)
{
   bool const ClearLogs = Conf.Get("ClearLogs", true);
   std::ios_base::openmode Mode = ClearLogs ? (std::ios_base::out | std::ios_base::trunc)
      : (std::ios_base::out | std::ios_base::app);

   std::string const EnergyLogFile = BasePath + ".e";
   int const EnergyLogLevel = Conf.Get("EnergyLogLevel", 0);
   Logger("EnergyLog").SetThreshold(EnergyLogLevel);
   if (EnergyLogLevel == 0 && ClearLogs)
      remove(EnergyLogFile.c_str());
   else if (EnergyLogLevel > 0)
   {
      EnergyLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(EnergyLogFile.c_str(), Mode));
      EnergyLog->precision(15);
      (*EnergyLog) << "#TotSweepNum #SweepNum #LSize #RSize "
                   << "#NStates #NMult #Trunc #DOS_norm #DOS_overlap #GF_real #Local_resid\n";
      Logger("EnergyLog").SetStream(*EnergyLog);
   }

   std::string const SweepLogFile = BasePath + ".sweep";
   int const SweepLogLevel = Conf.Get("SweepLogLevel", 0);
   Logger("SweepLog").SetThreshold(SweepLogLevel);
   if (SweepLogLevel == 0 && ClearLogs)
      remove(SweepLogFile.c_str());
   else if (SweepLogLevel > 0)
   {
      SweepLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(SweepLogFile.c_str(), Mode));
      SweepLog->precision(15);
      (*SweepLog) << "#TotSweepNum #SweepNum #NIter #AvNStates #MaxNStates "
                  << "#NMult #Freq #Broad #DOS_norm #DOS_overlap #DOS_func #GF_real #Trunc "
                  << "#Entropy #G_resid #Converged\n";
      Logger("SweepLog").SetStream(*SweepLog);
   }

   std::string const DensityLogFile = BasePath + ".density";
   int const DensityLogLevel = Conf.Get("DensityLogLevel", 0);
   Logger("DensityLog").SetThreshold(DensityLogLevel);
   if (DensityLogLevel == 0 && ClearLogs)
      remove(DensityLogFile.c_str());
   else if (DensityLogLevel > 0)
   {
      DensityLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(DensityLogFile.c_str(), Mode));
      DensityLog->precision(15);
      Logger("DensityLog").SetStream(*DensityLog);
   }

   // CpuLog
   // DiagLog

   this->ReadConfOptions(Conf);
}

void Solver::RestoreLogFiles(std::string const& BasePath, ConfList const& Conf)
{
   std::string const EnergyLogFile = BasePath + ".e";
   int const EnergyLogLevel = Conf.Get("EnergyLogLevel", 0);
   Logger("EnergyLog").SetThreshold(EnergyLogLevel);
   if (EnergyLogLevel > 0)
   {
      EnergyLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(EnergyLogFile.c_str(), std::ios_base::out | std::ios_base::app));
      EnergyLog->precision(15);
      Logger("EnergyLog").SetStream(*EnergyLog);
   }

   std::string const SweepLogFile = BasePath + ".sweep";
   int const SweepLogLevel = Conf.Get("SweepLogLevel", 0);
   Logger("SweepLog").SetThreshold(SweepLogLevel);
   if (SweepLogLevel > 0)
   {
      SweepLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(SweepLogFile.c_str(), std::ios_base::out | std::ios_base::app));
      SweepLog->precision(15);
      Logger("SweepLog").SetStream(*SweepLog);
   }

   std::string const DensityLogFile = BasePath + ".density";
   int const DensityLogLevel = Conf.Get("DensityLogLevel", 0);
   Logger("DensityLog").SetThreshold(DensityLogLevel);
   if (DensityLogLevel > 0)
   {
      DensityLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(DensityLogFile.c_str(), std::ios_base::out | std::ios_base::app));
      DensityLog->precision(15);
      Logger("DensityLog").SetStream(*DensityLog);
   }

   // CpuLog
   // DiagLog

   this->ReadConfOptions(Conf);
}

void Solver::ReadConfOptions(ConfList const& Conf)
{
   LanczosMixFactor = Conf.Get("LanczosMixFactor", 1.0);
   AxMixFactor = Conf.Get("AxMixFactor", 0.0);
   TruncateModifyScale = Conf.Get("TruncateModifyScale", false);
   FlushLogs = Conf.Get("FlushLogs", false);
   TRACE("SimpleCG configuration settings")(LanczosMixFactor)(AxMixFactor)(TruncateModifyScale)(FlushLogs);
}

void Solver::StartIteration()
{
   IterationNumMultiplies = 0;
   IterationNumStates = 0;
   IterationTruncation = 0;
   IterationEntropy = 0;
}

void Solver::EndIteration()
{
   double IterationDOS_norm = Minus1OverPi * IterationGF_norm;
   double IterationDOS_overlap = Minus1OverPi * IterationGF_overlap;
   msg_log(1, "EnergyLog") << TotalSweepNumber << ' '
                           << TotalSweepRecNumber << ' '
                           << this->LeftSize() << ' '
                           << this->RightSize() << ' '
                           << IterationNumStates << ' '
                           << IterationNumMultiplies << ' '
                           << IterationTruncation << ' '
                           << IterationEntropy << ' '
                           << IterationDOS_norm << ' '
                           << IterationDOS_overlap << ' '
                           << IterationGF_real << ' '
                           << IterationSolverResid << '\n';

   if (FlushLogs)
      msg_log(1, "EnergyLog") << std::flush;

   ++SweepNumIterations;
   SweepNumMultiplies += IterationNumMultiplies;
   SweepSumStates += IterationNumStates;
   SweepMaxStates = std::max(SweepMaxStates, IterationNumStates);
   //   SweepEnergy = std::min(SweepEnergy, IterationEnergy);
   SweepTruncation += IterationTruncation;
   SweepEntropy = std::max(SweepEntropy, IterationEntropy);

   SweepGF_norm = IterationGF_norm;
   SweepGF_overlap = IterationGF_overlap;
   SweepGF_real = IterationGF_real;
}

void Solver::StartSweep(bool IncrementSweepNumber, double Broad_)
{
   ++TotalSweepNumber;
   if (IncrementSweepNumber)
      ++TotalSweepRecNumber;

   SweepNumIterations = 0;
   SweepSumStates = 0;
   SweepMaxStates = 0;
   SweepNumMultiplies = 0;
   SweepTruncation = 0;
   SweepEntropy = 0;

   if (Broad_ != 0)
      Broadening = Broad_;
}

void Solver::EndSweep()
{
   TotalNumIterations += SweepNumIterations;
   TotalNumMultiplies += SweepNumMultiplies;

   double AverageNumStates = double(SweepSumStates) / double(SweepNumIterations);

   double ExpectA2 = this->ExpectationA2();
   double SweepGF_functional = this->Functional(ExpectA2);
   double GlobalResid = this->ExactResidualNorm(ExpectA2);
   GlobalResid = GlobalResid*GlobalResid;  // we want the square of the residual norm

   int Converged = 0;

   double SweepDOS_norm = Minus1OverPi * SweepGF_norm;
   double SweepDOS_overlap = Minus1OverPi * SweepGF_overlap;

   msg_log(1, "SweepLog") << TotalSweepNumber << ' '
                          << TotalSweepRecNumber << ' '
                          << SweepNumIterations << ' '
                          << AverageNumStates << ' '
                          << SweepMaxStates << ' '
                          << SweepNumMultiplies << ' '
                          << Frequency << ' '
                          << Broadening << ' '
                          << SweepDOS_norm << ' '
                          << SweepDOS_overlap << ' '
                          << SweepGF_functional << ' '
                          << SweepGF_real << ' '
                          << SweepTruncation << ' '
                          << SweepEntropy << ' '
                          << GlobalResid << ' '
                          << Converged << ' '
      //                          << ExpectA2
                          << std::endl;
   msg_log(1, "EnergyLog") << std::flush;
}

double Solver::ResidualNorm() const
{
   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     x_A_x.Left(),
                                     x.Center(),
                                     herm(x_A_x.Right())) + complex(0.0, Broadening) * x.Center();

   //   MatrixOperator AxNorm = Ax * (1.0 / norm_frob(Ax));
   //   MatrixOperator yNorm = yprime * (1.0 / norm_frob(yprime));
   //   TRACE(norm_frob(AxNorm - yNorm));
   return norm_frob(yprime-Ax) / norm_frob(yprime);

}

CenterWavefunction Solver::WavefunctionAx() const
{
   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     x_A_x.Left(),
                                     x.Center(),
                                     herm(x_A_x.Right())) + complex(0.0, Broadening) * x.Center();

   CenterWavefunction Res(x); Res.Center() = Ax;
   return Res;
}

void Solver::ExpandLeft()
{
   DEBUG_TRACE("expanding")(this->ResidualNorm());

   this->DebugCheckBasis();

   DEBUG_TRACE(A.LeftSize())(x.LeftSize());

   // expand x
   x.Left() = prod(x.Left(), x.Center());
   x.Center() = ExpandBasis2(x.Left());
   //TRACE(norm_frob_sq(x.Center()));

   // matrix elements
   x_A_x.PopLeft();
   x_A_x.PushLeft(operator_prod(herm(A.Left()),
                                     herm(x.Left()),
                                     x_A_x.Left(),
                                     x.Left()));

   x_y.PopLeft();
   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));

   // yprime is y transformed into the new basis
   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));
   //TRACE(norm_frob_sq(x.Center()));

   DEBUG_TRACE("Expanding done")(this->ResidualNorm());

   this->DebugCheckBasis();
}

void Solver::ExpandRight()
{
   // expand x
   x.Right() = prod(x.Center(), x.Right());
   x.Center() = ExpandBasis1(x.Right());

   // matrix elements
   x_A_x.PopRight();
   x_A_x.PushRight(operator_prod(A.Right(),
                                 x.Right(),
                                 x_A_x.Right(),
                                 herm(x.Right())));

   x_y.PopRight();
   x_y.PushRight(operator_prod(x.Right(), x_y.Right(), herm(y.Right())));

   // x is just y transformed into the new basis
   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));
   //TRACE(norm_frob_sq(x.Center()));

   this->DebugCheckBasis();
}

void Solver::ShiftLeft()
{
   A.RotateLeft();
   y.RotateLeft();
   x.RotateLeft();

   // new matrix elements
   x_A_x.PushRight(operator_prod(A.Right(),
                                 x.Right(),
                                 x_A_x.Right(),
                                 herm(x.Right())));
   x_A_x.PopLeft();

   x_y.PushRight(operator_prod(x.Right(), x_y.Right(), herm(y.Right())));
   x_y.PopLeft();

   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));

   this->DebugCheckBasis();
}

void Solver::ShiftRight()
{
   DEBUG_TRACE("ShiftRight")(this->ResidualNorm());

   MatrixOperator Ay = operator_prod(conj(A.Center()),
                                     x_A_x.Left(),
                                     yprime,
                                     herm(x_A_x.Right()));
   TRACE(inner_prod(Ay, yprime));

   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     x_A_x.Left(),
                                     x.Center(),
                                     herm(x_A_x.Right()));

   MPStateComponent ypp = prod(yprime, x.Right());
   MPStateComponent Axp = prod(Ax, x.Right());

   A.RotateRight();
   y.RotateRight();
   x.RotateRight();

   x_A_x.PushLeft(operator_prod(herm(A.Left()),
                                herm(x.Left()),
                                x_A_x.Left(),
                                x.Left()));
   x_A_x.PopRight();

   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));
   x_y.PopRight();

   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));

   MatrixOperator yp = scalar_prod(herm(x.Left()), ypp);
   DEBUG_TRACE(norm_frob(yp))(norm_frob(yprime-yp))(inner_prod(yprime,yp));

   MatrixOperator Axx = operator_prod(conj(A.Center()),
                                      x_A_x.Left(),
                                      x.Center(),
                                      herm(x_A_x.Right()));
   MatrixOperator Ap = scalar_prod(herm(x.Left()), Axp);
   DEBUG_TRACE(norm_frob(Ap))(norm_frob(Axx - Ap))(inner_prod(Axx,Ap));
   DEBUG_TRACE(norm_frob(Axx - yprime) / norm_frob(yprime))(norm_frob(Ap - yprime) / norm_frob(yprime));

   Ay = operator_prod(conj(A.Center()),
                      x_A_x.Left(),
                      yprime,
                      herm(x_A_x.Right()));
   TRACE(inner_prod(Ay, yprime));

   DEBUG_TRACE("ShiftRight done")(this->ResidualNorm());

   this->DebugCheckBasis();
}

void Solver::ShiftRightAndExpand()
{
   DEBUG_TRACE("ShiftRightAndExpand")(this->ResidualNorm());

   // Rotate the wavefunctions to the right.  A and y are fixed and do not
   // need the basis expanded

   A.RotateRight();
   y.RotateRight();

   DEBUG_TRACE(norm_frob_sq(y.Center()));
   DEBUG_TRACE(norm_frob(x.Center()));
   DEBUG_TRACE(inner_prod(x.Center(), yprime));

   x.PushLeft(prod(x.Center(), x.Right()));
   x.PopRight();
   x.Center() = ExpandBasis2(x.Left());

   // new matrix elements
   x_A_x.PushLeft(operator_prod(herm(A.Left()),
                                herm(x.Left()),
                                x_A_x.Left(),
                                x.Left()));
   x_A_x.PopRight();

   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));
   x_y.PopRight();

   // yprime is y transformed into the new basis
   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));

   DEBUG_TRACE("ShiftRightAndExpand done")(this->ResidualNorm());

   this->DebugCheckBasis();
}

void Solver::ShiftLeftAndExpand()
{
   // correction
   DEBUG_TRACE("ShiftLeftAndExpand")(this->ResidualNorm());

   // Rotate the wavefunctions to the left.  A and y are fixed and do not
   // need the basis expanded
   A.RotateLeft();
   y.RotateLeft();

   DEBUG_TRACE(norm_frob_sq(y.Center()));

   x.PushRight(prod(x.Left(), x.Center()));
   x.PopLeft();
   x.Center() = ExpandBasis1(x.Right());

   // new matrix elements
   x_A_x.PushRight(operator_prod(A.Right(),
                                 x.Right(),
                                 x_A_x.Right(),
                                 herm(x.Right())));
   x_A_x.PopLeft();

   x_y.PushRight(operator_prod(x.Right(),
                               x_y.Right(),
                               herm(y.Right())));
   x_y.PopLeft();

   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));

   DEBUG_TRACE("ShiftLeftAndExpand done")(this->ResidualNorm());

   this->DebugCheckBasis();
}

template <typename T1, typename T2>
void MultiplyMatrices(T1& x, T2 const& y)
{
   std::size_t size = size1(x);
   for (std::size_t i = 0; i < size; ++i)
   {
      x(i,i) = x(i,i) * y(i,i);
   }
}

TruncationInfo Solver::TruncateLeft(StatesInfo const& SInfo, double CFactor)
{
   MatrixOperator yp = prod(x_y.Left(), y.Center(), y.Center().TransformsAs());
   MatrixOperator Rho_y = scalar_prod(yp, herm(yp));

   MatrixOperator Rho_x = scalar_prod(x.Center(), herm(x.Center()));

   if (TruncateModifyScale)
   {
      MatrixOperator Rho_x_R = scalar_prod(herm(x.Center()), x.Center());
      SimpleOperator Alpha = trace_prod(herm(x_A_x.Right()), prod(Rho_x_R, x_A_x.Right()));
      MatrixOperator Weights = operator_prod(herm(Alpha), herm(x_A_x.Left()), x_A_x.Left());

      // Move to a basis in which Rho_x is diagonal.
      MatrixOperator J = Tensor::Regularize(Rho_x.Basis1());
      Rho_x = triple_prod(J, Rho_x, herm(J));
      MatrixOperator U = prod(DiagonalizeHermitian(Rho_x), J, J.TransformsAs());
      Weights = triple_prod(U, Weights, herm(U));
      // modify the density matrix by the weight
      for (std::size_t i = 0; i < Rho_x.Basis1().size(); ++i)
      {
         MultiplyMatrices(Rho_x(i,i), Weights(i,i));
      }
      // Transform Rho_x back to the usual basis
      Rho_x = triple_prod(herm(U), Rho_x, U);
   }

   MatrixOperator DensityMat = (LanczosMixFactor / trace(Rho_y)) * Rho_y
      + (1.0 / trace(Rho_x)) * Rho_x;

   if (AxMixFactor != 0)
   {
      // calculate Ax - TODO: this was already done by the solver, should be saved somewhere
      MatrixOperator Ax = operator_prod(conj(A.Center()),
                                        x_A_x.Left(),
                                        x.Center(),
                                        herm(x_A_x.Right()));
      MatrixOperator Rho_Ax = scalar_prod(Ax, herm(Ax));
      DensityMat += (AxMixFactor / trace(Rho_Ax)) * Rho_Ax;
   }

   if (CFactor != 0)
   {
      MatrixOperator Correction = operator_prod(x_A_x.Left(), Rho_y, herm(x_A_x.Left()));
      MatrixOperator xAxInv = scalar_prod(x_A_x.Left(), herm(x_A_x.Left()));
      Tensor::InvertIrregularHPD(xAxInv);
      Correction = triple_prod(xAxInv, Correction, herm(xAxInv));
      DensityMat += (CFactor / trace(Correction)) * Correction;
   }

   DensityMatrix<MatrixOperator> DM(DensityMat);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(),
                                                                                           DM.end(),
                                                                                           SInfo,
                                                                                           Info));
   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   // truncate
   x.Center() = U * x.Center();
   x.Left() = prod(x.Left(), herm(U));

   x_A_x.Left() = triple_prod(U, x_A_x.Left(), herm(U));
   x_y.Left() = U * x_y.Left();

   yprime = U * yprime;

   LeftCorrection = MatrixOperator();

#if defined(DEBUG)
   TRACE(this->ResidualNorm());

   // debug
   Ax = prod(U, Ax, Ident);
   TRACE(inner_prod(Ax, yprime));
   yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));
   TRACE(inner_prod(Ax, yprime));
   TRACE(norm_frob(yprime-Ax) / norm_frob(yprime));  // this is our residual norm

   MatrixOperator Axx = operator_prod(conj(A.Center()),
                                      x_A_x.Left(),
                                      x.Center(),
                                      herm(x_A_x.Right()));
   TRACE(norm_frob(Ax - Axx))(norm_frob(Ax))(norm_frob(Axx));

   DEBUG_TRACE("Intermediate")(this->ResidualNorm());
  x_A_x.PopLeft();
   x_A_x.PushLeft(operator_prod(herm(A.Left()),
                                     herm(x.Left()),
                                     x_A_x.Left(),
                                     x.Left()));

   x_y.PopLeft();
   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));

   // debug: calculate the normalized difference between Axx and yprime
   DEBUG_TRACE(norm_frob(Axx*(1.0 / norm_frob(Axx)) - yprime*(1.0 / norm_frob(yprime))));

#endif
   DEBUG_TRACE("Truncating done")(this->ResidualNorm());

   this->DebugCheckBasis();
   return Info;
}

TruncationInfo Solver::TruncateRight(StatesInfo const& SInfo, double CFactor)
{
   DEBUG_TRACE(this->ResidualNorm());

   // calculate Ax - TODO: this was already done by the solver, should be saved somewhere

   MatrixOperator yp = y.Center() * herm(x_y.Right());
   MatrixOperator Rho_y = scalar_prod(herm(yp), yp);
   MatrixOperator Rho_x = scalar_prod(herm(x.Center()), x.Center());

   if (TruncateModifyScale)
   {
      MatrixOperator Rho_x_L = scalar_prod(x.Center(), herm(x.Center()));
      SimpleOperator Alpha = trace_prod(herm(x_A_x.Left()), prod(Rho_x_L, x_A_x.Left()));
      MatrixOperator Weights = operator_prod(Alpha, x_A_x.Right(), herm(x_A_x.Right()));

      // Move to a basis in which Rho_x is diagonal.
      MatrixOperator J = Tensor::Regularize(Rho_x.Basis1());
      Rho_x = triple_prod(J, Rho_x, herm(J));
      MatrixOperator U = prod(DiagonalizeHermitian(Rho_x), J, J.TransformsAs());
      Weights = triple_prod(U, Weights, herm(U));
      // modify the density matrix by the weight
      for (std::size_t i = 0; i < Rho_x.Basis1().size(); ++i)
      {
         MultiplyMatrices(Rho_x(i,i), Weights(i,i));
      }
      // Transform Rho_x back to the usual basis
      Rho_x = triple_prod(herm(U), Rho_x, U);
   }

   MatrixOperator DensityMat = (LanczosMixFactor / trace(Rho_y)) * Rho_y
      + (1.0 / trace(Rho_x)) * Rho_x;

   if (AxMixFactor != 0)
   {
      // calculate Ax - TODO: this was already done by the solver, should be saved somewhere
      MatrixOperator Ax = operator_prod(conj(A.Center()),
                                        x_A_x.Left(),
                                        x.Center(),
                                        herm(x_A_x.Right()));
      MatrixOperator Rho_Ax = scalar_prod(herm(Ax), Ax);
      DensityMat += (AxMixFactor / trace(Rho_Ax)) * Rho_Ax;
   }

   if (CFactor != 0)
   {
      MatrixOperator Correction = operator_prod(x_A_x.Right(), Rho_y, herm(x_A_x.Right()));
      MatrixOperator xAxInv = scalar_prod(herm(x_A_x.Right()), x_A_x.Right());
      Tensor::InvertIrregularHPD(xAxInv);
      Correction = triple_prod(xAxInv, Correction, herm(xAxInv));
      DensityMat += (CFactor / trace(Correction)) * Correction;
   }

   //   DensityMat *= (1.0 / 4.0);

   DEBUG_TRACE("right")(x.Center().Basis1().total_dimension())(x.Center().Basis2().total_dimension());
   DEBUG_TRACE(x.LeftSize())(DensityMat.Basis1().total_dimension());

   //   DensityMat = (1.0 / trace(Rho_x)) * Rho_x;

   DensityMatrix<MatrixOperator> DM(DensityMat);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(),
                                                                                           DM.end(),
                                                                                           SInfo,
                                                                                           Info));
   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   // truncate

   //   //TRACE(prod(U, adjoint(U), Ident));

#if defined(DEBUG)
   MatrixOperator UU = prod(adjoint(U), U, Ident);
   TRACE(norm_frob(x.Center() - prod(x.Center(), UU, Ident)));
   TRACE(norm_frob(Ax - prod(Ax, UU, Ident)));
   TRACE(norm_frob(yprime - prod(yprime, UU, Ident)));
#endif

   x.Center() = x.Center() * herm(U);
   x.Right() = prod(U, x.Right());

   x_A_x.Right() = triple_prod(U, x_A_x.Right(), herm(U));
   x_y.Right() = prod(U, x_y.Right(), Ident);

   yprime =yprime * herm(U);

   RightCorrection = MatrixOperator();

   DEBUG_TRACE(this->ResidualNorm());

   this->DebugCheckBasis();
   return Info;
}

void Solver::DebugCheckBasis() const
{
   DEBUG_CHECK_EQUAL(x_A_x.Left().Basis1(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_A_x.Right().Basis1(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_A_x.Left().Basis2(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_A_x.Right().Basis2(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_y.Left().Basis1(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_y.Right().Basis1(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_y.Left().Basis2(), y.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_y.Right().Basis2(), y.Center().Basis2());
}

PStream::opstream& operator<<(PStream::opstream& out, Solver const& d)
{
   return out << d.x
              << d.A
              << d.y
              << d.yprime
              << d.x_A_x
              << d.x_y
              << d.Ident
              << d.LeftCorrection
              << d.RightCorrection

              << d.Frequency
              << d.Broadening

              << d.TotalSweepNumber
              << d.TotalSweepRecNumber
              << d.TotalNumIterations
              << d.TotalNumMultiplies

              << d.SweepNumIterations
              << d.SweepSumStates
              << d.SweepMaxStates
              << d.SweepNumMultiplies
              << d.SweepTruncation
              << d.SweepEntropy
              << d.SweepGF_norm
              << d.SweepGF_overlap
              << d.SweepGF_real

              << d.IterationNumMultiplies
              << d.IterationNumStates
              << d.IterationTruncation
              << d.IterationEntropy
              << d.IterationGF_norm
              << d.IterationGF_overlap
              << d.IterationGF_real
              << d.IterationSolverResid;
}

PStream::ipstream& operator>>(PStream::ipstream& in, Solver& d)
{
   return in >> d.x
             >> d.A
             >> d.y
             >> d.yprime
             >> d.x_A_x
             >> d.x_y
             >> d.Ident
             >> d.LeftCorrection
             >> d.RightCorrection

             >> d.Frequency
             >> d.Broadening

             >> d.TotalSweepNumber
             >> d.TotalSweepRecNumber
             >> d.TotalNumIterations
             >> d.TotalNumMultiplies

             >> d.SweepNumIterations
             >> d.SweepSumStates
             >> d.SweepMaxStates
             >> d.SweepNumMultiplies
             >> d.SweepTruncation
             >> d.SweepEntropy
             >> d.SweepGF_norm
             >> d.SweepGF_overlap
             >> d.SweepGF_real

             >> d.IterationNumMultiplies
             >> d.IterationNumStates
             >> d.IterationTruncation
             >> d.IterationEntropy
             >> d.IterationGF_norm
             >> d.IterationGF_overlap
             >> d.IterationGF_real
             >> d.IterationSolverResid;
}
