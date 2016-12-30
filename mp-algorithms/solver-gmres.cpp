// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/solver-gmres.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "solver-gmres.h"

typedef std::complex<double> complex;

PStream::opstream& operator<<(PStream::opstream& out, SolverGmres const& s)
{
   return out << static_cast<Solver const&>(s)
              << s.ASquared;
}

PStream::ipstream& operator>>(PStream::ipstream& in, SolverGmres& s)
{
   return in >> static_cast<Solver&>(s)
             >> s.ASquared;
}

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiply(SimpleOperator const& Op_,
                      MPStateComponent const& Left_,
                      MPStateComponent const& Right_, double Broadening_)
      : Op(Op_), Left(Left_), Right(Right_), Broadening(Broadening_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, Left, Psi, herm(Right)) + complex(0.0, Broadening) * Psi;
   }

   SimpleOperator Op;
   MPStateComponent Left, Right;
   double Broadening;
};

struct CompoundInversePrecondition
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   CompoundInversePrecondition(MPStateComponent const& Left_, MPStateComponent const& Right_,
                               MatrixOperator const& Ei, MatrixOperator const& Fi)
      : y_A_x_Left(Left_), y_A_x_Right(Right_), EInv(Ei), FInv(Fi) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      //      MatrixOperator Temp = operator_prod(herm(y_A_x_Left), r, y_A_x_Right);
      //      return triple_prod(EInv, Temp, herm(FInv));
      return operator_prod(herm(y_A_x_Left), triple_prod(EInv, r, herm(FInv)), y_A_x_Right);
   }

   MPStateComponent const& y_A_x_Left;
   MPStateComponent const& y_A_x_Right;
   MatrixOperator const& EInv;
   MatrixOperator const& FInv;
};

struct InvertScalarPrecondition
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   InvertScalarPrecondition(MPStateComponent const& Ei, MPStateComponent const& Fi)
      : iL(Ei), iR(Fi) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      return operator_prod(iL, r, herm(iR));
   }

   MPStateComponent const& iL;
   MPStateComponent const& iR;
};

struct TripleProdPrecondition
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   TripleProdPrecondition(MatrixOperator const& Ei, MatrixOperator const& Fi)
      : iL(Ei), iR(Fi) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      return triple_prod(iL, r, herm(iR));
   }

   MatrixOperator const& iL;
   MatrixOperator const& iR;
};

void SolverGmres::ReadConfOptions(ConfList const& Conf)
{
   Solver::ReadConfOptions(Conf);
   Tol = Conf.Get("GMRES::Tol", 0.0);
   SubspaceSize = Conf.Get("GMRES::SubspaceSize", 20);
   TRACE("GMRES conf options:")(Tol)(SubspaceSize);
}

double SolverGmres::Solve(int MaxIterations)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   int Iter = MaxIterations;
   int m = SubspaceSize;
   double Resid = Tol;

   std::string Preconditioner = "None";

   if (Preconditioner == "CompoundInverse")
   {
      MatrixOperator E = scalar_prod(herm(x_A_x.Left()), x_A_x.Left());
      MatrixOperator F = scalar_prod(herm(x_A_x.Right()), x_A_x.Right());
      //      E += MatrixOperator::make_identity(E.Basis1());
      //      F += MatrixOperator::make_identity(F.Basis1());

      MatrixOperator UE = Tensor::Regularize(E.Basis1());
      E = triple_prod(UE, E, herm(UE));
      InvertHPD(E);
      E = triple_prod(herm(UE), E, UE);

      MatrixOperator UF = Tensor::Regularize(F.Basis1());
      F = triple_prod(UF, F, herm(UF));
      InvertHPD(F);
      F = triple_prod(herm(UF), F, UF);

      GmRes(x.Center(),
            SuperblockMultiply(conj(A.Center()),
                               x_A_x.Left(),
                               x_A_x.Right(), Broadening),
            yprime,
            m, Iter, Tol,
            CompoundInversePrecondition(x_A_x.Left(), x_A_x.Right(), E, F));
   }
   else if (Preconditioner == "None")
   {
      GmRes(x.Center(),
            SuperblockMultiply(conj(A.Center()),
                               x_A_x.Left(),
                               x_A_x.Right(), Broadening),
            yprime,
            m, Iter, Resid,
            LinearAlgebra::Identity<MatrixOperator>());
   }
   else if (Preconditioner == "InvertScalar")
   {
      MPStateComponent E = x_A_x.Left();
      MPStateComponent F = x_A_x.Right();

      MatrixOperator UE = Tensor::Regularize(E.Basis1());
      for (std::size_t i = 0; i < E.size(); ++i)
      {
         if (is_scalar(E[i].TransformsAs()))
         {
            MatrixOperator x = triple_prod(UE, E[i], herm(UE));
            Tensor::InvertGeneral(x);
            E[i] = triple_prod(herm(UE), x, UE);
         }
         else
         {
            E[i] = MatrixOperator(E.Basis1(), E.Basis2(), E[i].TransformsAs());
         }
      }

      MatrixOperator UF = Tensor::Regularize(F.Basis1());
      for (std::size_t i = 0; i < F.size(); ++i)
      {
         if (is_scalar(F[i].TransformsAs()))
         {
            MatrixOperator x = triple_prod(UF, F[i], herm(UF));
            Tensor::InvertGeneral(x);
            F[i] = triple_prod(herm(UF), x, UF);
         }
         else
         {
            F[i] = MatrixOperator(F.Basis1(), F.Basis2(), F[i].TransformsAs());
         }
      }

      GmRes(x.Center(),
            SuperblockMultiply(conj(A.Center()),
                               x_A_x.Left(),
                               x_A_x.Right(), Broadening),
            yprime,
            m, Iter, Tol,
            InvertScalarPrecondition(E, F));
   }
   else if (Preconditioner == "InverseSqrt")
   {
      //      MatrixOperator E = scalar_prod(herm(x_A_x.Left()), x_A_x.Left());
      //      MatrixOperator F = scalar_prod(herm(x_A_x.Right()), x_A_x.Right());
      MatrixOperator E = scalar_prod(x_A_x.Left(), herm(x_A_x.Left()));
      MatrixOperator F = scalar_prod(x_A_x.Right(), herm(x_A_x.Right()));
      //      E += MatrixOperator::make_identity(E.Basis1());
      //      F += MatrixOperator::make_identity(F.Basis1());

      MatrixOperator UE = Tensor::Regularize(E.Basis1());
      E = triple_prod(UE, E, herm(UE));
      MatrixOperator ME = DiagonalizeHermitian(E);
      for (std::size_t i = 0; i < E.Basis1().size(); ++i)
      {
         for (std::size_t j = 0; j < size1(E(i,i)); ++j)
         {
            E(i,i)(j,j) = 1.0 / std::sqrt(E(i,i)(j,j));
         }
      }
      E = triple_prod(herm(ME), E, ME);
      E = triple_prod(herm(UE), E, UE);

      MatrixOperator UF = Tensor::Regularize(F.Basis1());
      F = triple_prod(UF, F, herm(UF));
      MatrixOperator MF = DiagonalizeHermitian(F);
      for (std::size_t i = 0; i < F.Basis1().size(); ++i)
      {
         for (std::size_t j = 0; j < size1(F(i,i)); ++j)
         {
            F(i,i)(j,j) = 1.0 / std::sqrt(F(i,i)(j,j));
         }
      }
      F = triple_prod(herm(MF), F, MF);
      F = triple_prod(herm(UF), F, UF);

      GmRes(x.Center(),
            SuperblockMultiply(conj(A.Center()),
                               x_A_x.Left(),
                               x_A_x.Right(), Broadening),
            yprime,
            m, Iter, Tol,
            TripleProdPrecondition(E, F));
   }
   else
   {
      PANIC("Unknown preconditioner")(Preconditioner);
   }

   IterationNumMultiplies = Iter;
   IterationSolverResid = Resid;
   //TRACE(Resid);
   IterationGF_norm = -Broadening * norm_frob_sq(x.Center());
   std::complex<double> g = inner_prod(x.Center(), yprime);
   // ** conjugation here, to account for the solver conjugate bug
   IterationGF_overlap = -imag(g);
   IterationGF_real = real(g);

   return IterationGF_overlap;
}

std::complex<double> SolverGmres::Overlap() const
{
   return expectation(x, A, y);
#if 0
   return inner_prod(operator_prod(A.Center(),
                                   x_A_x.Left(),
                                   x.Center(),
                                   herm(x_A_x.Right())) + complex(0.0, Broadening) * x.Center(),
                     yprime);
#endif
}

std::complex<double> SolverGmres::GreensFunction() const
{
   return inner_prod(yprime, x.Center());
}

double SolverGmres::Functional(double ExpectA2) const
{
#if defined(NDEBUG) || 1
   return 0;
#else
#error "This code needs updating now that the broadening is not included in A"
   CenterWavefunction xBar = CenterWavefunction(conj(x.AsLinearWavefunction()));
   CenterWavefunction xImag =
      CenterWavefunction(0.5 * (x.AsLinearWavefunction() - xBar.AsLinearWavefunction()));
   //   std::complex<double> ExpectA2Bar = expectation(x, ASquared, xBar);
   //   double SquarePart = 0.5 * (ExpectA2 - ExpectA2Bar.real());
   double SquarePart = expectation(xImag, ASquared, xImag).real();
   //   TRACE(ExpectA2)(expectation(x, ASquared, xBar))(expectation(x, ASquared, x))(SquarePart);

   TRACE(ExpectA2)(SquarePart)(inner_prod(x.Center(), yprime))(Frequency)(Broadening);
   //   TRACE(expectation(xBar, ASquared, xBar));

   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     x_A_x.Left(),
                                     x.Center(),
                                     herm(x_A_x.Right()));
   TRACE(norm_frob_sq(Ax));

   return Minus1OverPi*((1.0/Broadening)*SquarePart + 2*imag(inner_prod(x.Center(), yprime)));
#endif
}

double SolverGmres::ExpectationA2() const
{
   return expectation(x, ASquared, x).real();
}

double SolverGmres::ExactResidualNorm(double ExpectA2) const
{
#if 0
   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     x_A_x.Left(),
                                     x.Center(),
                                     herm(x_A_x.Right())) + Broadening*Broadening*x.Center();
   CenterWavefunction AxP(x);
   AxP.Center() = Ax;
   double yAx = 2 * overlap(AxP, y).real();
#endif

   double yNorm = norm_frob_sq(y.Center());
   double xNorm = norm_frob_sq(x.Center());
   std::complex<double> xyOverlap = inner_prod(x.Center(), yprime);
   std::complex<double> yAx = this->Overlap();
   std::complex<double> xAx = inner_prod(operator_prod(conj(A.Center()),
                                                       x_A_x.Left(),
                                                       x.Center(),
                                                       herm(x_A_x.Right())), x.Center());

   double Resid2 = yNorm + ExpectA2 + Broadening * Broadening * xNorm
      - 2 * yAx.real() - 2 * Broadening * xyOverlap.imag();
   double Norm = yNorm * (ExpectA2 + Broadening * Broadening * xNorm);
   TRACE_IF(Resid2 < 0)("Residual is less than zero?!?!")(Resid2);
   TRACE_IF(Norm < 0)("Normalization is less than zero?!?!")(Norm);
   return std::sqrt(Resid2 / std::sqrt(Norm));
}
