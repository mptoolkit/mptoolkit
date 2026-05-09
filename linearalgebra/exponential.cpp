// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/exponential.cpp
//
// Copyright (C) 2006-2016 Ian McCulloch <ian@qusim.net>
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
//
// implementation of the dense matrix exponential backend
//
// Created 2006-05-12 Ian McCulloch
//

#include "eigen.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

namespace LinearAlgebra
{

namespace Private
{

namespace
{

constexpr std::array<double, 5> ThetaBounds{
   1.495585217958292e-2, 2.539398330063230e-1, 9.504178996162932e-1,
   2.097847961257068, 5.371920351148152};

// Adaptive scaling-and-squaring Pade implementation, ported from uni20.
// This replaces the historical Fortran EXPOKIT zgpadm dependency.
class ExponentialMatrix
{
   public:
      using value_type = std::complex<double>;

      ExponentialMatrix() = default;

      ExponentialMatrix(std::size_t Rows, std::size_t Cols)
         : Rows_(Rows), Cols_(Cols), Data_(Rows * Cols) {}

      std::size_t rows() const { return Rows_; }
      std::size_t cols() const { return Cols_; }
      std::size_t size() const { return Data_.size(); }

      value_type& operator()(std::size_t Row, std::size_t Col)
      {
         return Data_[Row * Cols_ + Col];
      }

      value_type const& operator()(std::size_t Row, std::size_t Col) const
      {
         return Data_[Row * Cols_ + Col];
      }

      value_type* data() { return Data_.data(); }
      value_type const* data() const { return Data_.data(); }

   private:
      std::size_t Rows_ = 0;
      std::size_t Cols_ = 0;
      std::vector<value_type> Data_;
};

ExponentialMatrix
make_identity(std::size_t n)
{
   ExponentialMatrix Result(n, n);
   for (std::size_t i = 0; i < n; ++i)
   {
      Result(i, i) = 1.0;
   }
   return Result;
}

ExponentialMatrix
multiply(ExponentialMatrix const& lhs, ExponentialMatrix const& rhs)
{
   CHECK_EQUAL(lhs.cols(), rhs.rows());

   ExponentialMatrix Result(lhs.rows(), rhs.cols());
   for (std::size_t i = 0; i < lhs.rows(); ++i)
   {
      for (std::size_t k = 0; k < lhs.cols(); ++k)
      {
         std::complex<double> const Factor = lhs(i, k);
         if (Factor == std::complex<double>())
            continue;

         for (std::size_t j = 0; j < rhs.cols(); ++j)
         {
            Result(i, j) += Factor * rhs(k, j);
         }
      }
   }
   return Result;
}

ExponentialMatrix
add(ExponentialMatrix const& lhs, ExponentialMatrix const& rhs)
{
   CHECK_EQUAL(lhs.rows(), rhs.rows());
   CHECK_EQUAL(lhs.cols(), rhs.cols());

   ExponentialMatrix Result(lhs.rows(), lhs.cols());
   for (std::size_t i = 0; i < lhs.size(); ++i)
   {
      Result.data()[i] = lhs.data()[i] + rhs.data()[i];
   }
   return Result;
}

ExponentialMatrix
subtract(ExponentialMatrix const& lhs, ExponentialMatrix const& rhs)
{
   CHECK_EQUAL(lhs.rows(), rhs.rows());
   CHECK_EQUAL(lhs.cols(), rhs.cols());

   ExponentialMatrix Result(lhs.rows(), lhs.cols());
   for (std::size_t i = 0; i < lhs.size(); ++i)
   {
      Result.data()[i] = lhs.data()[i] - rhs.data()[i];
   }
   return Result;
}

ExponentialMatrix
scale(ExponentialMatrix const& Mat, std::complex<double> const& Scalar)
{
   ExponentialMatrix Result(Mat.rows(), Mat.cols());
   for (std::size_t i = 0; i < Mat.size(); ++i)
   {
      Result.data()[i] = Mat.data()[i] * Scalar;
   }
   return Result;
}

long double
matrix_one_norm(ExponentialMatrix const& Mat)
{
   long double Result = 0.0L;
   for (std::size_t j = 0; j < Mat.cols(); ++j)
   {
      long double ColumnSum = 0.0L;
      for (std::size_t i = 0; i < Mat.rows(); ++i)
      {
         ColumnSum += static_cast<long double>(std::abs(Mat(i, j)));
      }
      Result = std::max(Result, ColumnSum);
   }
   return Result;
}

void
swap_rows(ExponentialMatrix& Mat, std::size_t lhs, std::size_t rhs)
{
   if (lhs == rhs)
      return;

   for (std::size_t j = 0; j < Mat.cols(); ++j)
   {
      std::swap(Mat(lhs, j), Mat(rhs, j));
   }
}

ExponentialMatrix
solve_linear_system(ExponentialMatrix A, ExponentialMatrix B)
{
   CHECK_EQUAL(A.rows(), A.cols());
   CHECK_EQUAL(A.rows(), B.rows());

   std::size_t const n = A.rows();
   std::size_t const nrhs = B.cols();

   for (std::size_t k = 0; k < n; ++k)
   {
      std::size_t PivotRow = k;
      double PivotValue = std::abs(A(k, k));
      for (std::size_t i = k + 1; i < n; ++i)
      {
         double Candidate = std::abs(A(i, k));
         if (Candidate > PivotValue)
         {
            PivotValue = Candidate;
            PivotRow = i;
         }
      }

      CHECK(PivotValue != 0.0)("singular matrix in matrix exponential");

      swap_rows(A, k, PivotRow);
      swap_rows(B, k, PivotRow);

      std::complex<double> const Pivot = A(k, k);
      for (std::size_t i = k + 1; i < n; ++i)
      {
         std::complex<double> const Factor = A(i, k) / Pivot;
         if (Factor == std::complex<double>())
            continue;

         A(i, k) = {};
         for (std::size_t j = k + 1; j < n; ++j)
         {
            A(i, j) -= Factor * A(k, j);
         }
         for (std::size_t j = 0; j < nrhs; ++j)
         {
            B(i, j) -= Factor * B(k, j);
         }
      }
   }

   ExponentialMatrix X = B;
   for (int i = static_cast<int>(n) - 1; i >= 0; --i)
   {
      auto const Row = static_cast<std::size_t>(i);
      std::complex<double> const Pivot = A(Row, Row);
      for (std::size_t j = 0; j < nrhs; ++j)
      {
         std::complex<double> Value = X(Row, j);
         for (std::size_t k = Row + 1; k < n; ++k)
         {
            Value -= A(Row, k) * X(k, j);
         }
         X(Row, j) = Value / Pivot;
      }
   }

   return X;
}

std::complex<double>
to_scalar(double Value)
{
   return std::complex<double>(Value);
}

ExponentialMatrix
linear_combination(std::size_t Rows, std::size_t Cols,
                   std::initializer_list<std::pair<ExponentialMatrix const*, std::complex<double>>> Terms)
{
   ExponentialMatrix Result(Rows, Cols);
   for (auto const& Term : Terms)
   {
      ExponentialMatrix const& Mat = *Term.first;
      std::complex<double> const Coefficient = Term.second;
      for (std::size_t i = 0; i < Rows; ++i)
      {
         for (std::size_t j = 0; j < Cols; ++j)
         {
            Result(i, j) += Coefficient * Mat(i, j);
         }
      }
   }
   return Result;
}

ExponentialMatrix
solve_pade(ExponentialMatrix const& U, ExponentialMatrix const& V)
{
   ExponentialMatrix Numerator = add(V, U);
   ExponentialMatrix Denominator = subtract(V, U);
   return solve_linear_system(std::move(Denominator), std::move(Numerator));
}

ExponentialMatrix
pade3(ExponentialMatrix const& A)
{
   static constexpr std::array<double, 4> b{120.0, 60.0, 12.0, 1.0};
   std::size_t const n = A.rows();
   ExponentialMatrix const Identity = make_identity(n);
   ExponentialMatrix const A2 = multiply(A, A);

   ExponentialMatrix const Tmp =
      linear_combination(n, n, {{&A2, to_scalar(b[3])}, {&Identity, to_scalar(b[1])}});
   ExponentialMatrix const U = multiply(A, Tmp);

   ExponentialMatrix const V =
      linear_combination(n, n, {{&A2, to_scalar(b[2])}, {&Identity, to_scalar(b[0])}});

   return solve_pade(U, V);
}

ExponentialMatrix
pade5(ExponentialMatrix const& A)
{
   static constexpr std::array<double, 6> b{30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0};
   std::size_t const n = A.rows();
   ExponentialMatrix const Identity = make_identity(n);
   ExponentialMatrix const A2 = multiply(A, A);
   ExponentialMatrix const A4 = multiply(A2, A2);

   ExponentialMatrix const Tmp = linear_combination(n, n,
      {{&A4, to_scalar(b[5])}, {&A2, to_scalar(b[3])}, {&Identity, to_scalar(b[1])}});
   ExponentialMatrix const U = multiply(A, Tmp);

   ExponentialMatrix const V = linear_combination(n, n,
      {{&A4, to_scalar(b[4])}, {&A2, to_scalar(b[2])}, {&Identity, to_scalar(b[0])}});

   return solve_pade(U, V);
}

ExponentialMatrix
pade7(ExponentialMatrix const& A)
{
   static constexpr std::array<double, 8> b{
      17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0};
   std::size_t const n = A.rows();
   ExponentialMatrix const Identity = make_identity(n);
   ExponentialMatrix const A2 = multiply(A, A);
   ExponentialMatrix const A4 = multiply(A2, A2);
   ExponentialMatrix const A6 = multiply(A4, A2);

   ExponentialMatrix const Tmp = linear_combination(n, n,
      {{&A6, to_scalar(b[7])}, {&A4, to_scalar(b[5])}, {&A2, to_scalar(b[3])},
       {&Identity, to_scalar(b[1])}});
   ExponentialMatrix const U = multiply(A, Tmp);

   ExponentialMatrix const V = linear_combination(n, n,
      {{&A6, to_scalar(b[6])}, {&A4, to_scalar(b[4])}, {&A2, to_scalar(b[2])},
       {&Identity, to_scalar(b[0])}});

   return solve_pade(U, V);
}

ExponentialMatrix
pade9(ExponentialMatrix const& A)
{
   static constexpr std::array<double, 10> b{
      17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0,
      2162160.0,     110880.0,     3960.0,       90.0,        1.0};
   std::size_t const n = A.rows();
   ExponentialMatrix const Identity = make_identity(n);
   ExponentialMatrix const A2 = multiply(A, A);
   ExponentialMatrix const A4 = multiply(A2, A2);
   ExponentialMatrix const A6 = multiply(A4, A2);
   ExponentialMatrix const A8 = multiply(A6, A2);

   ExponentialMatrix const Tmp = linear_combination(n, n,
      {{&A8, to_scalar(b[9])}, {&A6, to_scalar(b[7])}, {&A4, to_scalar(b[5])},
       {&A2, to_scalar(b[3])}, {&Identity, to_scalar(b[1])}});
   ExponentialMatrix const U = multiply(A, Tmp);

   ExponentialMatrix const V = linear_combination(n, n,
      {{&A8, to_scalar(b[8])}, {&A6, to_scalar(b[6])}, {&A4, to_scalar(b[4])},
       {&A2, to_scalar(b[2])}, {&Identity, to_scalar(b[0])}});

   return solve_pade(U, V);
}

ExponentialMatrix
pade13(ExponentialMatrix const& A, ExponentialMatrix const& A2,
       ExponentialMatrix const& A4, ExponentialMatrix const& A6)
{
   static constexpr std::array<double, 14> b{
      64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
      1187353796428800.0, 129060195264000.0,   10559470521600.0,
      670442572800.0,     33522128640.0,       1323241920.0,
      40840800.0,         960960.0,            16380.0,
      182.0,              1.0};

   std::size_t const n = A.rows();
   ExponentialMatrix const Identity = make_identity(n);

   ExponentialMatrix const First = linear_combination(n, n,
      {{&A6, to_scalar(b[13])}, {&A4, to_scalar(b[11])}, {&A2, to_scalar(b[9])}});
   ExponentialMatrix Tmp = multiply(A6, First);
   ExponentialMatrix const Second = linear_combination(n, n,
      {{&A6, to_scalar(b[7])}, {&A4, to_scalar(b[5])}, {&A2, to_scalar(b[3])},
       {&Identity, to_scalar(b[1])}});
   Tmp = add(Tmp, Second);
   ExponentialMatrix const U = multiply(A, Tmp);

   ExponentialMatrix const Third = linear_combination(n, n,
      {{&A6, to_scalar(b[12])}, {&A4, to_scalar(b[10])}, {&A2, to_scalar(b[8])}});
   ExponentialMatrix V = multiply(A6, Third);
   ExponentialMatrix const Fourth = linear_combination(n, n,
      {{&A6, to_scalar(b[6])}, {&A4, to_scalar(b[4])}, {&A2, to_scalar(b[2])},
       {&Identity, to_scalar(b[0])}});
   V = add(V, Fourth);

   return solve_pade(U, V);
}

int
compute_scaling_exponent(ExponentialMatrix const& A4, ExponentialMatrix const& A6)
{
   long double const Norm4 = matrix_one_norm(A4);
   long double const Norm6 = matrix_one_norm(A6);
   CHECK(std::isfinite(Norm4))(Norm4);
   CHECK(std::isfinite(Norm6))(Norm6);

   long double const D4 = std::pow(Norm4, 0.25L);
   long double const D6 = std::pow(Norm6, 1.0L / 6.0L);
   long double const Eta = std::max(D4, D6);
   if (Eta == 0.0L)
      return 0;

   long double const Ratio = Eta / static_cast<long double>(ThetaBounds.back());
   if (Ratio <= 1.0L)
      return 0;

   long double const Exponent = std::log2(Ratio);
   if (Exponent <= 0.0L)
      return 0;
   return static_cast<int>(std::ceil(Exponent));
}

int
compute_prescaling_exponent(long double NormA, std::size_t Size)
{
   CHECK(std::isfinite(NormA));

   long double const Dimension = std::max(1.0L, static_cast<long double>(Size));
   long double const SafePowerNorm =
      0.5L * std::pow(static_cast<long double>(std::numeric_limits<double>::max()),
                      1.0L / 6.0L)
      / std::pow(Dimension, 5.0L / 6.0L);
   if (NormA <= SafePowerNorm)
      return 0;

   long double const Exponent = std::log2(NormA / SafePowerNorm);
   if (Exponent <= 0.0L)
      return 0;

   return static_cast<int>(std::ceil(Exponent));
}

ExponentialMatrix
expm(ExponentialMatrix const& Matrix, double t)
{
   CHECK_EQUAL(Matrix.rows(), Matrix.cols());

   if (Matrix.rows() == 0)
      return ExponentialMatrix();

   std::size_t const n = Matrix.rows();
   ExponentialMatrix A = scale(Matrix, std::complex<double>(t));

   long double const NormA = matrix_one_norm(A);
   CHECK(std::isfinite(NormA))(NormA);

   if (NormA == 0.0L)
      return make_identity(n);

   if (NormA <= static_cast<long double>(ThetaBounds[0]))
      return pade3(A);
   if (NormA <= static_cast<long double>(ThetaBounds[1]))
      return pade5(A);
   if (NormA <= static_cast<long double>(ThetaBounds[2]))
      return pade7(A);
   if (NormA <= static_cast<long double>(ThetaBounds[3]))
      return pade9(A);

   int const PreScaling = compute_prescaling_exponent(NormA, n);
   double const PreScale = std::ldexp(1.0, -PreScaling);
   ExponentialMatrix const PowerA =
      PreScaling == 0 ? A : scale(A, std::complex<double>(PreScale));

   ExponentialMatrix const A2 = multiply(PowerA, PowerA);
   ExponentialMatrix const A4 = multiply(A2, A2);
   ExponentialMatrix const A6 = multiply(A4, A2);

   int const AdditionalScaling = compute_scaling_exponent(A4, A6);
   int const s = PreScaling + AdditionalScaling;
   double const Scale = std::ldexp(1.0, -AdditionalScaling);
   ExponentialMatrix const ScaledA =
      AdditionalScaling == 0 ? PowerA : scale(PowerA, std::complex<double>(Scale));

   ExponentialMatrix A2Scaled;
   ExponentialMatrix A4Scaled;
   ExponentialMatrix A6Scaled;
   if (AdditionalScaling == 0)
   {
      A2Scaled = A2;
      A4Scaled = A4;
      A6Scaled = A6;
   }
   else
   {
      double const ScaleSq = Scale * Scale;
      double const ScalePow4 = ScaleSq * ScaleSq;
      double const ScalePow6 = ScalePow4 * ScaleSq;
      A2Scaled = scale(A2, std::complex<double>(ScaleSq));
      A4Scaled = scale(A4, std::complex<double>(ScalePow4));
      A6Scaled = scale(A6, std::complex<double>(ScalePow6));
   }

   ExponentialMatrix Result = pade13(ScaledA, A2Scaled, A4Scaled, A6Scaled);
   for (int k = 0; k < s; ++k)
   {
      Result = multiply(Result, Result);
   }

   return Result;
}

} // namespace

void Exponentiate(double t, int Size, std::complex<double> const* H, int ldH,
                  std::complex<double>* R, int ldR)
{
   ExponentialMatrix Input(Size, Size);
   for (int i = 0; i < Size; ++i)
   {
      for (int j = 0; j < Size; ++j)
      {
         Input(i, j) = H[i * ldH + j];
      }
   }

   ExponentialMatrix Result = expm(Input, t);
   for (int i = 0; i < Size; ++i)
   {
      std::memcpy(R + i * ldR, Result.data() + i * Size,
                  Size * sizeof(std::complex<double>));
   }
}

} // namespace Private

} // namespace LinearAlgebra
