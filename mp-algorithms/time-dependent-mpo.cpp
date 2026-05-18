// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/time-dependent-mpo.cpp
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
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

#include "time-dependent-mpo.h"
#include "lattice/infinite-parser.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

std::complex<double> const I(0.0, 1.0);
int const MaxImplementedMagnusOrder = 4;

std::pair<double, double>
LegendreAndDerivative(int Order, double x)
{
   if (Order == 0)
      return {1.0, 0.0};
   if (Order == 1)
      return {x, 1.0};

   double Pm2 = 1.0;
   double Pm1 = x;
   double P = Pm1;
   for (int n = 2; n <= Order; ++n)
   {
      P = ((2*n-1) * x * Pm1 - (n-1) * Pm2) / n;
      Pm2 = Pm1;
      Pm1 = P;
   }

   double dP = Order * (x * P - Pm2) / (x*x - 1.0);
   return {P, dP};
}

std::vector<std::vector<double>>
LagrangeCoefficients(std::vector<double> const& Nodes)
{
   int const Order = Nodes.size();
   std::vector<std::vector<double>> Coeff(Order, std::vector<double>(Order, 0.0));

   for (int i = 0; i < Order; ++i)
   {
      std::vector<double> Poly(1, 1.0);
      double Denom = 1.0;

      for (int j = 0; j < Order; ++j)
      {
         if (i == j)
            continue;

         std::vector<double> NewPoly(Poly.size()+1, 0.0);
         for (int k = 0; k < Poly.size(); ++k)
         {
            NewPoly[k] += -Nodes[j] * Poly[k];
            NewPoly[k+1] += Poly[k];
         }
         Poly = std::move(NewPoly);
         Denom *= Nodes[i] - Nodes[j];
      }

      for (int k = 0; k < Order; ++k)
         Coeff[k][i] = Poly[k] / Denom;
   }

   return Coeff;
}

std::vector<std::vector<double>>
CommutatorWeights(std::vector<double> const& Nodes,
                  std::vector<double> const& Weights)
{
   int const Order = Nodes.size();
   auto const Coeff = LagrangeCoefficients(Nodes);

   std::vector<std::vector<double>> Q(Order, std::vector<double>(Order, 0.0));
   for (int Row = 0; Row < Order; ++Row)
   {
      for (int Col = 0; Col < Order; ++Col)
      {
         double Power = Nodes[Row];
         for (int Degree = 0; Degree < Order; ++Degree)
         {
            Q[Row][Col] += Coeff[Degree][Col] * Power / double(Degree+1);
            Power *= Nodes[Row];
         }
      }
   }

   std::vector<std::vector<double>> Raw(Order, std::vector<double>(Order, 0.0));
   for (int Row = 0; Row < Order; ++Row)
      for (int Col = 0; Col < Order; ++Col)
         Raw[Row][Col] = Weights[Row] * Q[Row][Col]
                       - 0.5 * Weights[Row] * Weights[Col];

   std::vector<std::vector<double>> Result(Order, std::vector<double>(Order, 0.0));
   for (int Row = 0; Row < Order; ++Row)
   {
      for (int Col = 0; Col < Order; ++Col)
      {
         // Sign convention for the effective Hamiltonian.  For the two-point
         // Gauss rule this gives Gamma_12 = sqrt(3)/12, matching the existing
         // fourth-order implementation.
         Result[Row][Col] = 0.5 * (Raw[Col][Row] - Raw[Row][Col]);
      }
      Result[Row][Row] = 0.0;
   }

   return Result;
}

BasicTriangularMPO
WeightedSum(std::vector<BasicTriangularMPO> const& Operators,
            std::vector<double> const& Weights)
{
   BasicTriangularMPO Result = Weights[0] * Operators[0];
   for (int i = 1; i < Operators.size(); ++i)
      Result += Weights[i] * Operators[i];
   return Result;
}

} // namespace

void
ValidateMagnusOrder(int MagnusOrder)
{
   if (MagnusOrder < 2 || MagnusOrder % 2 != 0)
      throw std::runtime_error("Invalid Magnus order; expected a positive even integer.");
   if (MagnusOrder > MaxImplementedMagnusOrder)
      throw std::runtime_error("Magnus orders above " + std::to_string(MaxImplementedMagnusOrder)
                               + " are not implemented yet.");
}

int
DefaultMagnusQuadratureOrder(int MagnusOrder)
{
   ValidateMagnusOrder(MagnusOrder);

   return std::max(1, MagnusOrder / 2);
}

int
ResolveMagnusQuadratureOrder(int MagnusOrder, int MagnusQuadratureOrder)
{
   ValidateMagnusOrder(MagnusOrder);

   if (MagnusQuadratureOrder < 0)
      throw std::runtime_error("Invalid Magnus quadrature order; expected a non-negative integer.");

   int const Order = MagnusQuadratureOrder == 0
      ? DefaultMagnusQuadratureOrder(MagnusOrder)
      : MagnusQuadratureOrder;

   if (Order < 1)
      throw std::runtime_error("Invalid Magnus quadrature order; expected at least one point.");
   if (MagnusOrder >= 4 && Order < 2)
      throw std::runtime_error("Magnus orders of four or higher need at least two quadrature points.");

   return Order;
}

MagnusQuadrature
GaussLegendreMagnusQuadrature(int Order)
{
   if (Order < 1)
      throw std::runtime_error("Invalid Gauss-Legendre quadrature order; expected at least one point.");

   std::vector<double> Nodes(Order, 0.0);
   std::vector<double> Weights(Order, 0.0);

   double const Pi = std::acos(-1.0);
   int const RootCount = (Order + 1) / 2;
   for (int i = 0; i < RootCount; ++i)
   {
      double x = std::cos(Pi * (i + 0.75) / (Order + 0.5));
      for (int Iter = 0; Iter < 50; ++Iter)
      {
         auto const [P, dP] = LegendreAndDerivative(Order, x);
         double const dx = P / dP;
         x -= dx;
         if (std::abs(dx) <= 8.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(x)))
            break;
      }

      double const dP = LegendreAndDerivative(Order, x).second;
      double const Weight = 1.0 / ((1.0 - x*x) * dP * dP);

      Nodes[i] = 0.5 * (1.0 - x);
      Nodes[Order-1-i] = 0.5 * (1.0 + x);
      Weights[i] = Weight;
      Weights[Order-1-i] = Weight;
   }

   return {Nodes, Weights, CommutatorWeights(Nodes, Weights)};
}

BasicTriangularMPO
TimeDependentHamiltonianMPO(InfiniteLattice const& Lattice,
                            std::string const& HamOperator,
                            std::string const& TimeVar,
                            std::complex<double> Time,
                            std::complex<double> Timestep,
                            int MagnusOrder,
                            int MagnusQuadratureOrder)
{
   ValidateMagnusOrder(MagnusOrder);

   if (Timestep == 0.0)
      return ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, Time}});

   int const QuadratureOrder = ResolveMagnusQuadratureOrder(MagnusOrder, MagnusQuadratureOrder);
   MagnusQuadrature const Quadrature = GaussLegendreMagnusQuadrature(QuadratureOrder);

   std::vector<BasicTriangularMPO> Samples;
   Samples.reserve(QuadratureOrder);
   for (double Node : Quadrature.Nodes)
      Samples.push_back(ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, Time + Node * Timestep}}));

   BasicTriangularMPO Result = WeightedSum(Samples, Quadrature.Weights);

   if (MagnusOrder >= 4)
   {
      for (int i = 0; i < QuadratureOrder; ++i)
      {
         for (int j = i+1; j < QuadratureOrder; ++j)
         {
            double const Gamma = Quadrature.CommutatorWeights[i][j];
            if (Gamma != 0.0)
               Result += (-I * Timestep * Gamma) * commutator(Samples[i], Samples[j]);
         }
      }
   }

   return Result;
}
