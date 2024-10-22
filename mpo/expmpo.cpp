// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/basic_triangular_mpo.h
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
// Copyright (C) 2024 Ian McCulloch <ian@qusim.net>
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

#include "expmpo.h"

// These are not needed in C++17 (?), but C++14 requires the initializer
constexpr std::array<char const*, 3> ExpMpoSchemeTraits::Names;

double factorial(int x)
{
   if (x == 1) return 1;
   else return x * factorial(x-1);
}

// Calculate the aux_tensor_product of an array of pointers
OperatorComponent AuxTensorProdArray(std::vector<OperatorComponent const*> const& x)
{
   CHECK(x.size() >= 1);

   OperatorComponent Result = *x[0];
   for (int n = 1; n < x.size(); ++n)
   {
      Result = aux_tensor_prod(Result, *x[n]);
   }
   return Result;
}

// (1/number_of_permutations) permutations of aux_tensor_prod() of a vector of pointers to OperatorComponent.
// We use pointers so that duplicates can be identified
// This implements the operation {x[0] x[1] ... x[n]}, which is the sum of all permutations of x,
// divided by the number of permutations. We can take advantade of duplicates to reduce the number
// of permutations that we need.
OperatorComponent SumPerms(std::vector<OperatorComponent const*> x)
{
   std::sort(x.begin(), x.end());
   int n = 1;
   OperatorComponent Result = AuxTensorProdArray(x);
   while (std::next_permutation(x.begin(), x.end()))
   {
      Result += AuxTensorProdArray(x);
      ++n;
   }
   Result *= 1.0 / n;
   return Result;
}

// Constructs the sum of tensor products of all permutations of a list of operators, including D up to Order times.
// This corresponds to the series expansion
// DysonSeries(D, {X}) = X + {DX}/2! + {DDX}/3! + {DDDX}/4! + ....
// DysonSeries(D, {X,Y}) = {XY} + {DXY}/2! + {DDXY}/3! + ....
OperatorComponent DysonSeries(OperatorComponent const& D, int Order, std::initializer_list<OperatorComponent> list)
{
   std::vector<OperatorComponent const*> v;
   v.reserve(list.size()+Order);
   for (auto const& item : list)
   {
      v.push_back(&item);
   }

   OperatorComponent Result = SumPerms(v);
   for (int n = 1; n <= Order; ++n)
   {
      v.push_back(&D);
      Result += (1.0 / factorial(n)) * SumPerms(v);
   }
   return Result;
}

// Version that takes the set of operators as individual arguments
template <typename... Args>
OperatorComponent DysonSeries(OperatorComponent const& D, int Order, Args const&... list)
{
   return DysonSeries(D, Order, {std::ref(list)...});
}

ProductMPO expmpo_first_strict(BasicTriangularMPO const& HamMPO, std::complex<double> Tau)
{
   std::deque<OperatorComponent> Result;

   for (auto const& I : HamMPO)
   {
      // Get the indices for the middle rows and columns.
      std::set<int> MidRows, MidCols;
      for (int i = 1; i < I.Basis1().size()-1; ++i)
         MidRows.insert(i);
      for (int j = 1; j < I.Basis2().size()-1; ++j)
         MidCols.insert(j);

      // Get the Identity, A, B, C and D blocks:
      // (I C D)
      // (0 A B)
      // (0 0 I)
      OperatorComponent Ident = project_columns(project_rows(I, {0}), {0});
      OperatorComponent A = project_columns(project_rows(I, MidRows), MidCols);
      OperatorComponent B = project_columns(project_rows(I, MidRows), {(int) I.Basis2().size()-1});
      OperatorComponent C = project_columns(project_rows(I, {0}), MidCols);
      OperatorComponent D = project_columns(project_rows(I, {0}), {(int) I.Basis2().size()-1});

      // Calculate the four blocks of the optimal first-order evolution MPO.
      OperatorComponent WTopLeft = Ident + Tau*D + 0.5*Tau*Tau*aux_tensor_prod(D, D);
      //OperatorComponent WTopLeft = Ident + Tau * D + 0.5*Tau*Tau * aux_tensor_prod(D, D);
      OperatorComponent WBotLeft = Tau * B + 0.5*Tau*Tau * (aux_tensor_prod(B, D) + aux_tensor_prod(D, B));
      OperatorComponent WTopRight = C + 0.5*Tau * (aux_tensor_prod(C, D) + aux_tensor_prod(D, C));
      OperatorComponent WBotRight = A + 0.5*Tau * (aux_tensor_prod(A, D) + aux_tensor_prod(D, A)
                                                 + aux_tensor_prod(B, C) + aux_tensor_prod(C, B));

      // Construct the evolution MPO.
      OperatorComponent W = tensor_join(
         {{WTopLeft, WTopRight},
          {WBotLeft, WBotRight}}
      );

      Result.push_back(W);
   }

   return ProductMPO(GenericMPO(Result.begin(), Result.end()));
}

ProductMPO expmpo_first(BasicTriangularMPO const& HamMPO, std::complex<double> Tau, int N)
{
   std::deque<OperatorComponent> Result;

   for (auto const& I : HamMPO)
   {
      // Get the indices for the middle rows and columns.
      std::set<int> MidRows, MidCols;
      for (int i = 1; i < I.Basis1().size()-1; ++i)
         MidRows.insert(i);
      for (int j = 1; j < I.Basis2().size()-1; ++j)
         MidCols.insert(j);

      // Get the Identity, A, B, C and D blocks:
      // (I L D)
      // (0 A R)
      // (0 0 I)
      OperatorComponent Ident = project_columns(project_rows(I, {0}), {0});
      OperatorComponent A = project_columns(project_rows(I, MidRows), MidCols);
      OperatorComponent R = project_columns(project_rows(I, MidRows), {(int) I.Basis2().size()-1});
      OperatorComponent L = project_columns(project_rows(I, {0}), MidCols);
      OperatorComponent D = project_columns(project_rows(I, {0}), {(int) I.Basis2().size()-1});

      // Calculate the four blocks of the optimal first-order evolution MPO.
      OperatorComponent WTopLeft = DysonSeries(Tau*D, N, Ident);
      OperatorComponent WBotLeft = Tau*DysonSeries(Tau*D, N, R);
      OperatorComponent WTopRight = DysonSeries(Tau*D, N, L);
      OperatorComponent WBotRight = DysonSeries(Tau*D, N, A) + Tau*DysonSeries(Tau*D, N, {L,R});

      // Construct the evolution MPO.
      OperatorComponent W = tensor_join(
         {{WTopLeft, WTopRight},
          {WBotLeft, WBotRight}}
      );

      Result.push_back(W);
   }

   return ProductMPO(GenericMPO(Result.begin(), Result.end()));
}

ProductMPO expmpo_second(BasicTriangularMPO const& HamMPO, std::complex<double> Tau, int N)
{
   std::deque<OperatorComponent> Result;

   for (auto const& I : HamMPO)
   {
      // Get the indices for the middle rows and columns.
      std::set<int> MidRows, MidCols;
      for (int i = 1; i < I.Basis1().size()-1; ++i)
         MidRows.insert(i);
      for (int j = 1; j < I.Basis2().size()-1; ++j)
         MidCols.insert(j);

      // Get the Identity, A, R, L and D blocks:
      // (I R D)
      // (0 A L)
      // (0 0 I)
      OperatorComponent Ident = project_columns(project_rows(I, {0}), {0});
      OperatorComponent A = project_columns(project_rows(I, MidRows), MidCols);
      OperatorComponent R = project_columns(project_rows(I, MidRows), {(int) I.Basis2().size()-1});
      OperatorComponent L = project_columns(project_rows(I, {0}), MidCols);
      OperatorComponent D = project_columns(project_rows(I, {0}), {(int) I.Basis2().size()-1});

      // Calculate the four blocks of the optimal first-order evolution MPO.
      OperatorComponent W11 = DysonSeries(Tau*D, N, Ident);
      OperatorComponent W12 = DysonSeries(Tau*D, N, L);
      OperatorComponent W13 = DysonSeries(Tau*D, N, L, L);
      OperatorComponent W21 = Tau*DysonSeries(Tau*D, N, R);
      OperatorComponent W22 = DysonSeries(Tau*D, N, A) + Tau*DysonSeries(Tau*D, N, {L,R});
      OperatorComponent W23 = 2*DysonSeries(Tau*D, N, {A,L}) + Tau*DysonSeries(Tau*D, N, {L,L,R});
      OperatorComponent W31 = (0.5*Tau*Tau)*DysonSeries(Tau*D, N, {R,R});
      OperatorComponent W32 = Tau*DysonSeries(Tau*D, N, {A,R}) + (0.5*Tau*Tau)*DysonSeries(Tau*D, N, {L,R,R});
      OperatorComponent W33 = DysonSeries(Tau*D, N, {A,A}) + 2.0*Tau*DysonSeries(Tau*D, N, {A,L,R});

      // Construct the evolution MPO.
      OperatorComponent W = tensor_join(
         {{W11, W12, W13},
          {W21, W22, W23},
          {W31, W32, W33}}
      );

      Result.push_back(W);
   }

   return ProductMPO(GenericMPO(Result.begin(), Result.end()));
}


ProductMPO expmpo(BasicTriangularMPO const& x, ExpMpoScheme Scheme, ExpMpoSchemeParameters p)
{
   if (Scheme == ExpMpoScheme::first)
      return expmpo_first(x, 1.0, p.DysonOrder);
   else if (Scheme == ExpMpoScheme::first_strict)
      return expmpo_first_strict(x);
   else if (Scheme == ExpMpoScheme::second)
      return expmpo_second(x, 1.0, p.DysonOrder);
   PANIC("invalid expmpo scheme");
   __builtin_unreachable();
}
