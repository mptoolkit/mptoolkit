// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/aggregator.h
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

#if !defined(AGGREGATOR_H_KJY897689HUIH43Y8P9)
#define AGGREGATOR_H_KJY897689HUIH43Y8P9

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"

class Aggregator
{
   public:
      typedef MPStateComponent Component;
      typedef Component::OperatorType OperatorType;

      Aggregator(std::vector<MPWavefunction> Psi,
                 std::vector<MPOperator> Op,
                 bool ShowStates_,
                 int MaxStates,
                 double MinTrunc,
                 int Location = 0,
                 std::vector<double> Weights = std::vector<double>());

      OperatorType const& Center(int n) const
      {
         return Aggregate_center[n];
      }

      MPStateComponent const& OpLeft(int n) const { return Op_left[n]; }
      MPStateComponent const& OpRight(int n) const { return Op_right[n]; }
      SimpleOperator const& OpCenter(int n) const { return Op_center[n]; }

      std::complex<double> Expectation(MatrixOperator const& Psi1, int n) const;

      std::complex<double> Expectation(MatrixOperator const& Psi1, int n,
                                       MatrixOperator const& Psi2) const;

      MatrixOperator Apply(int n, MatrixOperator const& Psi) const;

   private:
      void ConstructLeft(std::vector<MPWavefunction> const& Psi,
                         std::vector<MPOperator> const& Op,
                         std::vector<OperatorType>& LeftMap,
                         std::vector<double> const& Weights,
                         int MaxStates,
                         double MinTrunc);

      void ConstructRight(std::vector<MPWavefunction> const& Psi,
                          std::vector<MPOperator> const& Op,
                          std::vector<OperatorType>& RightMap,
                          std::vector<double> const& Weights,
                          int MaxStates,
                          double MinTrunc);

      void RotateLeft(std::vector<MPWavefunction>& Psi, std::vector<MPOperator>& Op);
      void RotateRight(std::vector<MPWavefunction>& Psi, std::vector<MPOperator>& Op);

      bool ShowStates; // if true, show some info to std::cerr

      std::vector<OperatorType> Aggregate_center;
      std::vector<MPStateComponent> Op_left, Op_right;
      std::vector<SimpleOperator> Op_center;
      MPStateComponent Aggregate_left, Aggregate_right;

      QuantumNumber Ident;
};

#endif
