// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// benchmark/bench-tripleprod.cpp
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

#include "tensor/tensor.h"
#include "quantumnumbers/su2.h"
#include "linearalgebra/matrix_utility.h"
#include "common/proccontrol.h"

using namespace Tensor;

typedef std::complex<double> complex;

typedef double mytype;

typedef IrredTensor<LinearAlgebra::Matrix<mytype>, 
                    VectorBasis, 
                    VectorBasis> MatrixOperator;

double MakeTest(MatrixOperator const& L, MatrixOperator const& M, MatrixOperator const& N, int Repeat)
{
   double Start = ProcControl::GetCPUTime();
   for (int i = 0; i < Repeat; ++i)
   {
      MatrixOperator x = triple_prod(herm(L), M, N);
   }
   return ProcControl::GetCPUTime() - Start;
}

double DoTest(MatrixOperator const& L, MatrixOperator const& M, MatrixOperator const& N)
{
   MakeTest(L, M, N, 1);

   int Repeat = 1;

   double Time = 0;
   while (Time < 2)
   {
      Repeat *= 2;
      Time = MakeTest(L, M, N, Repeat);
   }

   return Time / Repeat;
}

double Test(int MatrixSize, int NumQuantum)
{
   SymmetryList SL("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(SL);

   VectorBasis b(SL);
   for (int i = 0; i < NumQuantum; ++i)
   {
      b.push_back(QN(i), MatrixSize);
   }

   MatrixOperator L(b, b, QN(0)), M(b, b, QN(0)), N(b, b, QN(0));
   for (unsigned i = 0; i < b.size(); ++i)
   {
      L(i,i) = LinearAlgebra::random_matrix<mytype>(MatrixSize, MatrixSize);
      M(i,i) = LinearAlgebra::random_matrix<mytype>(MatrixSize, MatrixSize);
      N(i,i) = LinearAlgebra::random_matrix<mytype>(MatrixSize, MatrixSize);
   }

   return DoTest(L, M, N);
}

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      std::cerr << "usage: benchtripleprod <matrix-size> <num-matrices>\n";
      return 1;
   }

   int MSize = boost::lexical_cast<int>(argv[1]);
   int NQuantum = boost::lexical_cast<int>(argv[2]);
   double Time = Test(MSize, NQuantum);
   double Flops = 2.0 * NQuantum * MSize*MSize*MSize;
   std::cout << "Time per iteration: " << Time 
             << "\nNumber of operations: " << Flops
             << "\nMFlop/s: " << (Flops / (Time * 1e6)) << '\n';
}
