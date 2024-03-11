// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// junk/lanczos.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
#include "quantumnumbers/u1.h"
#include "pheap/pheap.h"
#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include "matrixproduct/copyright.h"
#include <iostream>

#define FAST_B

int main(int argc, char** argv)
{
   if (argc < 5 || argc > 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-lanczos <lattice> <operator> <psi> <iterations> [<max-states>]\n";
      return 1;
   }

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[3], 655360);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string Operator = argv[2];
   int Iterations = boost::lexical_cast<int>(argv[4]);
   int MaxStates = 100000;
   if (argc == 6) MaxStates = boost::lexical_cast<int>(argv[5]);

   int HFMaxStates = MaxStates*4;
   int FMaxStates = MaxStates;

   MPOperator H = (*System)[Operator];

   std::vector<MPWavefunction> f(Iterations+1);
   std::vector<double>         a(Iterations+1);
   std::vector<double>         b(Iterations+1);
   std::vector<double>         NormSq(Iterations+1);

   f[0] = *Psi;

   double PsiNorm = norm_2(*Psi);

   std::cout.precision(14);

   // do n=0 iteration
   int n = 0;
   MPWavefunction Hfn = prod(H, f[n], f[n].TransformsAs(), MaxStates);
   NormSq[n] = norm_2_sq(f[n]);
   a[n] = overlap(Hfn, f[n]) / NormSq[n];
   b[n] = 0;
   double Overlap = 1;
   double Energy = a[0];
   std::cout << "   n                    a                  b^2               Energy     Norm^2    Overlap\n";

      std::cout << std::setw(4) << n << " "
                << std::setprecision(14)
                << std::setw(20) << a[n] << " " << std::setw(20) << b[n]
                << " " << std::setw(20) << Energy
                << std::setprecision(4)
                << " " << std::setw(10) << std::sqrt(NormSq[n]) << " " << std::setw(10)
                << Overlap << std::endl;

   ++n; // n=1 iteration
   f[n] = sum(Hfn, -a[n-1] * f[n-1], FMaxStates);
   Hfn = prod(H, f[n], f[n].TransformsAs(), HFMaxStates);
   NormSq[n] = norm_2_sq(f[n]);
   a[n] = overlap(Hfn, f[n]) / NormSq[n];
#if defined(FAST_B)
   b[n] = NormSq[n] / NormSq[n-1];
#else
   b[n] = overlap(f[n-1], Hfn) / NormSq[n-1];
#endif
   Overlap = overlap(*Psi, f[n]) / (PsiNorm * std::sqrt(NormSq[n]));

   LinearAlgebra::Matrix<double> M(Iterations+1, Iterations+1, 0.0);
   M(0,0) = a[0];
   M(1,1) = a[1];
   M(1,0) = M(0,1) = std::sqrt(b[1]);
   Energy = LinearAlgebra::EigenvaluesSymmetric(M)[0];

      std::cout << std::setw(4) << n << " "
                << std::setprecision(14)
                << std::setw(20) << a[n] << " " << std::setw(20) << b[n]
                << " " << std::setw(20) << Energy
                << std::setprecision(4)
                << " " << std::setw(10) << std::sqrt(NormSq[n]) << " " << std::setw(10)
                << Overlap << std::endl;

   ++n;
   double Now = ProcControl::GetCPUTime();
   for ( ; n <= Iterations; ++n)
   {
#if 1
      std::list<MPWavefunction> SumList;
      SumList.push_back(-a[n-1] * f[n-1]);
      SumList.push_back(-b[n-1] * f[n-2]);
      SumList.push_back(Hfn);
      f[n] = accumulate(SumList.begin(), SumList.end(), FMaxStates);
#else
      f[n] = Hfn - a[n-1] * f[n-1];
      f[n] = sum(f[n], -b[n-1] * f[n-2], FMaxStates);
#endif
      Hfn = prod(H, f[n], f[n].TransformsAs(), HFMaxStates);
      NormSq[n] = norm_2_sq(f[n]);
      a[n] = overlap(Hfn, f[n]) / NormSq[n];
#if defined(FAST_B)
      b[n] = NormSq[n] / NormSq[n-1];
#else
      b[n] = overlap(f[n-1], Hfn) / NormSq[n-1];
#endif
      Overlap = overlap(*Psi, f[n]) / (PsiNorm * std::sqrt(NormSq[n]));

      M(n,n) = a[n];
      M(n-1,n) = M(n,n-1) = std::sqrt(b[n]);

      Energy = LinearAlgebra::EigenvaluesSymmetric(M)[0];

      std::cout << std::setw(4) << n << " "
                << std::setprecision(14)
                << std::setw(20) << a[n] << " " << std::setw(20) << b[n]
                << " " << std::setw(20) << Energy
                << std::setprecision(4)
                << " " << std::setw(10) << std::sqrt(NormSq[n]) << " " << std::setw(10)
                << Overlap << std::endl;
      double Next = ProcControl::GetCPUTime();
      std::cout << "CPU time for sweep = " << Next-Now << std::endl;
      Now = Next;
   }

   LinearAlgebra::Vector<double> EValues;
   DiagonalizeSymmetric(M, EValues);
   std::cout.precision(14);
   std::cout << "Lanczos eigenvalues are: " << EValues << std::endl;

   // calculate the ground-state Lanczos vector
#if 1
   for (int i = 0; i < EValues.size(); ++i)
   {
      f[i] *= M(0,i) / std::sqrt(NormSq[i]);
   }
   *Psi.mutate() = accumulate(f.begin(), f.begin()+EValues.size(), MaxStates);
#else
   *Psi.mutate() = (M(0,0) / std::sqrt(NormSq[0])) * f[0];
   for (int i = 1; i < int(EValues.size())-1; ++i)
   {
      *Psi.mutate() = sum(*Psi, (M(0,i) / std::sqrt(NormSq[i])) * f[i], MaxStates * 2);
      std::cout << "Working..." << std::endl;
   }
   if (EValues.size() > 1)
   {
      int i = EValues.size()-1;
      *Psi.mutate() = sum(*Psi, (M(0,i) / std::sqrt(NormSq[i])) * f[i], MaxStates);
   }
#endif

   Psi.mutate()->normalize();

   pheap::ShutdownPersistent(Psi);
}
