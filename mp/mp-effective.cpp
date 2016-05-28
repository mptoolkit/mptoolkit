// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-effective.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include "linearalgebra/eigen.h"
#include <iostream>
#include "common/environment.h"
#include "interface/operator-parser.h"

using namespace LinearAlgebra;

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc < 4)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-effective <wavefunction> <lattice> <op1> <op2> <op3> ... \n";
      return 1;
   }

   int Verbose = 2;

   std::string WavefuncFile = argv[1];
   std::string LatticeFile = argv[2];
   std::vector<std::string> Operators(&argv[3], &argv[argc]);

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::OpenPersistent(WavefuncFile, CacheSize, true);
   MPWavefunction Psi = *Psi1;

   pvalue_ptr<OperatorList> System = pheap::ImportHeap(LatticeFile);

   Vector<double> Hdiag(Operators.size());
   for (unsigned i = 0; i < Operators.size(); ++i)
   {
      CHECK(System->HasOperator(Operators[i]))("Operator does not exist!")(Operators[i]);
      complex x = expectation(Psi, (*System)[Operators[i]], Psi);
      CHECK(std::abs(x.imag()) < 1E-10)("Operator is not Hermitian")(Operators[i]);
      Hdiag[i] = x.real();
   }

   TRACE(Hdiag);

   Matrix<double> Hmat(Operators.size(), Operators.size());
   for (unsigned i = 0; i < Operators.size(); ++i)
   {
      if (Verbose >= 2)
         std::cout << "Calculating row " << i << " of " << Operators.size() << std::endl;
      for (unsigned j = i; j < Operators.size(); ++j)
      {
	 MPOperator Op = (*System)[Operators[i]] * (*System)[Operators[j]];
	 Hmat(i,j) = expectation(Psi, Op, Psi).real();
	 Hmat(j,i) = Hmat(i,j);
      }
   }

   //TRACE(Hmat);

   Matrix<double> U = Hmat;
   Vector<double> Eigen = DiagonalizeHermitian(U);

   //TRACE(U);

   //TRACE(Eigen);
   
   Vector<double> HdiagInDiagBasis(Operators.size(), 0.0);
   for (unsigned i = 0; i < Operators.size(); ++i)
   {
      for (unsigned j = 0; j < Operators.size(); ++j)
      {
	 HdiagInDiagBasis[i] += U(i,j) * Hdiag[j];
      }
   }

   Vector<double> Alpha(Operators.size(), 0.0);
   double Resid = 1.0;
   for (unsigned i = 0; i < Operators.size(); ++i)
   {
      Alpha[i] = HdiagInDiagBasis[i] / Eigen[i];
      Resid -= HdiagInDiagBasis[i] * Alpha[i];
   }

   //TRACE(Alpha)(Resid);

   // transform Alpha back to the normal basis.
   Vector<double> AlphaOrig(Operators.size(), 0.0);
   for (unsigned i = 0; i < Operators.size(); ++i)
   {
      for (unsigned j = 0; j < Operators.size(); ++j)
      {
	 AlphaOrig[i] += U(j,i) * Alpha[j];
      }
   }

   TRACE(AlphaOrig);
   TRACE(AlphaOrig*((1.0/AlphaOrig[0])));

   pheap::Shutdown();
}
