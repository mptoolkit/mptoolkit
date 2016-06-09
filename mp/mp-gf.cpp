// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-gf.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/math_const.h"
#include "common/environment.h"

double const Minus1OverPi = -math_const::r_1_pi;
typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc < 8 || argc > 9)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-gf <lattice> <operator> <energy> <frequency> <broadening> <cv> <lv> [resid]\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], CacheSize, true);
   MPOperator Ham = (*System)[argv[2]];
   MPOperator Ident = (*System)["I"];
 
   pvalue_ptr<MPWavefunction> CvPtr = pheap::ImportHeap(argv[6]);
   pvalue_ptr<MPWavefunction> RhsPtr = pheap::ImportHeap(argv[7]);

   double Energy = boost::lexical_cast<double>(argv[3]);
   double Freq = boost::lexical_cast<double>(argv[4]);
   double Broad = boost::lexical_cast<double>(argv[5]);

   MPOperator A =  std::complex<double>(Energy + Freq, -Broad) * Ident - Ham;  // argh, we have a conjugation bug somewhere.  Actually, probably everywhere!  But at least we are consistent...
   MPOperator Part = (Energy + Freq) * Ident - Ham;
   MPOperator A2 = prod(Part, Part, Part.TransformsAs()) + (Broad*Broad) * Ident;
   
   MPWavefunction Cv = *CvPtr;
   MPWavefunction Rhs = *RhsPtr;
   MPWavefunction CvBar = conj(Cv); 
   TRACE(overlap(Cv, CvBar));
   MPWavefunction CvImag = complex(0.0,0.5) * (Cv - CvBar);

   //complex Overlap = overlap(CvImag, Rhs);
   double ImagExpectationA2 = expectation(CvImag, A2, CvImag).real();
   double ExpectationA2 = expectation(Cv, A2, Cv).real();
   
   complex Ovr = overlap(Rhs, Cv);

   double OvrPart = expectation(Rhs, A, Cv).real();
   double ResidNorm = ExpectationA2 + norm_2_sq(Rhs) - 2.0 * OvrPart;

   std::cout.precision(12);
   std::cout << "Correction vector norm^2 * broadening = " << (norm_2_sq(Cv)*Broad) << '\n';
   std::cout << "<lv|cv> = " << Ovr << '\n';
   std::cout << "DDMRG functional = " << ((1.0/Broad)*ImagExpectationA2 + 2*Ovr.imag()) << '\n';
   std::cout << "Residual norm^2 = " << ResidNorm << '\n';
   std::cout << "<cv|A^2|cv> = " << ExpectationA2 << '\n';
   std::cout << "norm_2_sq(Rhs) = " << norm_2_sq(Rhs) << '\n';
   std::cout << "Overlap part = " << OvrPart << '\n';

   if (argc == 9)
   {
      std::cout << "Writing residual to " << argv[8] << '\n';
      pvalue_ptr<MPWavefunction> Resid = new MPWavefunction(prod(A, Cv, Cv.TransformsAs()) - Rhs);
      pheap::ExportHeap(argv[8], Resid);
   }
}
