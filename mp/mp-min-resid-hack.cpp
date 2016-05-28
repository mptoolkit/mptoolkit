// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-min-resid-hack.cpp
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

#include "linearalgebra/eigen.h"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>

using namespace LinearAlgebra;

typedef std::complex<double> complex;

// read a matrix from a list of (row, col, value)
Matrix<complex> ReadHermitianMatrix(std::string const& FName)
{
   typedef std::map<std::pair<int, int>, complex> MatMapType;
   int i,j;
   complex z;
   int MaxRow = 0;
   int MaxCol = 0;
   MatMapType Mat;
   std::ifstream F(FName.c_str());
   while (F >> i >> j >> z)
   {
      MaxRow = std::max(MaxRow, i);
      MaxCol = std::max(MaxCol, j);
      Mat[std::make_pair(i,j)] = z;
   }

   Matrix<complex> Result(MaxRow+1, MaxCol+1, 0.0);
   for (MatMapType::const_iterator I = Mat.begin(); I != Mat.end(); ++I)
   {
      Result(I->first.first, I->first.second) = I->second;
      Result(I->first.second, I->first.first) = conj(I->second);
   }
   return Result;
}

Matrix<complex> ReadHermitianMatrixArgh(std::string const& FName)
{
   typedef std::map<std::pair<int, int>, complex> MatMapType;
   int i,j;
   complex z;
   int MaxRow = 0;
   int MaxCol = 0;
   MatMapType Mat;
   std::ifstream F(FName.c_str());
   while (F >> i >> j >> z)
   {
      i=i-1; j=j-1;
      MaxRow = std::max(MaxRow, i);
      MaxCol = std::max(MaxCol, j);
      Mat[std::make_pair(i,j)] = z;
   }

   Matrix<complex> Result(MaxRow+1, MaxCol+1, 0.0);
   for (MatMapType::const_iterator I = Mat.begin(); I != Mat.end(); ++I)
   {
      Result(I->first.first, I->first.second) = I->second;
      Result(I->first.second, I->first.first) = conj(I->second);
   }
   return Result;
}

int main(int argc, char** argv)
{
   if (argc != 9)
   {
      std::cerr << "usage: mp-min-resid <identity-mat-elements> <H-mat-elements> <H^2-mat-elements>"
         " <eta> <GroundstateEnergy> <w-min> <w-max> <w-count>\n";
      return 1;
   }

   std::string const IdentFileStr = argv[1];
   std::string const HFileStr = argv[2];
   std::string const H2FileStr = argv[3];

   double const Eta = boost::lexical_cast<double>(argv[4]);
   double const E = boost::lexical_cast<double>(argv[5]);
   double const wMin = boost::lexical_cast<double>(argv[6]);
   double const wMax = boost::lexical_cast<double>(argv[7]);
   int const wCount = boost::lexical_cast<int>(argv[8]);

   TRACE(Eta)(E);

   Matrix<complex> Ident = ReadHermitianMatrix(IdentFileStr);
   Matrix<complex> H = ReadHermitianMatrix(HFileStr);
   Matrix<complex> H2 = ReadHermitianMatrixArgh(H2FileStr);

   // The lanczos vector is (1,0,0,...) in this basis.
   // Actually, better make it a matrix 
   Matrix<complex> Lv(size1(H), 1, 0.0); Lv(0,0) = 1.0;

   // when we calculate the residual, use the proper matrix for H^2,
   // so we can play around with H2 in the solver and see what effect it
   Matrix<complex> RealH2 = H2;

   // Set H^2 = H*H in projected basis (for a test!)
#if 0
   Matrix<complex> IdentInv = Ident;
   InvertHPD(IdentInv);
   H2 = H*IdentInv*H;
#endif

   // transform into the orthonormal basis
#if 1
   Matrix<complex> X = Ident;
   CholeskyFactorizeUpper(X); // Ident = herm(X) * X
   // zero out the lower-triangular part
   for (unsigned i = 0; i < size1(X); ++i)
   {
      for (unsigned j = 0; j < i; ++j)
      {
         X(i,j) = 0.0;
      }
   }
   Matrix<complex> Xinv = X;
   InvertUpperTriangular(Xinv);

   // switch to ON basis
   Ident = herm(Xinv) * Ident * Xinv;
   H = herm(Xinv) * H * Xinv;
   H2 = herm(Xinv) * H2 * Xinv;
   RealH2 = herm(Xinv) * RealH2 * Xinv;
   Lv = X * Lv;
#endif

   std::cout.precision(14);

   //TRACE(norm_frob_sq(RealH2 - H*H));
   //TRACE(norm_frob_sq(RealH2 - H2));

   for (int i = 0; i < wCount; ++i)
   {
      double w = wMin + ((wMax - wMin) / wCount) * i;

      // Assemble the left-hand-side matrix (w+E-H)^2 + \eta^2,
      Matrix<complex> LHS = ((w+E)*(w+E) + Eta*Eta)*Ident + H2 - 2.0*(w+E)*H;

      //TRACE(EigenvaluesHermitian(LHS));

      // Assemble the right-hand-side matrix w+E-H - i\eta
      Matrix<complex> RHS = (w+E+complex(0.0, -Eta))*Ident - H;

      Matrix<complex> RHS_vec = RHS * Lv;

      // Linear solver for the correction vector
      Matrix<complex> Cv = LinearSolveHPD(LHS, RHS_vec);

      // residual norm
      double r2 = inner_prod(Cv, Matrix<complex>((LHS-H2+RealH2)*Cv)).real() 
         + inner_prod(Lv, Matrix<complex>(Ident*Lv)).real()
         - 2.0*inner_prod(Cv, Matrix<complex>(RHS*Lv)).real();

      // spectral function - and we have a conjugation bug somewhere
      complex G = conj(inner_prod(Cv, Matrix<complex>(Ident*Lv)));

      std::cout << std::setw(20) << w
                << "    " << std::setw(20) << G.real() 
                << "    " << std::setw(20) << G.imag() 
                << "    " << std::setw(20) << r2
                << std::endl;
   }
}
