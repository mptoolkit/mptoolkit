// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-min-resid.cpp
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

#include "mp/copyright.h"
#include "linearalgebra/eigen.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_opt_accum.h"
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>

namespace prog_opt = boost::program_options;

using namespace LinearAlgebra;

typedef std::complex<double> complex;

double const minus1pi = -1.0 / math_const::pi;

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

void DoMinimization(double w, double E,  double Eta, Matrix<complex> const& Ident,
                    Matrix<complex> const& H, Matrix<complex> const& H2,
                    Matrix<complex> const& RealH2, Matrix<complex> const& Lv,
                    bool ShowDOS, bool SolveConj)
{
   // Assemble the left-hand-side matrix (w+E-H)^2 + \eta^2,
   Matrix<complex> LHS = ((w+E)*(w+E) + Eta*Eta)*Ident + H2 - 2.0*(w+E)*H;

   // Assemble the right-hand-side matrix w+E-H - i\eta
   Matrix<complex> RHS = (w+E+complex(0.0, -Eta))*Ident - H;

   // the right hand side vector = RHS * Lv
   Matrix<complex> RHS_vec = RHS * Lv;

   // Linear solver for the correction vector
   Matrix<complex> Cv = LinearSolveHPD(LHS, RHS_vec);

   // residual norm
   double r2 = inner_prod(Cv, Matrix<complex>((LHS-H2+RealH2)*Cv)).real() 
      + inner_prod(Lv, Matrix<complex>(Ident*Lv)).real()
      - 2.0*inner_prod(Cv, Matrix<complex>(RHS*Lv)).real();

   // spectral function - and we have a conjugation bug somewhere
   complex G;
   if (SolveConj)
   {
      // solve the conjugate correction vector
      Matrix<complex> RHS_bar = (w+E+complex(0.0, Eta))*Ident - H;
      Matrix<complex> RHS_bar_vec = RHS_bar * Lv;
      Matrix<complex> Cv_bar = LinearSolveHPD(LHS, RHS_bar_vec);

      Matrix<complex> Cv_im = 0.5 * (Cv - Cv_bar);
      //      TRACE(Cv_im);
      G = complex(0.0, (-1.0/Eta)*(inner_prod(Cv_im, Matrix<complex>(LHS*Cv_im)).real()
                                   + 2.0*Eta*inner_prod(Cv_im,  Matrix<complex>(Ident*Lv)).real()));
      //      TRACE(inner_prod(Cv_im,  Matrix<complex>(Ident*Lv)).real());
      //      TRACE(G);
   }
   else
   {
      G = conj(inner_prod(Cv, Matrix<complex>(Ident*Lv)));
   }

   if (ShowDOS)
      std::cout << std::setw(20) << w
                << "    " << std::setw(20) << (minus1pi*G.imag())
                << "    " << std::setw(20) << r2
                << '\n';
   else
      std::cout << std::setw(20) << w
                << "    " << std::setw(20) << G.real() 
                << "    " << std::setw(20) << G.imag() 
                << "    " << std::setw(20) << r2
                << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      std::string IdentFileStr, HFileStr, H2FileStr;
      double Eta;
      double E;
      double wMin;
      double wMax;
      double wInterval;
      int wCount = 21;
      int Verbose = 0;
      bool UseProjectedH2 = false;
      bool UseOrthoBasis = false;
      bool ShowDOS = false;
      bool SolveConj = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("ident,i", prog_opt::value(&IdentFileStr),
          "file for the matrix elements of the identity operator [required]")
         ("ham,h", prog_opt::value(&HFileStr),
          "file for the matrix elements of H [required]")
         ("ham2,2", prog_opt::value(&H2FileStr),
          "file for the matrix elements of H^2")
         ("Broadening,B", prog_opt::value(&Eta),
          "Broadening")
         ("GroundstateEnergy,G", prog_opt::value(&E),
          "Groundstate energy of the Hamiltonian")
         ("projected", prog_opt::bool_switch(&UseProjectedH2),
          "use the square of the projected Hamiltonian in the minimization, but the "
          "real H^2 for the residual calculation")
         ("orthogonal", prog_opt::bool_switch(&UseOrthoBasis),
          "Orthogonalize the basis before the minimization [probably irrelevant]")
         ("dos,d", prog_opt::bool_switch(&ShowDOS),
          "show the normalized density of states, instead of the real and imag parts "
          "of the Green's function")
         ("min-freq,m", prog_opt::value(&wMin),
          "Minimum frequency to scan [required]")
         ("max-freq,x", prog_opt::value(&wMax),
          "Maximum frequency to scan [required]")
         ("count,c", prog_opt::value(&wCount),
          "Number of frequencies to calculate [default 21]")
         ("interval", prog_opt::value(&wInterval),
          "Alternative to --count, calculate the spectral function at"
          " frequencies separated by this interval")
         ("conjugate", prog_opt::bool_switch(&SolveConj),
          "Solve for the conjugate correction vector to get the DDMRG"
          " functional for the spectral function"
          " (apparantly has no effect)")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "extra debug output")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") || vm.count("ident") == 0 || vm.count("ham") == 0
          || vm.count("Broadening") == 0 || vm.count("GroundstateEnergy") == 0
          || vm.count("min-freq") == 0 || vm.count("max-freq") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-min-resid [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));
      
      if (Verbose >= 1) 
         std::cerr << "mp-min-resid: Loading matrix elements of identity operator from file " 
                   << IdentFileStr << '\n';
      Matrix<complex> Ident = ReadHermitianMatrix(IdentFileStr);
      if (Verbose >= 2)
         std::cerr << "I = " << Ident << '\n';

      Vector<double> IdentEigen = EigenvaluesHermitian(Ident);
      if (min(IdentEigen) < std::numeric_limits<double>::epsilon() * 10)
      {
         PANIC("Identity operator is not positive definite")(IdentEigen);
      }           

      if (Verbose >= 1) 
         std::cerr << "mp-min-resid: Loading matrix elements of Hamiltonian operator from file " 
                   << HFileStr << '\n';
      Matrix<complex> H = ReadHermitianMatrix(HFileStr);
      if (Verbose >= 2)
         std::cerr << "H = " << H << '\n';


      // when we calculate the residual, use the proper matrix for H^2 if possible,
      // so we can play around with H2 in the solver and see what effect it has
      Matrix<complex> H2, RealH2;
      if (vm.count("ham2"))
      {
         if (Verbose >= 1) 
            std::cerr << "mp-min-resid: Loading matrix elements of Hamiltonian"
               " squared operator from file " << H2FileStr << '\n';
         RealH2 = ReadHermitianMatrix(H2FileStr);
         if (Verbose >= 2)
            std::cerr << "H^2 = " << RealH2 << '\n';
      }
      else
      {
         std::cerr << "mp-min-resid: warning: no matrix elements for H^2 supplied.  "
            "The residual calculation will not be correct.\n";
      }

      if (UseProjectedH2 || vm.count("ham2") == 0)
      {
         if (Verbose >= 1)
            std::cerr << "mp-min-resid: using projected H for the matrix elements of H^2.\n";

         // Set H^2 = H*H in projected basis
         Matrix<complex> IdentInv = Ident;
         InvertHPD(IdentInv);
         H2 = H*IdentInv*H;
         if (vm.count("ham2") == 0)
            RealH2 = H2;
      }
      else
         H2 = RealH2;

      // The lanczos vector is (1,0,0,...) in this basis.
      // Actually, better make it a matrix 
      Matrix<complex> Lv(size1(H), 1, 0.0); Lv(0,0) = 1.0;

      if (UseOrthoBasis)
      {
         // find the transform for the orthonormal basis
         if (Verbose >= 1)
            std::cerr << "mp-min-resid: transforming to orthonormal basis.\n";
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

         // switch everything to the orthonormal basis
         Ident = herm(Xinv) * Ident * Xinv;
         H = herm(Xinv) * H * Xinv;
         H2 = herm(Xinv) * H2 * Xinv;
         RealH2 = herm(Xinv) * RealH2 * Xinv;
         Lv = X * Lv;
      }

      if (ShowDOS)
         std::cout << "#Frequency              #DOS                    #resid\n";
      else
         std::cout << "#Frequency              #Real                   "
            "#Imag                   #resid\n";

      if (vm.count("interval"))
      {
         for (double w = wMin; w < wMax+wInterval/10; w += wInterval)
         {
            DoMinimization(w,E,Eta,Ident,H,H2,RealH2,Lv,ShowDOS,SolveConj);
         }
      }
      else
      {
         for (int i = 0; i < wCount; ++i)
         {
            double w = wMin + (wMax - wMin) * i / (wCount-1);
            DoMinimization(w,E,Eta,Ident,H,H2,RealH2,Lv,ShowDOS,SolveConj);
         }
      }
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
