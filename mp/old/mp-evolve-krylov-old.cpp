// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-evolve-krylov-old.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/simplekrylov.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include <boost/program_options.hpp>
#include "linearalgebra/matrix_utility.h"
#include <iostream>

namespace prog_opt = boost::program_options;

Matrix<std::complex<double> > M;

Matrix<std::complex<double> >
BackSubstitute(Matrix<std::complex<double> > const& M)
{
   Matrix<std::complex<double> > Result(M);
   for (int i = Result.size1()-1; i > 0; --i)
   {
      for (int j = 0; j < i; ++j)
      {
         Result(i,all) -=

void SweepRight(SimpleKrylov& dmrg, bool TwoSite, int MaxStates, double MixFactor)
{
   dmrg.ExpandLeft();
   if (TwoSite) dmrg.ExpandRight();
   dmrg.ConstructKrylovBasis(M);
   //dmrg.MaximizeKrylovVectors();
   //dmrg.OrthogonalizeKrylovVectors();
   dmrg.TruncateLeft(MaxStates, MixFactor);
   std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ")\n";
   // sweep right
   while (dmrg.RightSize() > 1)
   {
      dmrg.ShiftRightAndExpand();
      if (TwoSite) dmrg.ExpandRight();
      dmrg.ConstructKrylovBasis(M);
      //dmrg.MaximizeKrylovVectors();
      //dmrg.OrthogonalizeKrylovVectors();
      dmrg.TruncateLeft(MaxStates, MixFactor);
      std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ")\n";
   }
}

void SweepLeft(SimpleKrylov& dmrg, bool TwoSite, int MaxStates, double MixFactor)
{
   dmrg.ExpandRight();
   if (TwoSite) dmrg.ExpandLeft();
   dmrg.ConstructKrylovBasis(M);
   //dmrg.MaximizeKrylovVectors();
   //dmrg.OrthogonalizeKrylovVectors();
   dmrg.TruncateRight(MaxStates, MixFactor);
   std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ")\n";
   // sweep left
   while (dmrg.LeftSize() > 1)
   {
      dmrg.ShiftLeftAndExpand();
      if (TwoSite) dmrg.ExpandLeft();
      dmrg.ConstructKrylovBasis(M);
      //dmrg.MaximizeKrylovVectors();
      //dmrg.OrthogonalizeKrylovVectors();
      dmrg.TruncateRight(MaxStates, MixFactor);
      std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ")\n";
   }
}

int main(int argc, char** argv)
{
   try
   {
      int NumKrylov = 4;
      int MaxStates = 100;
      double MixFactor = 0.01;
      bool TwoSite = false;
      int NumSweeps = 2;
      std::complex<double> Timestep = 0.1;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options");
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(),
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "wavefunction to evolve (required)")
         ("two-site,2", "modify 2 neighboring sites at once")
         ("krylov,k", prog_opt::value<int>(&NumKrylov), "Number of Krylov vectors")
         ("max-states,m", prog_opt::value<int>(&MaxStates), "Number of states to keep [default 100]")
         ("mix-factor,f", prog_opt::value<double>(&MixFactor), "Mixing coefficient for the density matrix [default 0.01]")
         ("sweeps,s", prog_opt::value<int>(&NumSweeps), "Number of half-sweeps to perform [default 2]")
         ("timestep,t", prog_opt::value<std::complex<double> >(&Timestep), "Time step [default 0.01]")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-evolve-krylov [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout << "Starting Krylov...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

      // Open the wavefunction
      pvalue_ptr<MPWavefunction> P = pheap::OpenPersistent(InputWavefunction, 655360);
      MPWavefunction Psi = *P;
      P = pvalue_ptr<MPWavefunction>();
      // make sure the wavefunction is normalized
      Psi.normalize();

      // Make sure the center matrix is at one edge
      if (Psi.LeftSize() != 1 && Psi.RightSize() != 1)
      {
         TRACE(Psi.LeftSize())(Psi.RightSize());
         std::cout << "The center matrix is not located at an edge.  Rotating..." << std::flush;
         if (Psi.LeftSize() > Psi.RightSize())
         {
            while (Psi.RightSize() > 1)
               Psi.RotateRight();
         }
         else
         {
            while (Psi.LeftSize() > 1)
               Psi.RotateLeft();
         }
         std::cout << "done" << std::endl;
      }

      // Set up the Hamiltonian
      std::string HamString;
      if (vm.count("Hamiltonian") == 1)
      {
         HamString = vm["Hamiltonian"].as<std::string>();
         Psi.Attributes()["Hamiltonian"] = HamString;
      }
      else
      {
         if (Psi.Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: No Hamiltonian specified.\n"
               "Specify a Hamiltonian either with --Hamiltonian parameter, "
               "or as an attribute of the initial wavefunction.\n";
            return 1;
         }
         HamString = Psi.Attributes()["Hamiltonian"].as<std::string>();
      }
      std::cout << "Hamiltonian: " << HamString << std::endl;
      MPOperator Hamiltonian = ParseOperator(HamString);

      // Now we can construct the actual SimpleKrylov object
      std::cout << "Constructing Krylov vectors..." << std::flush;
      TRACE(expectation(Psi, Hamiltonian, Psi));
      SimpleKrylov dmrg(Hamiltonian, std::vector<MPWavefunction>(NumKrylov+1, Psi));
      Psi = MPWavefunction(); // we don't need Psi anymore, it will take up space on disk
      std::cout << "done" << std::endl;

      TwoSite = vm.count("two-site");
      if (TwoSite)
         std::cout << "Optimizing two sites at a time" << std::endl;

      std::cout << "Density matrix mixing coefficient: " << MixFactor << std::endl;
      std::cout << "Number of Krylov vectors: " << NumKrylov << std::endl;
      std::cout << "Number of half-sweeps: " << NumSweeps << std::endl;
      std::cout << "Timestep: " << Timestep << std::endl;


      M = Matrix<std::complex<double> >(NumKrylov+1, NumKrylov+1, 0.0);

         TRACE(expectation(dmrg.Krylov[0], Hamiltonian, dmrg.Krylov[0]));
      dmrg.ConstructKrylovBasis(M);
         TRACE(expectation(dmrg.Krylov[0], Hamiltonian, dmrg.Krylov[0]));
      MPWavefunction OldPsi = dmrg.Krylov.back();

      for (int i = 0; i <= NumKrylov; ++i)
      {
         M(i,i) = 1.0;
      }
      Matrix<std::complex<double> > D = dmrg.OrthogonalityMatrix();
      CholeskyFactorizeLower(D);
      zero_upper_triangular(D);
      InvertLowerTriangular(D);
      M = 0.01 * M + 0.99 * D * M;

      for (int Sweeps = 0; Sweeps < NumSweeps; ++Sweeps)
      {
         TRACE(expectation(dmrg.Krylov[0], Hamiltonian, dmrg.Krylov[0]));
         if (dmrg.LeftSize() == 1)
            SweepRight(dmrg, TwoSite, MaxStates, MixFactor);
         else
            SweepLeft(dmrg, TwoSite, MaxStates, MixFactor);

         dmrg.DebugCheckBasis();

         if (Sweeps % 10 == 0)
         {
            Matrix<std::complex<double> > D = dmrg.OrthogonalityMatrix();
            TRACE(D);
            CholeskyFactorizeLower(D);
            zero_upper_triangular(D);
            TRACE(D);
            InvertLowerTriangular(D);
            TRACE(D)(M);
            M = 0.01 * M + 0.99 * D * M;
            TRACE(M);
         }

         Psi = dmrg.Krylov.back();
         double NormF = norm_frob_sq(Psi) + norm_frob_sq(OldPsi);
         double Overlap = (NormF - 2.0 * norm_frob(overlap(Psi, OldPsi))) / NormF;

         std::cout << "Last Krylov vector difference from last half-sweep = " << Overlap << '\n';
         OldPsi = Psi;

         TRACE(dmrg.OrthogonalityMatrix());
      }

      Matrix<std::complex<double> > subI(NumKrylov, NumKrylov);
      Matrix<std::complex<double> > subH(NumKrylov, NumKrylov);
      for (int i = 0; i < NumKrylov; ++i)
      {
         subI(i,i) = norm_frob_sq(dmrg.Krylov[i]);
         subH(i,i) = expectation(dmrg.Krylov[i], Hamiltonian, dmrg.Krylov[i]);
         for (int j = i+1; j < NumKrylov; ++j)
         {
            subI(i,j) = overlap(dmrg.Krylov[i], dmrg.Krylov[j]);
            subI(j,i) = conj(subI(i,j));
            subH(i,j) = expectation(dmrg.Krylov[i], Hamiltonian, dmrg.Krylov[j]);
            subH(j,i) = conj(subH(i,j));
         }
      }
      TRACE(subI)(subH);
      TRACE(EigenvaluesHermitian(subI));
      Matrix<std::complex<double> > subU = subI;
      Vector<double> eigenI = DiagonalizeHermitian(subU);
      for (std::size_t i = 0; i < eigenI.size(); ++i)
      {
         eigenI[i] = 1.0 / sqrt(eigenI[i]);
      }
      subU = diagonal_matrix(eigenI) * subU;
      Matrix<std::complex<double> > Hp = subU * subH * herm(subU);
      TRACE(Hp);
      TRACE(EigenvaluesHermitian(Hp));
      TRACE(subU * subI * herm(subU));

      std::cout << "Constructing evolved wavefunction..." << std::endl;
      //      dmrg.AdvanceTimestep(MaxStates);
      std::cout << "Finished." << std::endl;
      P = pvalue_ptr<MPWavefunction>(new MPWavefunction(dmrg.Krylov[0]));
      pheap::ShutdownPersistent(P);

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
