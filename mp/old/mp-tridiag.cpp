// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-tridiag.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/dmrg.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include <boost/program_options.hpp>
#include "common/prog_opt_accum.h"
#include <iostream>

namespace prog_opt = boost::program_options;

bool ShowStates = false; // for verbose output
bool ShowTitles = false;
bool ShowStoppingReason = false;

double MinTrunc = std::numeric_limits<double>::epsilon()*8;
int MaxStates = 100000;
double Eps = 1E-5;

class Aggregator
{
   public:
      typedef MPStateComponent Component;
      typedef Component::OperatorType OperatorType;

      Aggregator(std::vector<CenterWavefunction> const& Psi_, int RightLanczos_,
                 int LeftLanczos_, SplitOperator const& H_, int Location_);

      void Tridiagonalize(int NumIter, double Threshold);

   private:
      void RotateRight();
      void ConstructLeft();

      void RotateLeft();
      void ConstructRight();

      std::vector<CenterWavefunction> Psi;
      int RightLanczos, LeftLanczos;
      SplitOperator H;

      std::vector<OperatorType> LeftMap, RightMap;
      std::vector<double> Weights; // density matrix weights for each state

      CenterWavefunction Result;
      std::vector<OperatorType> Center;
      MPStateComponent H_left, H_right;

      QuantumNumber Ident;
};

void Aggregator::ConstructLeft()
{
   // Construct the mapping from the vectors to the result space
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      LeftMap[i] = operator_prod(herm(Result.Left()), LeftMap[i], Psi[i].Left());
   }

   // Construct the density matrix
   OperatorType Rho;
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      OperatorType ThisPart = triple_prod(LeftMap[i], scalar_prod(Psi[i].Center(), herm(Psi[i].Center())),
                                          herm(LeftMap[i]));
      Rho += (Weights[i] / norm_frob(ThisPart)) * ThisPart;
   }
   //Rho *= 1.0 / norm_frob(Rho);

   // Form the density matrix
   DensityMatrix<OperatorType> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                                                                           0,
                                                                                           MaxStates,
                                                                                           MinTrunc,
                                                                                           Info));
   if (ShowStates)
      std::cerr << "left density matrix at partition (" << Psi[0].LeftSize() << "," << Psi[0].RightSize()
                << "), states=" << Info.KeptStates() << ", trunc=" << Info.TruncationError() << '\n';

   // Truncate
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      LeftMap[i] = prod(U, LeftMap[i], Ident);
   }
   Result.Left() = prod(Result.Left(), herm(U));

   // construct matrix elements of H
   H_left = operator_prod(herm(H.Left()), herm(Result.Left()), H_left, Result.Left());
}

void Aggregator::ConstructRight()
{
   // Construct the mapping from the vectors to the result space
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      RightMap[i] = operator_prod(Result.Right(), RightMap[i], herm(Psi[i].Right()));
   }

   // Construct the density matrix
   OperatorType Rho;
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      OperatorType ThisPart = triple_prod(RightMap[i], scalar_prod(herm(Psi[i].Center()), Psi[i].Center()),
                                          herm(RightMap[i]));
      Rho += (Weights[i] / norm_frob(ThisPart)) * ThisPart;
   }
   //Rho *= 1.0 / norm_frob(Rho);

   // Form the density matrix
   DensityMatrix<OperatorType> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                                                                           0,
                                                                                           MaxStates,
                                                                                           MinTrunc,
                                                                                           Info));
   if (ShowStates)
      std::cerr << "right density matrix at partition (" << Psi[0].LeftSize() << "," << Psi[0].RightSize()
                << "), states=" << Info.KeptStates() << ", trunc=" << Info.TruncationError() << '\n';

   // Truncate
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      RightMap[i] = prod(U, RightMap[i], Ident);
   }
   Result.Right() = prod(U, Result.Right());

   // construct matrix elements of H
   H_right = operator_prod(H.Right(), Result.Right(), H_right, herm(Result.Right()));
}

void Aggregator::RotateLeft()
{
   H.RotateLeft();
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      Psi[i].RotateLeft();
   }
   Result.PushRight(Component::ConstructFullBasis1(Psi[0].Right().SiteBasis(),
                                                   Result.Right().Basis1()));
}

void Aggregator::RotateRight()
{
   H.RotateRight();
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      Psi[i].RotateRight();
   }
   Result.PushLeft(Component::ConstructFullBasis2(Result.Left().Basis2(),
                                                  Psi[0].Left().SiteBasis()));
}

Aggregator::Aggregator(std::vector<CenterWavefunction> const& Psi_,  int RightLanczos_,
                       int LeftLanczos_, SplitOperator const& H_, int Location_)
   : Psi(Psi_), RightLanczos(RightLanczos_), LeftLanczos(LeftLanczos_),
     H(H_), Weights(Psi_.size(), 1.0), Ident(Psi_[0].GetSymmetryList())
{
   // Rotate the wavefunctions to left most position
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      while (Psi[i].LeftSize() > 1)
         Psi[i].RotateLeft();
      while (H.LeftSize() > 1)
         H.RotateLeft();
   }

   // Initialize the result matrices
   OperatorType LeftVac = OperatorType::make_identity(Psi[0].LeftVacuumBasis());
   LeftMap = std::vector<OperatorType>(Psi.size(), LeftVac);

   H_left = make_vacuum_state(Psi[0].LookupLeft(0).Basis1()[0]);
   Result.PushLeft(Component::ConstructFullBasis2(Psi[0].Left().Basis1(),
                                                  Psi[0].Left().SiteBasis()));

   this->ConstructLeft();

   while (Psi[0].LeftSize() < Location_)
   {
      this->RotateRight();
      this->ConstructLeft();
   }

   // now shift all the way to the right and construct the right operators
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      while (Psi[i].RightSize() > 1)
         Psi[i].RotateRight();
      while (H.RightSize() > 1)
         H.RotateRight();
   }

   OperatorType RightVac = OperatorType::make_identity(Psi[0].RightVacuumBasis());
   RightMap = std::vector<OperatorType>(Psi.size(), RightVac);

   H_right = make_vacuum_state(Psi[0].GetSymmetryList());
   Result.PushRight(Component::ConstructFullBasis1(Psi[0].Right().SiteBasis(),
                                                    Psi[0].Right().Basis2()));

   this->ConstructRight();

   while (Psi[0].LeftSize() > Location_)
   {
      this->RotateLeft();
      this->ConstructRight();
   }

   // Now construct the center matrices
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      Center.push_back(triple_prod(LeftMap[i], Psi[i].Center(), herm(RightMap[i])));
   }
}

void Aggregator::Tridiagonalize(int MaxIter, double Threshold)
{
   if (ShowTitles)
      std::cerr << "  n                   Alpha                    Beta              LeftReal                  LeftImag\n";

   std::vector<std::complex<double> > LeftVector; // left lanczos vector in the Lanczos basis
   MatrixOperator Initial = Center[RightLanczos]; //triple_prod(LeftMap[0], RightLanczos.Center(), herm(RightMap[0]));
   MatrixOperator LeftL = Center[LeftLanczos];
   MatrixOperator f = Initial;
   std::vector<std::complex<double> > Alpha;
   std::vector<double> Beta;
   Alpha.push_back(0.0);
   Beta.push_back(norm_frob(f));
   f *= 1.0 / Beta.back();
   Initial = f; // we want Initial to be normalized properly
   LeftVector.push_back(inner_prod(LeftL, f));
   MatrixOperator fLast = MatrixOperator();
   std::cout << std::setw(3) << (Alpha.size()-1) << "   "
             << std::setw(20) << Alpha.back().real() << "   "
             << std::setw(20) << Beta.back() << "  "
             << std::setw(20) << LeftVector.back().real() << "  "
             << std::setw(20) << LeftVector.back().imag() << std::endl;

   MatrixOperator Hf = operator_prod(conj(H.Center()), H_left, f, herm(H_right));
   Alpha.push_back(inner_prod(f, Hf));
   Hf -= Alpha.back() * f + Beta.back() * fLast;
   fLast = f;
   f = Hf;
   Beta.push_back(norm_frob(f));
   f *= 1.0 / Beta.back();
   LeftVector.push_back(inner_prod(LeftL, f));
   std::cout << std::setw(3) << (Alpha.size()-1) << "   "
             << std::setw(20) << Alpha.back().real() << "   "
             << std::setw(20) << Beta.back() << "  "
             << std::setw(20) << LeftVector.back().real() << "  "
             << std::setw(20) << LeftVector.back().imag() << std::endl;
   double ActualThresh = norm_frob(inner_prod(Initial, f));
   while (int(Alpha.size()) <= MaxIter && ActualThresh < Threshold && Beta.back() > Eps)
   {
      Hf = operator_prod(conj(H.Center()), H_left, f, herm(H_right));
      Alpha.push_back(inner_prod(f, Hf));
      Hf -= Alpha.back() * f + Beta.back() * fLast;
      DEBUG_TRACE(inner_prod(Hf, fLast))(inner_prod(Hf, f));
      fLast = f;
      f = Hf;
      Beta.push_back(norm_frob(f));
      f *= 1.0 / Beta.back();
      LeftVector.push_back(inner_prod(LeftL, f));
      std::cout << std::setw(3) << (Alpha.size()-1) << "   "
                << std::setw(20) << Alpha.back().real() << "   "
                << std::setw(20) << Beta.back() << "  "
             << std::setw(20) << LeftVector.back().real() << "  "
             << std::setw(20) << LeftVector.back().imag() << std::endl;
      ActualThresh = norm_frob(inner_prod(Initial, f));
   }
   if (ShowStoppingReason)
   {
      if (int(Alpha.size()) >= MaxIter)
      {
         std::cerr << "Stopped on max-iter.\n";
         std::cerr << "Threshold = " << ActualThresh << '\n';
      }
      else if (ActualThresh >= Threshold)
      {
         std::cerr << "Stopped on orthogonality threshold of " << ActualThresh << '\n';
      }
      else if (Beta.back() <= Eps)
      {
         std::cerr << "Stopped on Beta below epsilon.\n";
      }
      else PANIC("WTF?");
   }
}

int main(int argc, char** argv)
{
   try
   {
      int Location = -1;
      int MaxIter = 100000;
      double Threshold = 1e-2;
      double GroundstateEnergy = 0;
      int Verbosity = 0;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(),
          "operator to use for the Hamiltonian (right Lanczos vector attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::vector<std::string> >(),
          "input wavefunction to generate the effective basis (zero or more)")
         ("right,r", prog_opt::value<std::string>(),
          "right-hand Lanczos vector |R> (defaults to the first --wavefunction)")
         ("left,l", prog_opt::value<std::string>(),
          "Left-hand Lanczos vector <L| (defaults to --right)")
         ("GroundstateEnergy,G", prog_opt::value(&GroundstateEnergy),
          "groundstate energy of the Hamiltonian (wavefunction attribute"
          " \"GroundstateEnergy\", not needed if --no-preamble)")
         ("max-states,m", prog_opt::value<int>(&MaxStates),
          ("Maximum number of states to keep in the effective basis [default "
           + boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("min-trunc,t", prog_opt::value<double>(&MinTrunc),
          ("Minimum desired truncation error per site of the effective basis [default "
           + boost::lexical_cast<std::string>(MinTrunc) + "]").c_str())
         ("bond,b", prog_opt::value<int>(&Location),
          "Generate the basis at this bond, valid is 1 .. L-1 [default L/2]")
         ("max-iter,i", prog_opt::value<int>(&MaxIter),
          "maximum number of iterations")
         ("threshold,s", prog_opt::value<double>(&Threshold),
          ("stopping criteria for the orthogonality of the Krylov vector <kn|k0> [default "
           + boost::lexical_cast<std::string>(Threshold) + "]").c_str())
         ("epsilon,e", prog_opt::value<double>(&Eps),
          ("stopping criteria for the magnitude of Beta [default "
           + boost::lexical_cast<std::string>(Eps) + "]").c_str())
         ("no-preamble", prog_opt::bool_switch(),
          "only show the tridiagonal coefficients, don't show the groundstate energy")
         ("verbose,v", prog_opt_ext::accum_value(&Verbosity), "increase verbosity (can be used more than once)")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || (vm.count("wavefunction") == 0 && vm.count("right") == 0))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-tridiag [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      ShowStates = Verbosity >= 2;
      ShowTitles = Verbosity >= 1;
      ShowStoppingReason = Verbosity >= 1;

      std::vector<std::string> InputWavefunctions;
      if (vm.count("wavefunction") != 0)
         InputWavefunctions = vm["wavefunction"].as<std::vector<std::string> >();

      // Load the wavefunctions
      std::string TempFile = getenv_or_default("MP_BINPATH", std::string());
      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      if (Verbosity >= 2)
      {
         std::cerr << "Using page size " << PageSize
                   << ", cache size " << CacheSize << '\n';
      }
      int TempFileDesc = ProcControl::CreateUniqueFile(TempFile);
      CHECK(TempFileDesc != -1);  // if this happens, then we failed to create a temporary file
      pheap::Initialize(TempFile, 1, PageSize, CacheSize);

      std::vector<CenterWavefunction> Psi;
      for (unsigned i = 0; i < InputWavefunctions.size(); ++i)
      {
         if (Verbosity >= 1)
            std::cerr << "Loading wavefunction: " << InputWavefunctions[i] << '\n';
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(InputWavefunctions[i]);
         Psi.push_back(*P);
      }

      // Load the Lanczos vectors
      int RightLanczos, LeftLanczos;
      if (vm.count("right") != 0)
      {
         std::string WavefuncName = vm["right"].as<std::string>();
         if (Verbosity >= 1)
            std::cerr << "Loading wavefunction: " << WavefuncName << '\n';
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(WavefuncName);
         RightLanczos = Psi.size();
         Psi.push_back(*P);
      }
      else
         RightLanczos = 0; // the default, if --right was not supplied

      if (vm.count("left") != 0)
      {
         std::string WavefuncName = vm["left"].as<std::string>();
         if (Verbosity >= 1)
            std::cerr << "Loading wavefunction: " << WavefuncName << '\n';
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(WavefuncName);
         LeftLanczos = Psi.size();
         Psi.push_back(*P);
      }
      else
         LeftLanczos = RightLanczos; // the default, if --left was not supplied

      // Set up the Hamiltonian
      std::string HamString;
      if (vm.count("Hamiltonian") == 1)
      {
         HamString = vm["Hamiltonian"].as<std::string>();
      }
      else
      {
         if (Psi[RightLanczos].Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: No Hamiltonian specified.\n"
               "Specify a Hamiltonian either with --Hamiltonian parameter, "
               "or as an attribute of the right Lanczos vector.\n";
            return 1;
         }
         HamString = Psi[RightLanczos].Attributes()["Hamiltonian"].as<std::string>();
      }
      if (Verbosity >= 1)
         std::cerr << "Hamiltonian: " << HamString << std::endl;
      SplitOperator Hamiltonian = ParseOperator(HamString);

      // Set up the location bond
      if (Location == -1)
      {
         Location = Psi[0].size() / 2;
      }

      // Adjust the truncation error for the number of bonds; no longer necessary as this is now per site
      // MinTrunc /= (Psi[0].size()-1);

      // Set up the std::cout flags
      std::cout.precision(16);
      std::cout.setf(std::ios::showpoint);
      std::cout.setf(std::ios::fixed);

      // Show the norm and groundstate energy as the first lines of output, unless --no-preamble was specified
      if (!vm["no-preamble"].as<bool>())
      {
         if (vm.count("GroundstateEnergy") == 0)
         {
            if (Psi[0].Attributes().count("GroundstateEnergy") == 0)
            {
               std::cerr << "fatal: no groundstate energy specified.\n";
               return 1;
            }
            else
               GroundstateEnergy = Psi[0].Attributes()["GroundstateEnergy"].as<double>();
         }

         if (ShowTitles)
            std::cerr << "Groundstate energy\n";
         std::cout << GroundstateEnergy << '\n';
      }
      Aggregator Ag(Psi, RightLanczos, LeftLanczos, Hamiltonian, Location);
      Ag.Tridiagonalize(MaxIter, Threshold);

      pheap::Shutdown();

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
