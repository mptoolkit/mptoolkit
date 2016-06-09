// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-local-minres.cpp
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
#include "matrixproduct/centerwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/dmrg.h"
#include "mp-algorithms/ddmrg_functional.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include <boost/program_options.hpp>
#include "common/prog_opt_accum.h"
#include "interface/inittemp.h"
#include <iostream>

namespace prog_opt = boost::program_options;

bool ShowStates = false; // for verbose output
bool ShowTitles = false;
bool ShowStoppingReason = false;

int MinStates = 1;
int MaxStates = 100000;
double TruncCutoff = 0;
double EigenCutoff = -1;
StatesInfo SInfo;
int Verbosity = 0;

typedef std::complex<double> complex;

struct H2Multiply
{
   H2Multiply(MPStateComponent const& H_left_, MPStateComponent const& H_right_,
              MPStateComponent const& H2_left_, MPStateComponent const& H2_right_,
              double FreqPlusEnergy_, double Broad_)
      : H_left(H_left_), H_right(H_right_), H2_left(H2_left_), H2_right(H2_right_),
        FreqPlusEnergy(FreqPlusEnergy_), Broad(Broad_) {}

   MatrixOperator operator()(MatrixOperator const& c) const
   {
      MatrixOperator Hc = operator_prod(H_left, c, herm(H_right));
      MatrixOperator H2c = operator_prod(H2_left, c, herm(H2_right));
      return H2c - 2*FreqPlusEnergy*Hc + (FreqPlusEnergy*FreqPlusEnergy + Broad*Broad)*c;
   }

   MPStateComponent const& H_left;
   MPStateComponent const& H_right;
   MPStateComponent const& H2_left;
   MPStateComponent const& H2_right;
   double FreqPlusEnergy, Broad;
};

struct HMultiply
{
   HMultiply(MPStateComponent const& H_left_, MPStateComponent const& H_right_,
              double FreqPlusEnergy_, double Broad_)
      : H_left(H_left_), H_right(H_right_),
        FreqPlusEnergy(FreqPlusEnergy_), Broad(Broad_), Neg(false) {}

   MatrixOperator operator()(MatrixOperator const& c) const
   {
      MatrixOperator Hc = operator_prod(H_left, c, herm(H_right));
      Neg = !Neg;
      return std::complex<double>(FreqPlusEnergy, Neg ? -Broad : Broad)*c - Hc;
   }

   MPStateComponent const& H_left;
   MPStateComponent const& H_right;
   double FreqPlusEnergy, Broad;
   mutable bool Neg;
};

class Aggregator
{
   public:
      typedef MPStateComponent Component;
      typedef Component::OperatorType OperatorType;

      Aggregator(std::vector<CenterWavefunction> const& Psi_, int NumCV_, 
                 SplitOperator const& H_, SplitOperator const& H2_,  int Location_);

      MatrixOperator RightLV() const { return Center[RightLanczos]; }

      MatrixOperator Wavefunction(int i) const { return Center[i]; }
   
      CenterWavefunction const& FullWavefunction() const { return Result; }

      complex ConstructCV(MatrixOperator& Guess, double FreqPlusEnergy, double Broadening, 
                          int MaxSubspaceSize, int MaxIter,
                          double& Resid, double StopDelta, bool Precondition = true);

   private:
      void RotateRight();
      void ConstructLeft();

      void RotateLeft();
      void ConstructRight();

      std::vector<CenterWavefunction> Psi;
      int NumCV; // the number of elements of Psi that represent CV's.  the remainder are Lanczos vectors.
      int RightLanczos;
      SplitOperator H;
      SplitOperator H2;

      std::vector<OperatorType> LeftMap, RightMap;
      std::vector<double> Weights; // density matrix weights for each state

      CenterWavefunction Result;

      std::vector<OperatorType> Center;
      MPStateComponent H_left, H_right;
      MPStateComponent H2_left, H2_right;

      MPStateComponent HL_left, HL_right; // matrix elements of H in mixed basis <cv|H|lv>

      QuantumNumber Ident;
};

void Aggregator::ConstructLeft()
{
   // Construct the mapping from the vectors to the result space
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      LeftMap[i] = operator_prod(herm(Result.Left()), LeftMap[i], Psi[i].Left());
   }
                      
   // Construct the density matrix, taking only the correction vectors, not the Lanczos vectors
   OperatorType Rho;
   for (int i = 0; i < NumCV; ++i)
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
                                                                                           SInfo,
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
   H2_left = operator_prod(herm(H2.Left()), herm(Result.Left()), H2_left, Result.Left());

   // Matrix elements connecting the LV and CV
   HL_left = operator_prod(herm(H.Left()), herm(Result.Left()), HL_left, Psi[RightLanczos].Left());
}

void Aggregator::ConstructRight()
{
   // Construct the mapping from the vectors to the result space
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      RightMap[i] = operator_prod(Result.Right(), RightMap[i], herm(Psi[i].Right()));
   }

   // Construct the density matrix, taking only the correction vectors, not the Lanczos vectors
   OperatorType Rho;
   for (int i = 0; i < NumCV; ++i)
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
                                                                                           SInfo,
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
   H2_right = operator_prod(H2.Right(), Result.Right(), H2_right, herm(Result.Right()));

   // Matrix elements connecting the LV and CV
   HL_right = operator_prod(H.Right(), Result.Right(), HL_right, herm(Psi[RightLanczos].Right()));
}

void Aggregator::RotateLeft()
{
   H.RotateLeft();
   H2.RotateLeft();
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
   H2.RotateRight();
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      Psi[i].RotateRight();
   }
   Result.PushLeft(Component::ConstructFullBasis2(Result.Left().Basis2(),
                                                  Psi[0].Left().SiteBasis()));
}

Aggregator::Aggregator(std::vector<CenterWavefunction> const& Psi_,  int NumCV_, 
                       SplitOperator const& H_, SplitOperator const& H2_, 
                       int Location_)
   : Psi(Psi_), NumCV(NumCV_), RightLanczos(NumCV_),
     H(H_), H2(H2_), Weights(Psi_.size(), 1.0), Ident(Psi_[0].GetSymmetryList())
{
   // Rotate the wavefunctions to left most position
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      while (Psi[i].LeftSize() > 1)
         Psi[i].RotateLeft();
   }
   while (H.LeftSize() > 1)
      H.RotateLeft();
   while (H2.LeftSize() > 1)
      H2.RotateLeft();

   // Initialize the result matrices
   OperatorType LeftVac = OperatorType::make_identity(Psi[0].LeftVacuumBasis());
   LeftMap = std::vector<OperatorType>(Psi.size(), LeftVac);

   H_left = make_vacuum_state(Psi[0].LookupLeft(0).Basis1()[0]);
   H2_left = make_vacuum_state(Psi[0].LookupLeft(0).Basis1()[0]);
   HL_left = make_vacuum_state(Psi[0].LookupLeft(0).Basis1()[0]);
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
   }
   while (H.RightSize() > 1)
      H.RotateRight();
   while (H2.RightSize() > 1)
      H2.RotateRight();

   OperatorType RightVac = OperatorType::make_identity(Psi[0].RightVacuumBasis());
   RightMap = std::vector<OperatorType>(Psi.size(), RightVac);

   H_right = make_vacuum_state(Psi[0].GetSymmetryList());
   H2_right = make_vacuum_state(Psi[0].GetSymmetryList());
   HL_right = make_vacuum_state(Psi[0].GetSymmetryList());
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

complex Aggregator::ConstructCV(MatrixOperator& Guess, double FreqPlusEnergy, double Broadening, 
                                int MaxSubspaceSize,
                                int MaxIter,
                                double& Resid, double StopDelta, bool Precondition)
{
   // the lanczos vector
   MatrixOperator Lanc = Psi[RightLanczos].Center();

   double LancNorm = norm_frob_sq(Lanc); // squared norm of the Lanczos vector, to calculate the residual

   // A*|lv>, projected onto our local Hilbert space, negated
   // Center[RightLanczos] is the Lanczos vector projected onto our basis
   MatrixOperator b = complex(-FreqPlusEnergy, Broadening)*Center[RightLanczos]
      + operator_prod(HL_left, Lanc, herm(HL_right));

   // Multiply functor
   H2Multiply Mult(H_left, H_right, H2_left, H2_right, FreqPlusEnergy, Broadening);
   HMultiply PreMult(H_left, H_right, FreqPlusEnergy, Broadening);

   int Iter = 1;
   // do the actual minimization
   if (Verbosity >= 2)
      std::cerr << "Starting minimization.\n";
   int m = MaxSubspaceSize;
   double OldResid = 0;
   if (Precondition)
      Resid = FunctionalMinimizeWithPre(Guess, Mult, b, m, PreMult);
   else
      Resid = FunctionalMinimize(Guess, Mult, b, m, LinearAlgebra::Identity<MatrixOperator>());
   Resid += LancNorm;

   while (fabs(Resid-OldResid) > StopDelta && Iter < MaxIter)
   {
      if (Verbosity >= 2)
         std::cerr << "Restarting minimization, resid=" << Resid << "\n";
      m = MaxSubspaceSize;
      OldResid = Resid;
      if (Precondition)
         Resid = FunctionalMinimizeWithPre(Guess, Mult, b, m, PreMult);
      else
         Resid = FunctionalMinimize(Guess, Mult, b, m, LinearAlgebra::Identity<MatrixOperator>());
      Resid += LancNorm;
      ++Iter;
   }
   if (Verbosity >= 2)
   {
      std::cerr << "Stopping minimizaion after " << Iter << " iterations, ";
      if (Iter >= MaxIter)
         std::cerr << "iteration count has hit limit.\n";
      else if (fabs(Resid-OldResid) > StopDelta)
         std::cerr << "change in residual at last iteration " << fabs(Resid-OldResid) 
                   << " is smaller than threshold.\n";
      else
         std::cerr << "the moon is in the wrong cycle.\n";
   }
   complex x = inner_prod(Center[RightLanczos], Guess);
   return x;
}

int main(int argc, char** argv)
{
   try
   {
      int Location = -1;
      double GroundstateEnergy = 0;
      std::vector<double> Frequencies;
      std::vector<std::string> FrequenciesAsString; // string version of the freq for the output file
      double Broadening;
      std::string OutputPrefix;
      int KrylovSize = 100;
      std::vector<std::string> LanczosFiles;
      std::vector<std::string> InputWavefunctions;
      bool NoPrecondition = false;
      

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(), 
          "operator to use for the Hamiltonian (rLanczos vector attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value(&InputWavefunctions),
          "input wavefunction to generate the effective basis (one or more)")
         ("lanczos,l", prog_opt::value(&LanczosFiles),
          "input Lanczos vector (required, does not form part of the effective basis)")
         ("Broadening,B", prog_opt::value(&Broadening),
          "Broadening")
         ("GroundstateEnergy,G", prog_opt::value(&GroundstateEnergy),
          "groundstate energy of the Hamiltonian (Lanczos vector attribute"
          " \"GroundstateEnergy\"")
         ("min-states", prog_opt::value<int>(&MinStates), 
          ("Minimum number of states to keep [default " 
           +boost::lexical_cast<std::string>(MinStates)+"]").c_str())
	 ("max-states,m", prog_opt::value<int>(&MaxStates), 
          ("Maximum number of states to keep [default "
           +boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("trunc,r", prog_opt::value(&TruncCutoff), 
          ("Cutoff truncation error per site [default "
           +boost::lexical_cast<std::string>(TruncCutoff)+"]").c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff), 
          ("Cutoff threshold for density matrix eigenvalues [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
         ("out,o", prog_opt::value(&OutputPrefix),
          "save the correction vectors using this filename prefix followed by .Frequency")
         ("Frequency,F", prog_opt::value(&FrequenciesAsString),
          "list of frequencies to evaluate the correction vectors, only used if --out is specified")
         ("bond,b", prog_opt::value<int>(&Location),
          "Generate the basis at this bond, valid is 1 .. L-1 [default L/2]")
         ("subspace-size,s", prog_opt::value<int>(&KrylovSize),
          ("maximum size of the Krylov subspace for the functional minimization [default"
           + boost::lexical_cast<std::string>(KrylovSize) + "]").c_str())
         ("precision,p", prog_opt::value(&CVPrecision),
          ("stop building the Krylov subspace when the input wavefunctions can be reproduced with this "
           "fidelity loss ["+boost::lexical_cast<std::string>(CVPrecision)+"]").c_str())
         ("verbose,v", prog_opt_ext::accum_value(&Verbosity),
          "increase verbosity (can be used more than once)")
	  ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("lanczos") == 0) 
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-tridiag [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // Convert the frequencies to double
      for (unsigned i = 0; i < FrequenciesAsString.size(); ++i)
         Frequencies.push_back(boost::lexical_cast<double>(FrequenciesAsString[i]));

      bool const Precondition = !NoPrecondition;

      ShowStates = Verbosity >= 1;
      ShowTitles = Verbosity >= 1;
      ShowStoppingReason = Verbosity >= 2;

      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      // Load the wavefunctions
      mp_pheap::InitializeTempPHeap(Verbosity >= 2);

      std::vector<CenterWavefunction> Psi;
      for (unsigned i = 0; i < InputWavefunctions.size(); ++i)
      {
         if (Verbosity >= 2)
            std::cerr << "Loading wavefunction: " << InputWavefunctions[i] << '\n';
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(InputWavefunctions[i]);
         Psi.push_back(*P);
      }
      int const NumCV = Psi.size();

      for (unsigned i = 0; i < LanczosFiles.size(); ++i)
      {
         if (Verbosity >= 2)
            std::cerr << "Loading Lanczos wavefunction: " << LanczosFiles[i] << '\n';
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(LanczosFiles[i]);
         Psi.push_back(*P);
      }

      // Set up the Hamiltonian
      std::string HamString;
      if (vm.count("Hamiltonian") == 1)
      {
         HamString = vm["Hamiltonian"].as<std::string>();
      }
      else
      {
         if (Psi[NumCV].Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: No Hamiltonian specified.\n"
               "Specify a Hamiltonian either with --Hamiltonian parameter, "
               "or as an attribute of the right Lanczos vector.\n";
            return 1;
         }
         HamString = Psi[NumCV].Attributes()["Hamiltonian"].as<std::string>();
      }
      if (Verbosity >= 2)
         std::cerr << "Hamiltonian: " << HamString << std::endl;
      MPOperator Ham = ParseOperator(HamString);
      MPOperator Ham2 = Ham*Ham;
      SplitOperator Hamiltonian = Ham;
      SplitOperator Hamiltonian2 = Ham2;

      // Set up the location bond
      if (Location == -1)
      {
         Location = Psi[0].size() / 2;
      }

      // Adjust the truncation error for the number of bonds; no longer necessary as this is now per site
      // MinTrunc /= (Psi[0].size()-1);

      // Set up the std::cout flags - do we want this?
      //      std::cout.setf(std::ios::showpoint);
      //      std::cout.setf(std::ios::fixed);

      if (vm.count("GroundstateEnergy") == 0)
      {
         if (Psi[NumCV].Attributes().count("GroundstateEnergy") == 0)
         {
            std::cerr << "fatal: no groundstate energy specified.\n";
            return 1;
         }
         else
            GroundstateEnergy = Psi[NumCV].Attributes()["GroundstateEnergy"].as<double>();
      }

      Aggregator Ag(Psi, NumCV, Hamiltonian, Hamiltonian2, Location);
      MatrixOperator Guess = Ag.Wavefunction(0);  // our guess of the first correction vector
      std::cout << "         #Frequency                #G_Real               #G_Imag                #Resid\n";
      for (unsigned i = 0; i < Frequencies.size(); ++i)
      {
         double Resid = 0;
         complex x = Ag.ConstructCV(Guess, Frequencies[i]+GroundstateEnergy, Broadening, KrylovSize, 
                                    MaxIter, Resid, StopDelta, Precondition);
         std::cout << std::setw(20) << Frequencies[i] << "  "  
                   << std::setw(20) << x.real() << "  "
                   << std::setw(20) << x.imag() << "  "
                   << std::setw(20) << Resid << std::endl;
         if (!OutputPrefix.empty())
         {
            std::string FName = OutputPrefix + ".F" + FrequenciesAsString[i];
            if (Verbosity >= 2)
               std::cerr << "Saving wavefunction: " << FName << '\n';

            CenterWavefunction Full = Ag.FullWavefunction();
            Full.Center() = Guess;
            pvalue_ptr<MPWavefunction> Psi = new MPWavefunction(Full.AsLinearWavefunction());
            pheap::ExportHeap(FName, Psi);
         }
      }

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
