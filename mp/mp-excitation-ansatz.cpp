// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-excitation-ansatz.cpp
//
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "common/environment.h"
#include "common/formatting.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/rangelist.h"
#include "common/terminal.h"
#include "common/unique.h"
#include "interface/inittemp.h"
#include "lattice/infinite-parser.h"
#include "lattice/infinitelattice.h"
#include "linearalgebra/arpack_wrapper.h"
#include "mp-algorithms/excitation-ansatz.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      int Quiet = 0;
      std::string InputFileLeft;
      std::string InputFileRight;
      std::string HamStr;
      std::string KStr = "0";
      double Tol = 1e-10;
      int NumEigen = 1;
      std::string QuantumNumber;
      int Rotate = 0;
      std::string String;
      std::string OutputPrefix;
      int OutputDigits = -1;
      bool Random = false;
      bool Streaming = false;
      bool NoStreaming = false;

      EASettings Settings;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "Operator to use for the Hamiltonian (if unspecified, use wavefunction attribute \"Hamiltonian\")")
         ("momentum,k", prog_opt::value(&KStr), FormatDefault("The momentum, in units of pi: this can be a single number, or a range of the form start:end:step or start:end,num", KStr).c_str())
         ("ky", prog_opt::value(&Settings.ky), "Target this value of the y-momentum [2D cylinders]")
         ("alpha", prog_opt::value(&Settings.Alpha), FormatDefault("Energy parameter to penalize states with the wrong y-momentum [2D Cylinders]", Settings.Alpha).c_str())
         ("numeigen,n", prog_opt::value<int>(&NumEigen),
          FormatDefault("The number of lowest eigenvalues to calculate", NumEigen).c_str())
         ("quantumnumber,q", prog_opt::value(&QuantumNumber),
          "The quantum number sector for the excitation [default identity]")
         ("rotate,r", prog_opt::value(&Rotate), "Rotate the right boundary by this many sites to the left")
         ("string", prog_opt::value(&String),
          "Use this string MPO representation for the cylinder translation operator")
         ("output,o", prog_opt::value(&OutputPrefix), "Prefix for saving output files (will not save if not specified)")
         ("digits", prog_opt::value(&OutputDigits), "Manually use this number of decimal places for the filenames")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Error tolerance for the eigensolver", Tol).c_str())
         ("gmrestol", prog_opt::value(&Settings.GMRESTol),
          FormatDefault("Error tolerance for the GMRES algorithm", Settings.GMRESTol).c_str())
         ("seed", prog_opt::value<unsigned long>(), "Random number generator seed")
         ("random", prog_opt::bool_switch(&Random), "Use a random initial state for each momentum (otherwise, use the previous result as an initial guess)")
         ("streaming", prog_opt::bool_switch(&Streaming), "Store the left and right strips by reference to the input files")
         ("no-streaming", prog_opt::bool_switch(&NoStreaming), "Store the left and right strips into the output file [default]")
         ("quiet", prog_opt_ext::accum_value(&Quiet), "Hide column headings, use twice to hide momentum")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&InputFileLeft), "psi")
         ("psi2", prog_opt::value(&InputFileRight), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> [psi-right]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      if ((vm.count("ky") || vm.count("alpha")) && vm.count("string") == 0)
      {
         std::cerr << "fatal: Please specify the cylinder translation operator using --string." << std::endl;
         return 1;
      }

      if (Streaming && NoStreaming)
      {
         std::cerr << "fatal: Cannot use --streaming and --no-streaming simultaneously!" << std::endl;
         return 1;
      }
      else if (!Streaming && !NoStreaming)
         NoStreaming = true; // This is the current default behavior.

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      Settings.Verbose = Verbose;

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX)
         : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      mp_pheap::InitializeTempPHeap();

      RangeList KList(KStr);

      pvalue_ptr<MPWavefunction> InPsiLeft = pheap::ImportHeap(InputFileLeft);
      InfiniteWavefunctionLeft PsiLeft = InPsiLeft->get<InfiniteWavefunctionLeft>();

      if (HamStr.empty())
      {
         if (InPsiLeft->Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: no Hamiltonian specified, use -H option or set wavefunction attribute Hamiltonian." << std::endl;
            return 1;
         }
         HamStr = InPsiLeft->Attributes()["Hamiltonian"].as<std::string>();
      }

      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO;
      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

      if (vm.count("string"))
      {
         std::tie(Settings.StringOp, std::ignore) = ParseProductOperatorAndLattice(String);
      }

      // The number of lattice unit cells in Psi.
      int LatticeUCsPerPsiUC = PsiLeft.size() / Lattice.GetUnitCell().size();
      // Scale the momentum by the number of lattice unit cells in the unit cell of PsiLeft.
      Settings.k = KList.get_start() * LatticeUCsPerPsiUC;

      if (OutputDigits == -1)
      {
         OutputDigits = std::max(formatting::digits(KList.get_start()), formatting::digits(KList.get_step()));
         if (vm.count("ky"))
            OutputDigits = std::max(OutputDigits, formatting::digits(Settings.ky));
      }

      // Initialize the effective Hamiltonian.
      HEff H;
      InfiniteWavefunctionRight PsiRight;
      if (vm.count("psi2"))
      {
         pvalue_ptr<MPWavefunction> InPsiRight = pheap::ImportHeap(InputFileRight);
         if (InPsiRight->is<InfiniteWavefunctionRight>())
            PsiRight = InPsiRight->get<InfiniteWavefunctionRight>();
         else if (InPsiRight->is<InfiniteWavefunctionLeft>())
         {
            if (Streaming)
            {
               std::cerr << "fatal: psi-right must be an InfiniteWavefunctionRight if streaming is enabled." << std::endl;
               return 1;
            }

            PsiRight = InfiniteWavefunctionRight(InPsiRight->get<InfiniteWavefunctionLeft>());
         }
         else
         {
            std::cerr << "fatal: psi-right must be an InfiniteWavefunctionLeft or InfiniteWavefunctionRight." << std::endl;
            return 1;
         }
      }
      else
      {
         // This condition could be relaxed, in that case we would only save
         // the left boundary by reference, but we assume that the user want
         // both boundaries saved by reference.
         if (Streaming)
         {
            std::cerr << "fatal: psi-right must be specified if streaming is enabled." << std::endl;
            return 1;
         }

         PsiRight = InfiniteWavefunctionRight(PsiLeft);
      }

      // Save the original read-only input wavefunctions for the output files.
      InfiniteWavefunctionLeft PsiLeftOriginal = PsiLeft;
      InfiniteWavefunctionRight PsiRightOriginal = PsiRight;

      // Shift the left boundary if we want a different quantum number for the excitation.
      QuantumNumbers::QuantumNumber I(HamMPO[0].GetSymmetryList()); // The identity quantum number.
      QuantumNumbers::QuantumNumber Q = I;
      if (vm.count("quantumnumber"))
      {
         Q = QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), QuantumNumber);
         inplace_qshift(PsiLeft, Q);
      }

      // Rotate PsiRight: making sure that Rotate is in [0, PsiRight.size()-1].
      Rotate = numerics::divp(Rotate, PsiRight.size()).rem;
      PsiRight.rotate_left(Rotate);

      // Turn off momentum targeting if ky is unspecified.
      if (vm.count("ky") == 0)
         Settings.Alpha = 0.0;

      H = HEff(PsiLeft, PsiRight, HamMPO, Settings);

      // Print column headers.
      if (Quiet == 0)
      {
         if (vm.count("string"))
            std::cout << "#kx/pi              ";
         else
            std::cout << "#k/pi               ";
         if (NumEigen > 1)
            std::cout << "#n        ";
         if (vm.count("string"))
            std::cout << "#ky/pi              #|Ty|               ";
         std::cout << "#E" << std::endl;
      }
      std::cout << std::left;

      std::vector<std::complex<double>> Guess;

      // Calculate the excitation spectrum for each k desired.
      for (double const k : KList)
      {
         H.SetK(k * LatticeUCsPerPsiUC);

         PackHEff PackH = PackHEff(H);
         std::vector<std::complex<double>> EVectors;

         LinearAlgebra::Vector<std::complex<double>> EValues
            = LinearAlgebra::DiagonalizeARPACK(PackH, PackH.size(), NumEigen, LinearAlgebra::WhichEigenvalues::SmallestReal,
                                               Guess.data(), Tol, &EVectors, 0, true, Verbose);

         if (!Random)
         {
            // Save the smallest eigenvector as the initial guess for the next iteration.
            Guess = std::vector<std::complex<double>>(PackH.size());
            std::copy(&(EVectors[0]), &(EVectors[PackH.size()]), Guess.data());
         }

         // Print results for this k.
         auto E = EValues.begin();
         for (int i = 0; i < NumEigen && E != EValues.end(); ++i, ++E)
         {
            int Index = i * PackH.size();
            if (Quiet < 2)
               std::cout << std::setw(20) << k;
            if (NumEigen > 1)
               std::cout << std::setw(10) << i;
            if (vm.count("string") && Quiet < 2)
            {
               std::deque<MatrixOperator> XDeque = PackH.unpack(&(EVectors[Index]));
               std::cout << std::setw(20) << std::arg(H.Ty(XDeque))/math_const::pi
                         << std::setw(20) << std::abs(H.Ty(XDeque));
            }
            std::cout << std::setw(20) << std::real(*E) + Settings.Alpha << std::endl;

            // Save wavefunction.
            if (OutputPrefix != "")
            {
               std::vector<WavefunctionSectionLeft> WindowVec = H.ConstructWindowVec(PackH.unpack(&(EVectors[Index])));
               EAWavefunction PsiEA(PsiLeftOriginal, WindowVec, PsiRightOriginal, Q, I, 0, Rotate, H.ExpIK);

               if (Streaming)
               {
                  PsiEA.set_left_filename(InputFileLeft);
                  PsiEA.set_right_filename(InputFileRight);
               }

               MPWavefunction Wavefunction;
               Wavefunction.Wavefunction() = std::move(PsiEA);
               Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
               Wavefunction.SetDefaultAttributes();
               Wavefunction.Attributes()["Prefix"] = OutputPrefix;
               Wavefunction.Attributes()["ExcitationEnergy"] = std::real(*E) + Settings.Alpha;
               Wavefunction.Attributes()["Hamiltonian"] = HamStr;

               std::string FName = OutputPrefix;
               if (vm.count("ky") == 0)
                  FName += ".k" + formatting::format_digits(k, OutputDigits);
               else
                  FName += ".kx" + formatting::format_digits(k, OutputDigits)
                         + ".ky" + formatting::format_digits(Settings.ky, OutputDigits);
               if (NumEigen > 1)
                  FName += ".n" + std::to_string(i);

               pvalue_ptr<MPWavefunction> PsiPtr(new MPWavefunction(Wavefunction));
               pheap::ExportHeap(FName, PsiPtr);
            }
         }
      }
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      pheap::Cleanup();
      return 1;
   }
}
