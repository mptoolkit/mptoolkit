// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-create.cpp
//
// Copyright (C) 2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "lattice/infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "linearalgebra/arpack_wrapper.h"
#include "mps/packunpack.h"

namespace prog_opt = boost::program_options;

struct CMultiply
{
   CMultiply(StateComponent const& E_, StateComponent const& F_)
      : P(E_.Basis2(), F_.Basis2(), E_[0].TransformsAs()), E(E_), F(F_)
   {
   }

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator C = P.unpack(In);
      C = operator_prod(E, C, herm(F));
      P.pack(C, Out);
   }

   PackMatrixOperator P;
   StateComponent const& E;
   StateComponent const& F;
};

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OutputFilename;
      std::string LeftFilename;
      std::string RightFilename;
      std::string HamStr;
      std::string QuantumNumber;
      bool Random = false;
      bool BoundaryLambda = false;
      bool Force = false;
      double GMRESTol = 1e-13;
      double Tol = 1e-15;
      bool Streaming = false;
      bool NoStreaming = false;
      int NumEigen = 1;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("output,o", prog_opt::value(&OutputFilename), "Output IBC filename [required]")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "Operator to use for the Hamiltonian (if unspecified, use wavefunction attribute \"Hamiltonian\")")
         ("quantumnumber,q", prog_opt::value(&QuantumNumber),
          "Shift the right boundary by this quantum number [default identity]")
         ("tol", prog_opt::value(&Tol), FormatDefault("Tolerance of the eigensolver", Tol).c_str())
         ("gmrestol", prog_opt::value(&GMRESTol), FormatDefault("Error tolerance for the GMRES algorithm", GMRESTol).c_str())
         ("streaming", prog_opt::bool_switch(&Streaming), "Store the left and right strips by reference to the input files")
         ("no-streaming", prog_opt::bool_switch(&NoStreaming), "Store the left and right strips into the output file [default]")
         ("random", prog_opt::bool_switch(&Random), "Use a random lambda matrix instead of minimizing the Hamiltonian")
         ("boundary-lambda", prog_opt::bool_switch(&BoundaryLambda), "Use the left boundary's lambda matrix instead of minimizing the Hamiltonian")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi-left", prog_opt::value(&LeftFilename), "psi-left")
         ("psi-right", prog_opt::value(&RightFilename), "psi-right")
         ;

      prog_opt::positional_options_description p;
      p.add("psi-left", 1);
      p.add("psi-right", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi-left") == 0 || vm.count("output") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> [psi-right] -o <psi-out>" << std::endl;
         std::cerr << desc << std::endl;
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

      pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      pvalue_ptr<MPWavefunction> InPsiLeft = pheap::ImportHeap(LeftFilename);
      InfiniteWavefunctionLeft PsiLeft = InPsiLeft->get<InfiniteWavefunctionLeft>();

      InfiniteWavefunctionRight PsiRight;
      if (vm.count("psi-right"))
      {
         pvalue_ptr<MPWavefunction> InPsiRight = pheap::ImportHeap(RightFilename);
         if (InPsiRight->is<InfiniteWavefunctionRight>())
            PsiRight = InPsiRight->get<InfiniteWavefunctionRight>();
         else if (InPsiRight->is<InfiniteWavefunctionLeft>())
         {
            if (!InPsiRight->is<InfiniteWavefunctionRight>() && Streaming)
            {
               std::cerr << "fatal: right_psi must be an InfiniteWavefunctionRight if streaming is enabled." << std::endl;
               return 1;
            }

            PsiRight = InfiniteWavefunctionRight(InPsiRight->get<InfiniteWavefunctionLeft>());
         }
         else
         {
            std::cerr << "fatal: right_psi must be an InfiniteWavefunctionLeft or InfiniteWavefunctionRight." << std::endl;
            return 1;
         }
      }
      else
      {
         // This condition could be relaxed, in that case we would only save
         // the left boundary by reference, but we assume that the user wants
         // both boundaries saved by reference.
         if (Streaming)
         {
            std::cerr << "fatal: psi-right must be specified if streaming is enabled." << std::endl;
            return 1;
         }

         PsiRight = InfiniteWavefunctionRight(PsiLeft);
         // In this case, we shouldn't need to solve for the ground state.
         BoundaryLambda = true;
      }

      // Get the QShift for the right boundary.
      QuantumNumbers::QuantumNumber QRight(PsiRight.GetSymmetryList());

      // If we specified the quantum number as an option.
      if (vm.count("quantumnumber"))
         QRight = QuantumNumbers::QuantumNumber(PsiRight.GetSymmetryList(), QuantumNumber);
      // If the bases of the two boundary unit cells have only one quantum
      // number sector, manually ensure that they match.
      else if (PsiLeft.Basis1().size() == 1 && PsiRight.Basis1().size() == 1)
      {
         QRight = delta_shift(PsiLeft.Basis1()[0], adjoint(PsiRight.Basis1()[0]));
         // If QRight is non-Abelian, we cannot shift by it.
         if (QRight.degree() != 1)
         {
            std::cerr << "fatal: The left and right boundaries do not have any possible common Abelian quantum number sectors." << std::endl;
            return 1;
         }
      }

      if (Verbose > 0)
         std::cout << "Shifting PsiRight by " << QRight << std::endl;
      inplace_qshift(PsiRight, QRight);

      // Check that we have common quantum number sectors in the two middle bases.
      if (PackMatrixOperator(PsiLeft.Basis1(), PsiRight.Basis1()).size() == 0)
      {
         std::cerr << "fatal: The left and right boundaries do not have any quantum number sectors in common. "
                      "You may need to set a different quantum number shift of the right boundary with --quantumnumber." << std::endl;
         return 1;
      }

      WavefunctionSectionLeft Window;

      // To join the wavefunctions we need a matrix to represent the
      // wavefunction in the left/right bipartite basis. We can do that by
      // minimizing the energy, or we can make a random state. We need a random
      // state anyway, as input for the Lanczos algorithm.

      // NOTE: This is between the Basis1 of PsiLeft and the Basis2 of
      // PsiRight, since we implicitly use the qshifted version of PsiLeft
      // here: once the wavefunction is loaded, PsiLeft can be qshifted explicitly.
      MatrixOperator C = MakeRandomMatrixOperator(PsiLeft.Basis1(), PsiRight.Basis1());
      C *= 1.0 / norm_frob(C);

      if (Random)
      {
         // We don't need to do anything here.
      }
      else if (BoundaryLambda)
      {
         if (PsiLeft.Basis1() != PsiRight.Basis1())
         {
            std::cerr << "fatal: The left and right boundaries must have the same middle basis if using --boundary. "
                         "You may need to set the quantum number shift of the right boundary explicitly with --quantumnumber." << std::endl;
            return 1;
         }
         // We use lambda_l because we implicitly qshift PsiLeft (compare above).
         C = PsiLeft.lambda_l();
      }
      else
      {
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
         BasicTriangularMPO HamMPO, HamMPOLeft, HamMPORight;
         std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

         HamMPOLeft = HamMPO;
         if (HamMPOLeft.size() < PsiLeft.size())
            HamMPOLeft = repeat(HamMPOLeft, PsiLeft.size() / HamMPOLeft.size());

         std::cout << "Solving fixed-point Hamiltonian..." << std::endl;
         StateComponent BlockHamL = Initial_E(HamMPOLeft, PsiLeft.Basis1());
         std::complex<double> LeftEnergy = SolveHamiltonianMPO_Left(BlockHamL, PsiLeft, HamMPOLeft, GMRESTol, Verbose);
         std::cout << "Starting energy (left eigenvalue) = " << LeftEnergy << std::endl;

         HamMPORight = HamMPO;
         if (HamMPORight.size() < PsiRight.size())
            HamMPORight = repeat(HamMPORight, PsiRight.size() / HamMPORight.size());

         StateComponent BlockHamR = Initial_F(HamMPORight, PsiRight.Basis2());
         std::complex<double> RightEnergy = SolveHamiltonianMPO_Right(BlockHamR, PsiRight, HamMPORight, GMRESTol, Verbose);
         std::cout << "Starting energy (right eigenvalue) = " << RightEnergy << std::endl;
         BlockHamR = delta_shift(BlockHamR, PsiRight.qshift());

         CMultiply Mult(BlockHamL, BlockHamR);
         std::vector<std::vector<std::complex<double>>> OutputVectors(NumEigen, std::vector<std::complex<double>>(Mult.P.size()));

         auto Eigs = LinearAlgebra::DiagonalizeARPACK(Mult, Mult.P.size(), NumEigen, LinearAlgebra::WhichEigenvalues::SmallestReal,
                                                      nullptr, Tol, OutputVectors.data(), 0, true, Verbose-1);

         C = Mult.P.unpack(OutputVectors[0].data());

         for (auto e : Eigs)
         {
            std::cout << "Energy = " << formatting::format_complex(remove_small_imag(e)) << std::endl;
         }
      }

      // Make the window from our centre matrix.
      Window = WavefunctionSectionLeft(C);

      inplace_qshift(PsiRight, adjoint(QRight));

      IBCWavefunction ResultPsi(PsiLeft, Window, PsiRight, PsiLeft.qshift(), QRight);

      if (Streaming)
      {
         ResultPsi.set_left_filename(LeftFilename);
         ResultPsi.set_right_filename(RightFilename);
      }

      MPWavefunction Result(ResultPsi);

      // Attributes
      Result.SetDefaultAttributes();

      // Set the Hamiltonian if we specified it.
      if (!HamStr.empty())
         Result.Attributes()["Hamiltonian"] = HamStr;

      // History log
      Result.AppendHistoryCommand(EscapeCommandline(argc, argv));

      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Result));
      pheap::ShutdownPersistent(P);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      if (e.Why == "File exists")
         std::cerr << "Note: Use --force (-f) option to overwrite." << std::endl;
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      return 1;
   }
}
