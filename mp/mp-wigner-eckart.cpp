// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-wigner-eckart.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

//
// Project a wavefunction using the Wigner-Eckart theorem
//
// The Regularize option is bugged - somehow the C_right and C_old
// matrices end up with incompatible bases, perhaps due to sorting of quantum numbers?

#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "tensor/wigner_eckart.h"
#include "tensor/regularize.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::Projection;

namespace prog_opt = boost::program_options;

StateComponent wigner_eckart(StateComponent const& A,
                             WignerEckartBasis<VectorBasis> const& W1,
                             WignerEckartBasis<VectorBasis> const& W2)
{
   // Get the projected local basis
   BasisList AbelianLocalBasis(W1.AbelianBasis().GetSymmetryList());
   for (unsigned i = 0; i < A.LocalBasis().size(); ++i)
   {
      QuantumNumbers::ProjectionList pl = enumerate_projections(A.LocalBasis()[i]);
      for (unsigned pi = 0; pi < pl.size(); ++pi)
      {
         AbelianLocalBasis.push_back(map_projection_to_quantum(pl[pi], AbelianLocalBasis.GetSymmetryList()));
      }
   }

   StateComponent Result(AbelianLocalBasis, W1.AbelianBasis(), W2.AbelianBasis());
   int k = 0;
   for (unsigned i = 0; i < A.LocalBasis().size(); ++i)
   {
      QuantumNumbers::ProjectionList pl = enumerate_projections(A.LocalBasis()[i]);
      for (unsigned pi = 0; pi < pl.size(); ++pi)
      {
         Result[k++] = wigner_eckart(A[i], pl[pi], W1, W2);
      }
   }
   return Result;
}

InfiniteWavefunctionLeft
wigner_project(InfiniteWavefunctionLeft const& Psi, SymmetryList const& FinalSL)
{
   InfiniteWavefunctionLeft Result;

   VectorBasis b1 = Psi.Basis1();
   WignerEckartBasis<VectorBasis> W2(b1, FinalSL);

   Result.setBasis1(W2.AbelianBasis());

   // get the identity projection
   QuantumNumbers::ProjectionList PL = enumerate_projections(Psi.lambda_l().TransformsAs());
   CHECK_EQUAL(PL.size(), 1U);
   Projection IdentP = PL[0];

   QuantumNumbers::ProjectionList QPL = enumerate_projections(Psi.qshift());
   Result.QShift = map_projection_to_quantum(QPL[0], FinalSL);

   for (int i = 0; i < Psi.size(); ++i)
   {

      StateComponent C = Psi[i];
      WignerEckartBasis<VectorBasis> W1 = W2;
      W2 = WignerEckartBasis<VectorBasis>(C.Basis2(), FinalSL);

      Result.push_back_lambda(wigner_eckart(Psi.lambda(i), IdentP, W1, W1));
      Result.push_back(wigner_eckart(C, W1, W2));
   }

   Result.push_back_lambda(wigner_eckart(Psi.lambda_r(), IdentP, W2, W2));

   Result.setBasis2(W2.AbelianBasis());
   Result.LogAmplitude = Psi.LogAmplitude;

   Result.check_structure();

   return Result;
}

InfiniteWavefunctionRight
wigner_project(InfiniteWavefunctionRight const& Psi, SymmetryList const& FinalSL)
{
   InfiniteWavefunctionRight Result;

   VectorBasis b1 = Psi.Basis1();
   WignerEckartBasis<VectorBasis> W2(b1, FinalSL);

   Result.setBasis1(W2.AbelianBasis());

   // get the identity projection
   QuantumNumbers::ProjectionList PL = enumerate_projections(Psi.lambda_l().TransformsAs());
   CHECK_EQUAL(PL.size(), 1U);
   Projection IdentP = PL[0];

   QuantumNumbers::ProjectionList QPL = enumerate_projections(Psi.qshift());
   Result.QShift = map_projection_to_quantum(QPL[0], FinalSL);

   for (int i = 0; i < Psi.size(); ++i)
   {

      StateComponent C = Psi[i];
      WignerEckartBasis<VectorBasis> W1 = W2;
      W2 = WignerEckartBasis<VectorBasis>(C.Basis2(), FinalSL);

      Result.push_back_lambda(wigner_eckart(Psi.lambda(i), IdentP, W1, W1));
      Result.push_back(wigner_eckart(C, W1, W2));
   }

   Result.push_back_lambda(wigner_eckart(Psi.lambda_r(), IdentP, W2, W2));

   Result.setBasis2(W2.AbelianBasis());
   //Result.LogAmplitude = Psi.LogAmplitude;

   Result.check_structure();

   return Result;
}

WavefunctionSectionLeft
wigner_project(WavefunctionSectionLeft const& Psi, SymmetryList const& FinalSL)
{
   PANIC("not implemented.");
}

IBCWavefunction
wigner_project(IBCWavefunction const& Psi, SymmetryList const& FinalSL)
{
   QuantumNumber LeftQShift = map_projection_to_quantum(enumerate_projections(Psi.left_qshift())[0], FinalSL);
   QuantumNumber RightQShift = map_projection_to_quantum(enumerate_projections(Psi.right_qshift())[0], FinalSL);

   return IBCWavefunction(wigner_project(Psi.left(), FinalSL),
                          wigner_project(Psi.window(), FinalSL),
                          wigner_project(Psi.right(), FinalSL),
                          LeftQShift,
                          RightQShift,
                          Psi.window_offset(),
                          Psi.window_left_sites(),
                          Psi.window_right_sites());
}

EAWavefunction
wigner_project(EAWavefunction const& Psi, SymmetryList const& FinalSL)
{
   QuantumNumber LeftQShift = map_projection_to_quantum(enumerate_projections(Psi.left_qshift())[0], FinalSL);
   QuantumNumber RightQShift = map_projection_to_quantum(enumerate_projections(Psi.right_qshift())[0], FinalSL);

   std::vector<WavefunctionSectionLeft> WindowVec = Psi.window_vec();

   for (auto& Window : WindowVec)
      Window = wigner_project(Window, FinalSL);

   return EAWavefunction(wigner_project(Psi.left(), FinalSL),
                         WindowVec,
                         wigner_project(Psi.right(), FinalSL),
                         LeftQShift,
                         RightQShift,
                         Psi.left_index(),
                         Psi.right_index(),
                         Psi.exp_ik(),
                         Psi.gs_overlap());
}

FiniteWavefunctionLeft
wigner_project(FiniteWavefunctionLeft const& Psi, SymmetryList const& FinalSL)
{
   // TODO: projection quantum number
   FiniteWavefunctionLeft Result;

   VectorBasis b1 = Psi.Basis1();
   WignerEckartBasis<VectorBasis> W2(b1, FinalSL);

   Result.setBasis1(W2.AbelianBasis());

   // get the identity projection
   QuantumNumbers::ProjectionList PL = enumerate_projections(Psi.lambda_l().TransformsAs());
   CHECK_EQUAL(PL.size(), 1U);
   Projection IdentP = PL[0];

   for (int i = 0; i < Psi.size(); ++i)
   {

      StateComponent C = Psi[i];
      WignerEckartBasis<VectorBasis> W1 = W2;
      W2 = WignerEckartBasis<VectorBasis>(C.Basis2(), FinalSL);

      Result.push_back_lambda(wigner_eckart(Psi.lambda(i), IdentP, W1, W1));
      Result.push_back(wigner_eckart(C, W1, W2));
   }

   Result.push_back_lambda(wigner_eckart(Psi.lambda_r(), IdentP, W2, W2));

   Result.setBasis2(W2.AbelianBasis());

   Result.check_structure();

   return Result;
}

// functor to use the visitor pattern with wavefunction types
struct ApplyWignerEckart : public boost::static_visitor<WavefunctionTypes>
{
   ApplyWignerEckart(SymmetryList const& FinalSL_)
      : FinalSL(FinalSL_) {}

   template <typename T>
   T operator()(T const& Psi) const
   {
      return wigner_project(Psi, FinalSL);
   }

   SymmetryList FinalSL;
};

int main(int argc, char** argv)
{
   try
   {
      bool Force = false;
      std::string SList;
      std::string InputFile;
      std::string OutputFile;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("force,f", prog_opt::bool_switch(&Force),
          "overwrite the output file, if it exists")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("symmetrylist", prog_opt::value(&SList), "slist")
         ("inpsi", prog_opt::value(&InputFile), "inpsi")
         ("outpsi", prog_opt::value(&OutputFile), "outpsi")
         ;

      prog_opt::positional_options_description p;
      p.add("symmetrylist", 1);
      p.add("inpsi", 1);
      p.add("outpsi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("inpsi") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <symmetry-list> <input-psi> [output-psi]\n";
         std::cerr << desc << '\n';
         std::cerr << "This tool maps a non-abelian symmetry into a set of abelian projections.\n"
            "The only projection currently defined is from SU(2) to U(1), and partial projections are\n"
            "not allowed - that is, if there is more than one SU(2) symmetry, then ALL of them must be\n"
            "projected.  The final symmetry list must be idential to the original symmetry list, but\n"
            "with each SU(2) symmetry replaced by a U(1) symmetry, in the same order (the name doesn't matter).\n"
            ;
         return 1;
      }

      pvalue_ptr<MPWavefunction> InputPsi;

      if (OutputFile.empty())
      {
         // re-use the input file as the output file
         InputPsi = pheap::OpenPersistent(InputFile, mp_pheap::CacheSize());
      }
      else
      {
         // create a new file for output
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         // and load the input wavefunction
         InputPsi = pheap::ImportHeap(InputFile);
      }

      SymmetryList FinalSL = SymmetryList(SList);

      // If we are overwriting the old file, copy the old history and attributes
      MPWavefunction Result;
      if (OutputFile.empty())
      {
         Result = MPWavefunction(InputPsi->Attributes(), InputPsi->History());
      }
      Result.AppendHistoryCommand(EscapeCommandline(argc, argv));
      Result.Wavefunction() = boost::apply_visitor(ApplyWignerEckart(FinalSL), InputPsi->Wavefunction());
      Result.SetDefaultAttributes();

      pvalue_ptr<MPWavefunction> OutputPsi = new MPWavefunction(Result);

      pheap::ShutdownPersistent(OutputPsi);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
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
