// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-reorder-symmetry.cpp
//
// Copyright (C) 2004-2022 Ian McCulloch <ian@qusim.net>
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

#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "interface/inittemp.h"
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>
#include "common/terminal.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

LinearWavefunction ReorderSymmetry(LinearWavefunction const& Psi, SymmetryList const& NewSL)
{
   LinearWavefunction Result(NewSL);
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result.push_front(CoerceSymmetryList(*I, NewSL));
   }
   return Result;
}

InfiniteWavefunctionLeft ReorderSymmetry(InfiniteWavefunctionLeft const& Psi, SymmetryList const& NewSL)
{
   InfiniteWavefunctionLeft Result;
   for (int i = 0; i < Psi.size(); ++i)
   {
      Result.push_back(CoerceSymmetryList(Psi[i], NewSL));
      Result.push_back_lambda(CoerceSymmetryList(Psi.lambda(i), NewSL));
   }

   Result.push_back_lambda(CoerceSymmetryList(Psi.lambda_r(), NewSL));

   Result.QShift = CoerceSymmetryList(Psi.qshift(), NewSL);

   Result.setBasis1(Result.lambda_l().Basis1());
   Result.setBasis2(Result.lambda_r().Basis2());
   Result.LogAmplitude = Psi.LogAmplitude;

   return Result;
}

InfiniteWavefunctionRight ReorderSymmetry(InfiniteWavefunctionRight const& Psi, SymmetryList const& NewSL)
{
   InfiniteWavefunctionRight Result;
   for (int i = 0; i < Psi.size(); ++i)
   {
      Result.push_back(CoerceSymmetryList(Psi[i], NewSL));
      Result.push_back_lambda(CoerceSymmetryList(Psi.lambda(i), NewSL));
   }

   Result.push_back_lambda(CoerceSymmetryList(Psi.lambda_r(), NewSL));

   Result.QShift = CoerceSymmetryList(Psi.qshift(), NewSL);

   Result.setBasis1(Result.lambda_l().Basis1());
   Result.setBasis2(Result.lambda_r().Basis2());
   //Result.LogAmplitude = Psi.LogAmplitude;

   return Result;
}

WavefunctionSectionLeft ReorderSymmetry(WavefunctionSectionLeft const& Psi, SymmetryList const& NewSL)
{
   WavefunctionSectionLeft Result;
   for (int i = 0; i < Psi.size(); ++i)
   {
      Result.push_back(CoerceSymmetryList(Psi[i], NewSL));
      Result.push_back_lambda(CoerceSymmetryList(Psi.lambda(i), NewSL));
   }

   Result.push_back_lambda(CoerceSymmetryList(Psi.lambda_r(), NewSL));

   Result.setBasis1(Result.lambda_l().Basis1());
   Result.setBasis2(Result.lambda_r().Basis2());
   return Result;
}

IBCWavefunction ReorderSymmetry(IBCWavefunction const& Psi, SymmetryList const& NewSL)
{
   return IBCWavefunction(ReorderSymmetry(Psi.left(), NewSL),
                          ReorderSymmetry(Psi.window(), NewSL),
                          ReorderSymmetry(Psi.right(), NewSL),
                          CoerceSymmetryList(Psi.left_qshift(), NewSL),
                          CoerceSymmetryList(Psi.right_qshift(), NewSL),
                          Psi.window_offset(),
                          Psi.window_left_sites(),
                          Psi.window_right_sites());
}

EAWavefunction ReorderSymmetry(EAWavefunction const& Psi, SymmetryList const& NewSL)
{
   std::vector<WavefunctionSectionLeft> WindowVec = Psi.window_vec();

   for (auto& Window : WindowVec)
      Window = ReorderSymmetry(Window, NewSL);

   return EAWavefunction(ReorderSymmetry(Psi.left(), NewSL),
                         WindowVec,
                         ReorderSymmetry(Psi.right(), NewSL),
                         CoerceSymmetryList(Psi.left_qshift(), NewSL),
                         CoerceSymmetryList(Psi.right_qshift(), NewSL),
                         Psi.left_index(),
                         Psi.right_index(),
                         Psi.exp_ik(),
                         Psi.gs_overlap());
}

FiniteWavefunctionLeft ReorderSymmetry(FiniteWavefunctionLeft const& Psi, SymmetryList const& NewSL)
{
   FiniteWavefunctionLeft Result;
   for (int i = 0; i < Psi.size(); ++i)
   {
      Result.push_back(CoerceSymmetryList(Psi[i], NewSL));
      Result.push_back_lambda(CoerceSymmetryList(Psi.lambda(i), NewSL));
   }

   Result.push_back_lambda(CoerceSymmetryList(Psi.lambda_r(), NewSL));

   Result.setBasis1(Result.lambda_l().Basis1());
   Result.setBasis2(Result.lambda_r().Basis2());
   return Result;
}

struct ApplyReorderSymmetry : public boost::static_visitor<WavefunctionTypes>
{
   ApplyReorderSymmetry(SymmetryList const& FinalSL_)
      : FinalSL(FinalSL_) {}

   template <typename T>
   T operator()(T const& Psi) const
   {
      return ReorderSymmetry(Psi, FinalSL);
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
         ("inpsi", prog_opt::value(&InputFile), "psi1")
         ("outpsi", prog_opt::value(&OutputFile), "psi1")
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
         std::cerr << "This program reorders the quantum numbers in a symmetry list.  New quantum numbers can\n"
            "be added, and they are set to the scalar quantum number.  Quantum numbers can also be removed, as long\n"
            "as they have degree 1 (ie, they are abelian, or in the abelian subset of a non-abelian symmetry).\n"
            "The symmetry-list is not allowed to be empty, so if removing all symmetries, use a symmetry-list of\n"
            "Null:Null\n";

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

      Result.Wavefunction() = boost::apply_visitor(ApplyReorderSymmetry(FinalSL), InputPsi->Wavefunction());
      Result.check_structure();

      Result.AppendHistoryCommand(EscapeCommandline(argc, argv));
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
