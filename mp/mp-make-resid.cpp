// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-make-resid.cpp
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
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

int main(int argc, char** argv)
{
   try
   {

      double GroundstateEnergy = 0;
      double Frequency = 0;
      double Broadening = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(),
          "operator to use for the Hamiltonian"
          " (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "initial correction vector wavefunction (required)")
         ("lanczos,l", prog_opt::value<std::string>(),
          "Lanczos vector for the right hand side (required)")
         ("out,o", prog_opt::value<std::string>(),
          "output filename for the residual (required)")
         ("GroundstateEnergy,G", prog_opt::value(&GroundstateEnergy),
          "groundstate energy of the Hamiltonian"
          " (wavefunction attribute \"GroundstateEnergy\")")
         ("Frequency,F", prog_opt::value(&Frequency),
          "frequency of the correction vector (wavefunction attribute \"Frequency\")")
         ("Broadening,B", prog_opt::value(&Broadening),
          "broadening constant eta (wavefunction attribute \"Broadening\")")
         ("broadening-propto-energy", prog_opt::bool_switch(),
          "use eta*(H-E) as the broadening term")
         ("verbose", prog_opt::bool_switch(&Verbose),

         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || !vm.count("wavefunction") || !vm.count("out")
          || !vm.count("lanczos"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-make-resid [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << '\n';

      std::string LanczosVector = vm["lanczos"].as<std::string>();
      std::cout << "Lanczos vector: " << LanczosVector << '\n';


   if (argc < 5 || argc > 7)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-make-resid\n";
      return 1;
   }

   std::string InPsi = argv[3];
   std::string OutPsi = argv[4];

   pvalue_ptr<MPWavefunction> Psi;

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   if (InPsi == OutPsi)
      Psi = pheap::OpenPersistent(InPsi.c_str(), CacheSize);
   else
   {
      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      pheap::Initialize(OutPsi, 1, PageSize, CacheSize);
      Psi = pheap::ImportHeap(InPsi.c_str());
   }

   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string Operator = argv[2];

   MPOperator Op = (*System)[Operator];

   int NStates = DefaultMaxStates;
   if (argc >= 6)
      NStates = boost::lexical_cast<int>(argv[5]);

   QuantumNumbers::QuantumNumber Q;
   if (argc >= 7)
   {
      Q = QuantumNumbers::QuantumNumber(Op.GetSymmetryList(), argv[6]);
   }
   else
   {
      QuantumNumbers::QuantumNumberList QList = transform_targets(Op.TransformsAs(), Psi->TransformsAs());
      if (QList.size() != 1)
      {
         std::cerr << "quantum number of operator-wavefunction product is ambiguous,\n"
                   << "  choices are ";
         std::copy(QList.begin(), QList.end(), std::ostream_iterator<QuantumNumber>(std::cerr, " "));
         std::cerr << '\n';
         exit(1);
      }
      Q = QList[0];
   }

   *Psi.mutate() = prod(Op, *Psi, Q, NStates);

   pheap::ShutdownPersistent(Psi);
}
