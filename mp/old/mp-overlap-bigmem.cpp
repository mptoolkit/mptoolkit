// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-overlap-bigmem.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/environment.h"
#include <fstream>
#include <sys/types.h>
#include <sys/wait.h>

std::vector<std::string> Files;
std::vector<MPWavefunction> Wavefunctions;

std::vector<std::pair<int, int> > PairsNeeded;

int main(int argc, char** argv)
{
   if (argc < 2 || (argc == 2 && argv[1] != std::string("-")))
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-overlap <psi1> <psi2> [<psi3> ...]\n";
      std::cerr << "or use '-' as the first filename to read files from standard input.\n";
      return 1;
   }

   if (argc == 2 && argv[1] == std::string("-"))
   {
      std::string s;
      while (std::cin >> s)
         Files.push_back(s);
   }
   else
   {
      std::copy(&argv[1], &argv[argc], std::back_inserter(Files));
   }

   // first file
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> p = pheap::OpenPersistent(Files[0], CacheSize, true);
   Wavefunctions.push_back(*p);

   for (unsigned i = 1; i < Files.size(); ++i)
   {
      p = pheap::ImportHeap(Files[i]);
      Wavefunctions.push_back(*p);
   }
   p = pvalue_ptr<MPWavefunction>();

   for (unsigned i = 0; i < Files.size(); ++i)
   {
      for (unsigned j = i+1; j < Files.size(); ++j)
      {
         PairsNeeded.push_back(std::make_pair(i,j));
      }
   }

   int NumThreads = getenv_or_default("MP_NUMTHREADS", 1);
   {
      int Offset = 0;
      std::vector<pid_t> ForkList;
      for (int i = 0; i < NumThreads; ++i)
      {
         pid_t pid = fork();
         if (pid == 0)
         {
            std::ofstream out(("file" + boost::lexical_cast<std::string>(Offset)).c_str());
            out.precision(14);
            while (Offset < int(PairsNeeded.size()))
            {
               int x = PairsNeeded[Offset].first;
               int y = PairsNeeded[Offset].second;
               std::complex<double> val = overlap(Wavefunctions[x], Wavefunctions[y]);
               out << Files[x] << ' ' << Files[y] << ' ' << val << std::endl;
               Offset += NumThreads;
            }
            out.flush();
            _exit(0);
         }
         else
            ForkList.push_back(pid);

         ++Offset;
      }

      for (int i = 0; i < NumThreads; ++i)
      {
         int status;
         waitpid(ForkList[i], &status, 0);
      }
   }

   pheap::Shutdown();
}
