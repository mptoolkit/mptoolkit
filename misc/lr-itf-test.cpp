// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/lr-itf-test.cpp
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// Test program for LR-ITF model

#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <vector>


int main(int argc, char** argv)
{
   if (argc != 4)
   {
      std::cerr << "usage: lr-itf-test <file> <nlegs> <length>\n";
      std::cout << "File format is i j Alpha Beta r2\n"
                << "for an interaction alpha*Sz(i)*Sz(j) (with quirks!)\n";
      return 1;
   }

   std::ifstream In(argv[1]);
   int NLegs = boost::lexical_cast<int>(argv[2]);

   int const N = boost::lexical_cast<int>(argv[3]);
   std::vector<std::vector<double> > SiteEnergy(NLegs, std::vector<double>(N, 0.0));
   int i, j, r2;
   double Alpha, Beta;
   double Energy = 0;
   while (In >> i >> j >> Alpha >> Beta >> r2)
   {
      std::cout << "i=" << i << " j=" << j << " Alpha=" << Alpha << " Beta=" << Beta << " r2=" << r2
                << '\n';
      if (j < NLegs)
      {
         std::cout << "Found local term " << Alpha << " * Sz(" << i << ") * Sz(" << j << ")\n";
         std::cout << "Ferromagnetic energy contribution is " << -Alpha << "\n";
         Energy += Alpha;
         SiteEnergy[i][j] += Alpha;
      }
      else
      {
         // sum 1 + beta + beta^2 + ....
         double Sum = 1.0 / (1.0 - Beta);
         if (j-NLegs > i)
         {
            std::cout << "Found a term that will act on the first unit cell also.\n";
            std::cout << "Ferromagnetic energy contribution is " << -Alpha*Sum << "\n";
            Energy += Alpha*Sum;
            double E = Alpha;
            for (int k = j-NLegs; k < N; k += NLegs)
            {
               SiteEnergy[i][k] += E;
               E *= Beta;
            }
         }
         else
         {
            std::cout << "Ferromagnetic energy contribution is " << -Alpha*Sum << "\n";
            Energy += Alpha*Sum;
            double E = Alpha;
            for (int k = j; k < N; k += NLegs)
            {
               SiteEnergy[i][k] += E;
               E *= Beta;
            }
         }
      }
   }

   std::cout << "\nDone.  Ferromagnetic energy is " << Energy << "\n";
   std::cout << "Energy contribution in the unit cell is\n"
             << "#i  #j    #E\n";
   std::cout.precision(12);
   for (i = 0; i < NLegs; ++i)
   {
      for (j = 0; j < N; ++j)
      {
         std::cout << i << " " << j << " " << std::setw(16) << SiteEnergy[i][j] << '\n';
      }
   }
}
