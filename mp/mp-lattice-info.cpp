// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-lattice-info.cpp
//
// Copyright (C) 2015-2020 Ian McCulloch <ian@qusim.net>
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

#include "lattice/infinitelattice.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"

namespace prog_opt = boost::program_options;

void DescribeLattice(InfiniteLattice const& L, int Verbose)
{
   std::cout << "Description: " << L.description() << '\n';
   for (auto const& x : L.authors())
   {
      std::cout << "Author: " << x.first << " <" << x.second << ">\n";
   }
   std::cout << "Command line: " << L.command_line() << '\n';
   std::cout << "Date: " << L.timestamp() << '\n';
   std::cout << "SymmetryList: " << L.GetSymmetryList() << '\n';
   std::cout << "Unit cell size: " << L.GetUnitCell().size() << '\n';
   for (int i = 0; i < L.GetUnitCell().size(); ++i)
   {
      std::cout << "   site [" << i << "] is: "
                << L.GetUnitCell()[i].Description() << '\n';
      if (Verbose > 0)
      {
         // show details of the lattice site
         std::cout << "      basis states:\n";
         for (unsigned n = 0; n < L.GetUnitCell()[i].Basis1().size(); ++n)
         {
            std::string State = "|" + L.GetUnitCell()[i].Basis1()[n].first + '>';
            std::cout << "         " << std::setw(10) << std::left << State
                      << " transforms as " << L.GetUnitCell()[i].Basis1()[n].second << '\n';
         }
         std::cout << "      local constants:\n";
         if (L.GetUnitCell()[i].arg_empty())
            std::cout << "         (none)\n";
         for (LatticeSite::const_argument_iterator I = L.GetUnitCell()[i].begin_arg();
              I != L.GetUnitCell()[i].end_arg(); ++I)
         {
            std::cout << "         " << std::setw(10) << std::left << I->first
                      << " = " << formatting::format_complex(I->second) << '\n';
         }
         std::cout << "      local operators:\n";
         for (LatticeSite::const_operator_iterator I = L.GetUnitCell()[i].begin_operator();
              I != L.GetUnitCell()[i].end_operator(); ++I)
         {
            std::cout << "         " << std::setw(10) << std::left << I->first
                      << " - " << I->second.description_or_none() << '\n';
            if (Verbose > 1)
            {
               // show the actual operator
               std::cout << "                      = " << I->second << '\n';
            }
         }
         std::cout << "      local functions:\n";
         if (L.GetUnitCell()[i].function_empty())
            std::cout << "         (none)\n";
         for (LatticeSite::const_function_iterator I = L.GetUnitCell()[i].begin_function();
              I != L.GetUnitCell()[i].end_function(); ++I)
         {
            std::cout << "         " << I->first << I->second << '\n';
         }
      }
   }

   std::cout << "\nUnit cell operators:\n";
   for (UnitCell::const_operator_iterator I = L.GetUnitCell().begin_operator();
        I != L.GetUnitCell().end_operator(); ++I)
   {
      std::cout << "   " << std::setw(10) << std::left << I->first
                << " - " << I->second.description() << '\n';
      if (Verbose > 0)
      {
         std::cout << "                -transforms: "
                   << std::setw(10) << std::left << I->second.TransformsAs()
                   << "  commutes: " << I->second.Commute() << '\n';
      }
   }

   std::cout << "\nUnit Cell functions:\n";
   if (L.GetUnitCell().function_empty())
      std::cout << "   (none)\n";
   for (UnitCell::const_function_iterator I = L.GetUnitCell().begin_function();
        I != L.GetUnitCell().end_function(); ++I)
   {
      if (!I->second.description().empty())
         std::cout << "   " << std::setw(10) << std::left << I->first << " - "
                   << I->second.description() << '\n';

      std::cout << "   " << I->first << I->second << '\n';
   }

   std::cout << "\nLattice constants:\n";
   if (L.arg_empty())
      std::cout << "   (none)\n";
   for (InfiniteLattice::const_argument_iterator I = L.begin_arg();
        I != L.end_arg(); ++I)
   {
      std::cout << "   " << std::setw(10) << std::left << I->first
                << " = " << formatting::format_complex(I->second) << '\n';
   }

   std::cout << "\nLattice operators:\n";
   for (InfiniteLattice::const_operator_iterator I = L.begin_operator();
        I != L.end_operator(); ++I)
   {
      std::cout << "   " << std::setw(10) << std::left << I->first
                << " - " << I->second.description() << '\n';
   }

   std::cout << "\nLattice functions:\n";
   if (L.function_empty())
      std::cout << "   (none)\n";
   for (InfiniteLattice::const_function_iterator I = L.begin_function();
        I != L.end_function(); ++I)
   {
      if (!I->second.description().empty())
         std::cout << "   " << std::setw(10) << std::left << I->first << " - "
                   << I->second.description() << '\n';

      std::cout << "   " << I->first << I->second << '\n';
   }
}

int main(int argc, char** argv)
{
   try
   {
      std::string LatticeName;
      int Verbose = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "show extra detail (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lattice", prog_opt::value<std::string>(&LatticeName), "lattice")
         ;
      prog_opt::positional_options_description p;
      p.add("lattice", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);
      if (vm.count("help") || !vm.count("lattice"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <lattice>\n";
         std::cerr << desc << "\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<InfiniteLattice> Lattice
         = pheap::OpenPersistent(LatticeName,
                                 mp_pheap::CacheSize(), true);

      DescribeLattice(*Lattice, Verbose);

   }
   catch(std::exception& e)
   {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
   }
   catch(...)
   {
      std::cerr << "Exception of unknown type!\n";
   }
   return 0;
}
