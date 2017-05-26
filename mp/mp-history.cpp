// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-history.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "quantumnumbers/all_symmetries.h"
#include "wavefunction/mpwavefunction.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "common/terminal.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      bool Reverse = false;
      std::string Filename;
      std::string Message;
      std::string Command;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("reverse,r", prog_opt::bool_switch(&Reverse), "reverse order, newest first")
         ("message,m", prog_opt::value(&Message), "add a new history entry as a note")
         ("command,c", prog_opt::value(&Command), "add a new history entry as a command")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&Filename), "psi")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", "mp-history");
         std::cerr << "usage: " << basename(argv[0]) << " [options] <wavefunction>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (vm.count("message") || vm.count("command"))
      {
         pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(Filename, mp_pheap::CacheSize());
         if (!Message.empty())
            Psi.mutate()->AppendHistoryNote(Message);
         if (!Command.empty())
            Psi.mutate()->AppendHistoryNote(Command);
         pheap::ShutdownPersistent(Psi);
      }
      else
      {
         pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(Filename, mp_pheap::CacheSize(), true);
         if (Reverse)
         {
            Psi->History().print_newest_first(std::cout);
         }
         else
         {
            Psi->History().print(std::cout);
         }
      }
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
