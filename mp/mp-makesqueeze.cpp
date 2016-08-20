// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-makesqueeze.cpp
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
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "common/math_const.h"
#include "pheap/pheap.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <cmath>
#include "common/environment.h"

namespace prog_opt = boost::program_options;
using std::sin;

// the int count here is the number of filled components
std::map<int, MPOperator>
AddSqueezedPart(std::map<int, MPOperator> const& In, MPOperator const& Filled, MPOperator const& Empty,
                int MaxFilled)
{
   std::map<int, MPOperator> Result;
   for (std::map<int, MPOperator>::const_iterator I = In.begin(); I != In.end(); ++I)
   {
      if (I->first <= MaxFilled)
         Result[I->first] += prod(I->second, Empty, I->second.TransformsAs());
      if (I->first+1 <= MaxFilled)
         Result[I->first+1] += prod(I->second, Filled, I->second.TransformsAs());
   }
   return Result;
}

std::vector<MPOperator>
ConstructSqueezedCorrelator(MPOperator Op1,
                            int FirstSite,
                            int Distance,
                            int MaxLag,
                            OperatorAtSite<OperatorList const, int> const& OpProjFilled,
                            OperatorAtSite<OperatorList const, int> const& OpProjEmpty,
                            OperatorAtSite<OperatorList const, int> const& Op2)
{
   std::vector<MPOperator> Result(Distance);
   std::map<int, MPOperator> BasicSet;
   Result[0] = prod(Op1, Op2(FirstSite+1), QuantumNumber(Op1.GetSymmetryList()));
   BasicSet[0] = Op1;
   for (int i = 0; i < Distance+MaxLag; ++i)
   {
      std::cout << "Distance " << (i+1) << std::endl;
      BasicSet = AddSqueezedPart(BasicSet, OpProjFilled(FirstSite+i+1), OpProjEmpty(FirstSite+i+1),
                                 Distance-1);
      for (std::map<int, MPOperator>::const_iterator I = BasicSet.begin(); I != BasicSet.end(); ++I)
      {
         Result[I->first] += prod(I->second, Op2(FirstSite+i+2), QuantumNumber(Op1.GetSymmetryList()));
      }
   }
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      std::string LatticeFile;
      std::string Operator1;
      std::string Operator2;
      std::string ProjFilled;
      std::string ProjEmpty;
      std::string OutName;
      int FirstSite;
      int Length;
      int MaxLag;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("lattice,l", prog_opt::value(&LatticeFile),
          "lattice file [required]")
         ("operator1,1", prog_opt::value(&Operator1),
          "name of the left operator [required]")
         ("operator2,2", prog_opt::value(&Operator2),
          "name of the right operator [required]")
         ("filled-proj,p", prog_opt::value(&ProjFilled),
          "name of the projector onto a filled site [required]")
         ("empty-proj,e", prog_opt::value(&ProjEmpty),
          "name of the projector onto an empty site [required]")
         ("out,o", prog_opt::value(&OutName),
          "name of the output operator [required]")
         ("first-site,f", prog_opt::value(&FirstSite),
          "first site for the correlators")
         ("length,g", prog_opt::value(&Length),
          "length to calculate the correlator over")
         ("lag,m", prog_opt::value(&MaxLag),
          "extra length permissible for the squeezed space")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || !vm.count("lattice"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-make-k [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // load the lattice file
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<OperatorList> System = pheap::OpenPersistent(LatticeFile, CacheSize);
      // get the lattice size
      //int LatticeSize = System->size();

      OperatorAtSite<OperatorList const, int> Op1(*System, Operator1);
      OperatorAtSite<OperatorList const, int> Op2(*System, Operator2);
      OperatorAtSite<OperatorList const, int> OpProjFilled(*System, ProjFilled);
      OperatorAtSite<OperatorList const, int> OpProjEmpty(*System, ProjEmpty);
      OperatorAtSite<OperatorList, int> OpOut(*System.mutate(), OutName);

      MPOperator Base = Op1(FirstSite);

      std::vector<MPOperator> Operators = ConstructSqueezedCorrelator(Op1(FirstSite),
                                                                      FirstSite, Length, MaxLag,
                                                                      OpProjFilled, OpProjEmpty, Op2);

      std::cout << "Adding " << OutName << "(0)" << std::endl;
      OpOut(0) = prod(Op1(FirstSite), Op2(FirstSite), QuantumNumber(Op1(FirstSite).GetSymmetryList()));
      for (unsigned k = 0; k < Operators.size(); ++k)
      {
         std::cout << "Adding " << OutName << "(" << (k+1) << ")" << std::endl;
         OpOut(k+1) = Operators[k];
      }

      // save the lattice
      pheap::ShutdownPersistent(System);

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
