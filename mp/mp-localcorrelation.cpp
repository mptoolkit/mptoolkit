// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-localcorrelation.cpp
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
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"
#include "interface/inittemp.h"
#include "common/proccontrol.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string StringOp = "I";
      bool ShowRealPart = false, ShowImagPart = false;
      std::string LatticeStr, PsiStr, Op1, Op2;
      int FirstSite, LastSite;
      bool IncludeFirst, Fermionic;
      bool Quiet = false;


      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("real,r", prog_opt::bool_switch(&ShowRealPart),
          "display only the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImagPart),
          "display only the imaginary part of the result")
         ("string,s", prog_opt::value(&StringOp),
          "calculate a string correlation with the given operator S, as operator1 "
          "\\otimes S \\otimes S .... S \\otimes operator2")
         ("includefirst,d", prog_opt::bool_switch(&IncludeFirst),
          "when calculating a string correlation, apply the string operator also "
          "to the first site, as operator1*S")
         ("fermionic,f", prog_opt::bool_switch(&Fermionic),
          "take the correlator to be fermionic; equivalent to \"--string P --includefirst\"")
         ("quiet,q", prog_opt::bool_switch(&Quiet),
          "suppress warning messages")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lat", prog_opt::value(&LatticeStr), "lat")
         ("psi", prog_opt::value(&PsiStr), "psi")
         ("op1", prog_opt::value(&Op1), "op1")
         ("first", prog_opt::value(&FirstSite), "first")
         ("op2", prog_opt::value(&Op2), "op2")
         ("last", prog_opt::value(&LastSite), "last")
         ;

      prog_opt::positional_options_description p;
      p.add("lat", 1);
      p.add("psi", 1);
      p.add("op1", 1);
      p.add("first", 1);
      p.add("op2", 1);
      p.add("last", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (Fermionic)
      {
         if (IncludeFirst || StringOp != "I")
         {
            std::cerr << "mp-localcorrelation: error: cannot combine string correlators with --fermionic.\n\n";
            return 1;
         }
         else
         {
            // a fermionic operator is equivalent to the string correlator
            // operator1*P \otimes P \otimes .... \otimes P \otimes operator2
            // which is equivalent to "--includefirst --string P" options.
            IncludeFirst = true;
            StringOp = "P";  // TODO: this should actually be the SignOperator() of the site operator at Op1
         }
      }

      if (vm.count("help") > 0 || vm.count("last") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-localcorrelation [options] <lattice> <psi> <operator1> <first> <operator2> <last>\n";
         std::cerr << "This uses directly local operators, not lattice operators!\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // to reproduce the old behaviour, if both ShowRealPart and ShowImagPart are false
      // we need to show both components.  Other tools would revert to showing the complex
      // number in standard C++ style (real,imag), but we never did this in the past for
      // mp-localcorrelation (we always used two columns for the real and imag, the
      // same format as if both --real and --imag were specified), and lets not change it now.
      if (!ShowRealPart && !ShowImagPart)
         ShowRealPart = ShowImagPart = true;

      bool FermionicWarning = false;  // set to true if we have already warned the user about
      // a possible fermionic problem

      mp_pheap::InitializeTempPHeap();

      pvalue_ptr<OperatorList> System = pheap::ImportHeap(LatticeStr);
      pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(PsiStr);

      Lattice Lat = System->GetLattice();

      CenterWavefunction Psi = *Psi1;
      Psi1 = pvalue_ptr<MPWavefunction>();

      for (int i = 1; i < FirstSite; ++i)
         Psi.RotateRight();

      std::cout.precision(12);

      // This holds the E matrices for the left system
      typedef std::map<int, pvalue_handle<MatrixOperator> > OpMapType;
      OpMapType OpMap;
      for (int i = FirstSite; i < LastSite; ++i)
      {
         // Update all the existing E matrices
         SiteBlock::const_iterator ThisStringOp = Lat[i].find(StringOp);
         if (ThisStringOp == Lat[i].end())
         {
            if (!Quiet)
               std::cerr << "mp-localcorrelation: warning: string operator " << StringOp
                         << " doesn't exist at site number " << i << ", assuming identity operator instead.\n";

            ThisStringOp = Lat[i].find("I");
            CHECK(ThisStringOp != Lat[i].end());  // the identity operator should always exist
         }

         SimpleOperator ThisString = ThisStringOp->second;

         for (OpMapType::iterator mI = OpMap.begin(); mI != OpMap.end(); ++mI)
         {
            mI->second = new MatrixOperator(operator_prod(herm(ThisString),
                                                          herm(Psi.Left()),
                                                          *mI->second.load(),
                                                          Psi.Left()));
         }

         // Add another E matrix for the new site, only if the left operator exists in the block
         SiteBlock::const_iterator I = Lat[i].find(Op1);
         if (I != Lat[i].end())
         {
            // we can do a check here that the operator is a fermionic operator, the user has
            // specified --fermionic too.  Otherwise, the result will be a fermionic string
            // correlation (which is still useful, but possibly not what the user wanted).
            if (I->second.Commute().Value_ == LatticeCommute::Fermionic
                && !Fermionic && !FermionicWarning)
            {
               FermionicWarning = true;
               if (!Quiet)
                  std::cerr << "mp-localcorrelation: warning: operator1 is fermionic, but --fermionic is not"
                     " specified.\n";
            }

            SimpleOperator MyOp = I->second;
            if (IncludeFirst)
               MyOp = MyOp * ThisString;// if this causes an error in the non-abelian case, it is user's fault;-)


            MatrixOperator LeftIdentity = MatrixOperator::make_identity(Psi.Left().Basis1());
            OpMap[i] = new MatrixOperator(operator_prod(herm(MyOp),
                                                        herm(Psi.Left()),
                                                        LeftIdentity,
                                                        Psi.Left()));
         }

         // For the right operator, construct the F matrix and the expectation value
         I = Lat[i+1].find(Op2);
         if (I != Lat[i+1].end())
         {
            SimpleOperator MyOp = I->second;
            MatrixOperator Ident = MatrixOperator::make_identity(Psi.Right().Basis2());
            MatrixOperator F = operator_prod(MyOp, Psi.Right(), Ident, herm(Psi.Right()));
            for (OpMapType::iterator mI = OpMap.begin(); mI != OpMap.end(); ++mI)
            {
               MatrixOperator E = *mI->second.load();
               std::complex<double> Res;
               if (!is_transform_target(adjoint(F.TransformsAs()), E.TransformsAs(), Psi.Center().TransformsAs()))
               {
                  // the quantum number's don't match up, the correlation is identically zero.
                  Res = 0.0;
               }
               else
               {
                  Res = inner_prod(Psi.Center(),
                                   triple_prod(*mI->second.load(),
                                               Psi.Center(),
                                               herm(F)));
               }
               std::cout << std::setw(5) << Lat.coordinate_at_site(mI->first) << "   "
                         << std::setw(5) << Lat.coordinate_at_site(i+1);
               if (ShowRealPart)
                  std::cout << "   " << std::setw(18) << Res.real();
               if (ShowImagPart)
                  std::cout << "   " << std::setw(18) << Res.imag();
               std::cout << '\n';
            }
         }

         if (i != LastSite-1)
            Psi.RotateRight();
      }

      pheap::Shutdown();

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
