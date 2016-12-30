// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-imake-transverse-lattice.cpp
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

#include "mpo/triangularoperator.h"
#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "mp/copyright.h"
#include "common/proccontrol.h"

#include "models/spin-su2.h"
#include "models/spin-u1.h"
#include "models/spin-z2.h"
#include "models/spin.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"
#include "models/hubbard-u1u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/bosehubbard-spinless.h"
#include "models/bosehubbard-spinless-u1.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "tensor/tensor_exponential.h"
#include "mpo/mpoperator.h"

#include "interface/inittemp.h"

#include "mpo/lattice.h"
#include "mpo/mpoperatorlist.h"

#include <boost/iterator/transform_iterator.hpp>

namespace prog_opt = boost::program_options;

MPOpComponent ConstructFromOperatorComponent(OperatorComponent const& x)
{
   MPOpComponent Result(x.LocalBasis1(), x.Basis1(), x.Basis2());
   for (OperatorComponent::const_iterator I = iterate(x); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
         for (SimpleRedOperator::const_iterator S = J->begin(); S != J->end(); ++S)
         {
            Result.set_operator(J.index1(), J.index2(), *S);
         }
      }
   }
   return Result;
}

LinearOperator ConstructFromMPOperator(MPOperator const& Op)
{
   std::vector<MPOpComponent> MyOp(boost::make_transform_iterator(Op.begin(), &ConstructFromOperatorComponent),
                                   boost::make_transform_iterator(Op.end(), &ConstructFromOperatorComponent));

   return MPOpCompressed(MyOp.begin(), MyOp.end());
}

SiteBlock MakeSiteBlockFromBasis(BasisList const& BL)
{
   SiteBlock Block;
   SiteBasis MyBasis(BL.GetSymmetryList());
   for (unsigned i = 0; i < BL.size(); ++i)
   {
      MyBasis.push_back(boost::lexical_cast<std::string>(i), BL[i]);
   }
   SiteOperator I = SiteOperator::Identity(MyBasis);
   Block["I"] = I;
   return Block;
}

int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true);
   try
   {
      // Timestep options
      double Timestep = 0.01;
      double NumTimesteps = 20;

      // Hamiltonian options
      std::string HamStr;
      double Lambda = 1.0;
      double J = 1.0;
      double J2 = 0.0;
      double B = 0;
      double D = 0;
      double U = 0;
      double Jz = 0.0;
      double tprime = 1.0;
      double delta = 0.0;
      double Theta = 0.0;
      double Beta = 0.0;
      half_int Spin = 0.5;
      int NLegs = 1;

      // wavefunction options
      std::string FName;

      std::string LatticeFile;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian.  Valid choices: itf, itf-z2, xxx-su2, xxx-u1, xxx, tj-zigzag-u1su2, "
          "tj-zigzag-u1, sf-zigzag-u1, klm-u1su2, klm-u1, bh, bh2, bh-u1, bh2-u1")
         ("wavefunction,w", prog_opt::value(&FName),
          "Initial wavefunction (required)")
         ("lattice,l", prog_opt::value(&LatticeFile),
          "Lattice file to construct (required)")
         ("timestep,t", prog_opt::value(&Timestep),
          "Timestep for 2nd order S-T evolution")
         ("numsteps,n", prog_opt::value(&NumTimesteps),
          "Number of timesteps to perform")
         ("spin", prog_opt::value(&Spin),
          FormatDefault("spin (for xxx,xxz,xyz hamiltonians)", Spin).c_str())
         ("J", prog_opt::value(&J),
          FormatDefault("nearest-neighbor exchange J (for xxx,itf, etc)", J).c_str())
         ("J2", prog_opt::value(&J2),
          FormatDefault("next-nearest-neighbor exchange J2 (for xxx)", J2).c_str())
         ("D", prog_opt::value(&D),
          FormatDefault("single-ion anisotropy (for xxx-u1 and xxx)", D).c_str())
         ("U", prog_opt::value(&U),
          FormatDefault("coulomb repulsion", U).c_str())
         ("B", prog_opt::value(&B),
          FormatDefault("magnetic field (for xxx)", B).c_str())
         ("Jz", prog_opt::value(&Jz),
          FormatDefault("Jz coupling (for Kondo)", Jz).c_str())
         ("nlegs", prog_opt::value(&NLegs),
          FormatDefault("Number of legs (for triangular ladder)", NLegs).c_str())
         ("tprime", prog_opt::value(&tprime),
          FormatDefault("next-nearest-neighbor hopping t' (for tj-zigzag, sf-zigzag)", tprime).c_str())
         ("delta", prog_opt::value(&delta),
          FormatDefault("Zigzag ladder potential imbalance (for tj-zigzag, sf-zigzag)", delta).c_str())
         ("theta", prog_opt::value(&Theta),
          FormatDefault("theta (for biquadratic xxx)", Theta).c_str())
         ("Beta", prog_opt::value(&Beta),
          FormatDefault("Beta (for biquadratic xxx)", Beta).c_str())
         ("lambda", prog_opt::value(&Lambda),
          FormatDefault("transverse field strength (for itf hamiltonian)", Lambda).c_str())
          ;


      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || HamStr.empty())
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-imake-transverse-lattice [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pheap::Initialize(LatticeFile, 1, PageSize, CacheSize);

      SimpleOperator H;
      SimpleOperator Obs;
      BasisList B1, B2;
      if (HamStr == "itf")
      {
         std::cout << "Hamiltonian is transverse-field Ising, J=" << J << ", Lambda=" << Lambda << "\n";
         SiteBlock Site = CreateSpinSite(0.5);
         H = J * 4.0 * tensor_prod(Site["Sz"], Site["Sz"]) + Lambda * (tensor_prod(Site["I"], Site["Sx"]) + tensor_prod(Site["Sx"], Site["I"]));
         B1 = Site.Basis1().Basis();
         B2 = Site.Basis1().Basis();
         Obs = Site["Sx"];
      }
      else if (HamStr == "itf-z2")
      {
         std::cout << "Hamiltonian is transverse-field Ising with Z2, J=" << J << ", Lambda=" << Lambda << "\n";
         SiteBlock Site = CreateZ2SpinSite(0.5);
         H = J * 4.0 * tensor_prod(Site["Sz"], Site["Sz"]) + Lambda * (tensor_prod(Site["I"], Site["Sx"]) + tensor_prod(Site["Sx"], Site["I"]));
         B1 = Site.Basis1().Basis();
         B2 = Site.Basis1().Basis();
         Obs = Site["Sx"];
      }

      InfiniteWavefunction Psi;
      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::ImportHeap(FName);
      //      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      Psi = *PsiPtr;

      StateComponent Phi = Psi.Psi.get_front();
      MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
      MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
      Phi = prod(prod(LambdaInvSqrt, Phi), LambdaSqrt);

      OperatorComponent MARaw, MBRaw;
      std::tie(MARaw, MBRaw) = decompose_tensor_prod(Exponentiate(H * std::complex<double>(0.0, Timestep)), B1, B2);
      OperatorComponent MA = aux_tensor_prod(MARaw, MBRaw);
      OperatorComponent MB = aux_tensor_prod(MBRaw, MARaw);
      OperatorComponent MATransverse = exchange(MA);
      OperatorComponent MBTransverse = exchange(MB);

      OperatorComponent MAHermRaw, MBHermRaw;
      std::tie(MAHermRaw, MBHermRaw) = decompose_tensor_prod(Exponentiate(H * std::complex<double>(0.0, -Timestep)), B1, B2);
      OperatorComponent MAHerm = aux_tensor_prod(MAHermRaw, MBHermRaw);
      OperatorComponent MBHerm = aux_tensor_prod(MBHermRaw, MAHermRaw);
      OperatorComponent MAHermTransverse = exchange(MAHerm);
      OperatorComponent MBHermTransverse = exchange(MBHerm);


#if defined(SECOND_ORDER_ST)
      OperatorComponent MARawHalf, MBRawHalf;
      std::tie(MARawHalf, MBRawHalf) = decompose_tensor_prod(Exponentiate(H * std::complex<double>(0.0, 0.5*Timestep)), B1, B2);
      OperatorComponent MAHalf = aux_tensor_prod(MARawHalf, MBRawHalf);
      OperatorComponent MBHalf = aux_tensor_prod(MBRawHalf, MARawHalf);
      OperatorComponent MATransverseHalf = exchange(MAHalf);
      OperatorComponent MBTransverseHalf = exchange(MBHalf);

      OperatorComponent MAHermRawHalf, MBHermRawHalf;
      std::tie(MAHermRawHalf, MBHermRawHalf) = decompose_tensor_prod(Exponentiate(H * std::complex<double>(0.0, -0.5*Timestep)), B1, B2);
      OperatorComponent MAHermHalf = aux_tensor_prod(MAHermRawHalf, MBHermRawHalf);
      OperatorComponent MBHermHalf = aux_tensor_prod(MBHermRawHalf, MAHermRawHalf);
      OperatorComponent MAHermTransverseHalf = exchange(MAHermHalf);
      OperatorComponent MBHermTransverseHalf = exchange(MBHermHalf);
#endif

      // Assemble the operator
      std::vector<OperatorComponent> TransOp1, TransOp2;

      // Add the left boundary and the conjugate tensors (backwards time evolution)
      TransOp1.push_back(RotateToOperatorLeftBoundary(Phi));
      TransOp2.push_back(RotateToOperatorLeftBoundary(Phi));

#if defined(SECOND_ORDER_ST)
      TransOp1.push_back(MBHermTransverseHalf);
      TransOp2.push_back(MAHermTransverseHalf);
#else
         TransOp1.push_back(MBHermTransverse);
         TransOp2.push_back(MAHermTransverse);
#endif

      for (int i = 0; i < NumTimesteps-1; ++i)
      {
         TransOp1.push_back(MAHermTransverse);
         TransOp2.push_back(MBHermTransverse);

         TransOp1.push_back(MBHermTransverse);
         TransOp2.push_back(MAHermTransverse);
      }

      TransOp1.push_back(MAHermTransverse);
      TransOp2.push_back(MBHermTransverse);

#if defined(SECOND_ORDER_ST)
      TransOp1.push_back(MBHermTransverseHalf);
      TransOp2.push_back(MAHermTransverseHalf);
#endif

      // the observable
      std::vector<OperatorComponent> TransOp1Obs = TransOp1;
      TransOp1Obs.back() = prod(TransOp1Obs.back(), Obs);

#if defined(SECOND_ORDER_ST)
      TransOp1.push_back(MBTransverseHalf);
      TransOp1Obs.push_back(MBTransverseHalf);
      TransOp2.push_back(MATransverseHalf);
#else
      TransOp1.push_back(MBTransverse);
      TransOp1Obs.push_back(MBTransverse);
      TransOp2.push_back(MATransverse);
#endif

      TransOp1.push_back(MATransverse);
      TransOp1Obs.push_back(MATransverse);
      TransOp2.push_back(MBTransverse);

      // now the forwards time evolution section
      for (int i = 0; i < NumTimesteps-1; ++i)
      {
         TransOp1.push_back(MBTransverse);
         TransOp1Obs.push_back(MBTransverse);
         TransOp2.push_back(MATransverse);

         TransOp1.push_back(MATransverse);
         TransOp1Obs.push_back(MATransverse);
         TransOp2.push_back(MBTransverse);
      }

#if defined(SECOND_ORDER_ST)
      TransOp1.push_back(MBTransverseHalf);
      TransOp1Obs.push_back(MBTransverseHalf);
      TransOp2.push_back(MATransverseHalf);
#endif

      // the right boundary operators
      TransOp1.push_back(RotateToOperatorRightBoundary(Phi));
      TransOp1Obs.push_back(RotateToOperatorRightBoundary(Phi));
      TransOp2.push_back(RotateToOperatorRightBoundary(Phi));

      // make the MPOperators

      MPOperator TransOperator1 = MPOperator(TransOp1.begin(), TransOp1.end());
      MPOperator TransOperator1WithObs = MPOperator(TransOp1Obs.begin(), TransOp1Obs.end());
      MPOperator TransOperator2 = MPOperator(TransOp2.begin(), TransOp2.end());

      // make a lattice file to represent this lattice

      SiteBlock SBlock = MakeSiteBlockFromBasis(TransOperator1[0].LocalBasis1());
      Lattice MyLattice(SBlock);
      for (unsigned i = 1; i < TransOperator1.size(); ++i)
      {
         MyLattice = join(MyLattice, MakeSiteBlockFromBasis(TransOperator1[i].LocalBasis1()));
      }

      OperatorList MyOpList(MyLattice);
      MyOpList["Op1"] = ConstructFromMPOperator(TransOperator1);
      MyOpList["Op2"] = ConstructFromMPOperator(TransOperator2);
      MyOpList["Obs"] = ConstructFromMPOperator(TransOperator1WithObs);

      MPOperator OpFull(TransOperator1.size());
      MPOperator ObsFull(TransOperator1.size());
      for (int i = 0; i < TransOperator1.size(); ++i)
      {
         OpFull[i] = aux_tensor_prod(TransOperator2[i], TransOperator1[i]);
         ObsFull[i] = aux_tensor_prod(TransOperator2[i], TransOperator1WithObs[i]);
      }

      MyOpList["Op"] =  ConstructFromMPOperator(OpFull);
      MyOpList["Ob"] =  ConstructFromMPOperator(ObsFull);

      pvalue_ptr<OperatorList> MyOpListPtr = new OperatorList(MyOpList);

      pheap::ShutdownPersistent(MyOpListPtr);

      ProcControl::Shutdown();
      return 0;
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
