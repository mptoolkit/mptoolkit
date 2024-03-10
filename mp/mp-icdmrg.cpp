// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-icdmrg.cpp
//
// Copyright (C) 2015-2023 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// single-site algorithm for optimizing an MPS while keeping it
// in canonical form

#include "mpo/basic_triangular_mpo.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include "common/formatting.h"
#include "common/prog_options.h"
#include <iostream>
#include "common/environment.h"
#include "common/unique.h"
#include "common/statistics.h"
#include "common/prog_opt_accum.h"
#include "mp-algorithms/gmres.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mp-algorithms/stateslist.h"
#include "wavefunction/operator_actions.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

#if !defined(NDEBUG)
#include "mp-algorithms/triangular_mpo_solver.h"
#endif

#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "mp-algorithms/eigensolver.h"

namespace prog_opt = boost::program_options;


void Optimize()
{
   // left orthogonalize
   // right orthogonalize
   // find fixed point Hamiltonian
   // optimize wavefunction
   
}




int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true, false);
   try
   {
      int NumIter = 20;
      int MinIter = 4;
      int MinStates = 1;
      std::string States = "";
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;
      std::string FName;
      std::string HamStr;
      bool Force = false;
      bool TwoSite = true;
      bool OneSite = false;
      int WavefuncUnitCellSize = 0;
      double MixFactor = 0.02;
      double RandomMixFactor = 0.0;
      bool NoFixedPoint = false;
      bool NoOrthogonalize = false;
      bool Create = false;
      bool ExactDiag = false;
      double FidelityScale = 1.0;
      int Verbose = 0;
      bool Quiet = false;
      bool DoRandom = false; // true if we want to start an iteration from a random centre matrix
      std::string TargetState;
      std::vector<std::string> BoundaryState;
      double EvolveDelta = 0.0;
      double InitialFidelity = 1E-7;
      //double ArnoldiTol = 1E-14;
      double GMRESTol = 1E-13;    // tolerance for GMRES for the initial H matrix elements.
                                  // ** 2016-01-25: 1e-14 seems too small, failed to converge with final tol 1.5e-14, increasing to 2e-14
                                  // ** 2016-02-23: increasing again to 3e-14 on an example that fails to converge
                                  // ** 2016-03-16: increased to 1e-13.  It perhaps needs to scale with the unit cell size?
      double MaxTol = 4E-4;  // never use an eigensolver tolerance larger than this
      double MinTol = 1E-16; // lower bound for the eigensolver tolerance - seems we dont really need it
      double HMix = 0;  // Hamiltonian length-scale mixing factor
      double ShiftInvertEnergy = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian, of the form lattice:operator")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
         ("force,f", prog_opt::bool_switch(&Force), "Allow overwriting output files")
         //("two-site", prog_opt::bool_switch(&TwoSite), "Modify two sites at once (default)") // not currently implemented
         //("one-site", prog_opt::bool_switch(&OneSite), "Modify one site at a time")
#if defined(ENABLE_ONE_SITE_SCHEME)
         ("onesiteboundary", prog_opt::bool_switch(&UseOneSiteScheme), "Modify one site also at the boundary")
#endif
         ("states,m", prog_opt::value(&States),
          FormatDefault("number of states, or a StatesList", States).c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
          FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("mix-factor", prog_opt::value(&MixFactor),
          FormatDefault("Mixing coefficient for the density matrix", MixFactor).c_str())
         ("random-mix-factor", prog_opt::value(&RandomMixFactor),
          FormatDefault("Random mixing for the density matrix", RandomMixFactor).c_str())
         ("hmix", prog_opt::value(&HMix),
          FormatDefault("Hamiltonian mixing factor", HMix).c_str())
         ("evolve", prog_opt::value(&EvolveDelta),
          "Instead of Lanczos, do imaginary time evolution with this timestep")
         ("random,a", prog_opt::bool_switch(&Create),
          "Create a new wavefunction starting from a random state")
         ("unitcell,u", prog_opt::value(&WavefuncUnitCellSize),
          "Only if --create is specified, the size of the wavefunction unit cell")
         ("startrandom", prog_opt::bool_switch(&DoRandom),
          "Start the first iDMRG iteration from a random centre matrix")
         ("exactdiag,e", prog_opt::bool_switch(&ExactDiag),
          "Start from an effective exact diagonalization of the unit cell")
         ("target,q", prog_opt::value(&TargetState),
          "the target quantum number per unit cell")
         ("boundary", prog_opt::value(&BoundaryState),
          "use this boundary quantum number for initializing the unit cell "
          "(useful for integer spin chains, can be used multiple times)")
         ("create,b", prog_opt::bool_switch(&NoFixedPoint),
          "Construct a new wavefunction from a random state or single-cell diagonalization")
         ("no-orthogonalize", prog_opt::bool_switch(&NoOrthogonalize),
          "Don't orthogonalize the wavefunction before saving")
         ("maxiter", prog_opt::value<int>(&NumIter),
          FormatDefault("Maximum number of Lanczos iterations per step (Krylov subspace size)", NumIter).c_str())
         ("miniter", prog_opt::value<int>(&MinIter),
          FormatDefault("Minimum number of Lanczos iterations per step", MinIter).c_str())
         ("maxtol", prog_opt::value(&MaxTol),
          FormatDefault("Maximum tolerance of the eigensolver", MaxTol).c_str())
         ("fidelityscale", prog_opt::value(&FidelityScale),
          FormatDefault("The tolerance of the eigensolver is min(maxtol, fidelityscale * sqrt(fidelity))",
                        FidelityScale).c_str())
         ("initialfidelity", prog_opt::value(&InitialFidelity),
          FormatDefault("Initial value for the fidelity to set the eigensolver tolerance, for the first iteration",
                        InitialFidelity).c_str())
         ("seed", prog_opt::value<unsigned long>(), "random seed")
         ("gmrestol", prog_opt::value(&GMRESTol),
          FormatDefault("tolerance for GMRES linear solver for the initial H matrix elements", GMRESTol).c_str())
	 ("solver", prog_opt::value<std::string>(),
	  "Eigensolver to use.  Supported values are lanzcos [default], arnoldi, arnoldi-lowest, shift-invert")
	 ("shift-invert-energy", prog_opt::value(&ShiftInvertEnergy),
	  "For the shift-invert solver, the target energy")

         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity (can be used more than once)")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::positional_options_description p;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      std::cout << "Starting iDMRG.  Hamiltonian = " << HamStr << '\n';
      std::cout << "Wavefunction = " << FName << std::endl;

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX)
         : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      if (vm.count("one-site"))
         TwoSite = !OneSite;

      bool StartFromFixedPoint = !NoFixedPoint; // we've reversed the option

      if (!StartFromFixedPoint && !ExactDiag)
         Create = true;  // start from random state

      // The main MPWavefunction object.  We use this for initialization (if we are starting from
      // an existing wavefunction), and it will be the final wavefunction that we save to disk.
      MPWavefunction Wavefunction;

      // The parameters for the iDMRG that we need to initialize
      LinearWavefunction Psi;
      QuantumNumber QShift;

      RealDiagonalOperator R;
      MatrixOperator UR;

      // Initialize the filesystem

      if (ExactDiag || Create)
      {
         pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
      }
      else
      {
         long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
         pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize);
         Wavefunction = *PsiPtr;

         InfiniteWavefunctionLeft StartingWavefunction = Wavefunction.get<InfiniteWavefunctionLeft>();

         std::tie(Psi, R) = get_left_canonical(StartingWavefunction);
         UR = MatrixOperator::make_identity(R.Basis2());
         QShift = StartingWavefunction.qshift();
      }

      // Hamiltonian
      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO;

      // get the Hamiltonian from the attributes, if it wasn't supplied
      if (HamStr.empty())
      {
         if (Wavefunction.Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: no Hamiltonian specified, use -H option or set wavefunction attribute Hamiltonian.\n";
            return 1;
         }
         HamStr = Wavefunction.Attributes()["Hamiltonian"].as<std::string>();
      }
      else
         Wavefunction.Attributes()["Hamiltonian"] = HamStr;

      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);
      int const UnitCellSize = HamMPO.size(); //Lattice.GetUnitCell().size();
      if (WavefuncUnitCellSize == 0)
         WavefuncUnitCellSize = UnitCellSize;

      optimize(HamMPO);

      // load the wavefunction
      if (ExactDiag)
      {
         QShift = QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), TargetState);
         std::cout << "Target quantum number = " << QShift << '\n';

         std::vector<BasisList> BL = ExtractLocalBasis1(HamMPO.data());
         std::vector<BasisList> FullBL = BL;
         while (int(FullBL.size()) < WavefuncUnitCellSize)
            std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));
         if (WavefuncUnitCellSize != int(FullBL.size()))
         {
            std::cout << "mp-idmrg: fatal: the wavefunction unit cell must be a multiple of the Hamiltonian unit cell.\n";
            return 1;
         }
         std::cout << "Creating exact diagonalization basis.  Wvaefunction unit cell size = "
                   << WavefuncUnitCellSize << '\n';

         QuantumNumbers::QuantumNumberList BoundaryQ;
         if (BoundaryState.empty())
         {
            BoundaryQ.push_back(QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList()));
         }
         else
         {
            for (unsigned i = 0; i < BoundaryState.size(); ++i)
            {
               std::cout << "Adding boundary quantum number " << BoundaryState[i] << '\n';
               BoundaryQ.push_back(QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), BoundaryState[i]));
            }
         }

         QuantumNumbers::QuantumNumberList LeftBoundary;
         for (unsigned i = 0; i < BoundaryQ.size(); ++i)
         {
            LeftBoundary.push_back(transform_targets(QShift, BoundaryQ[i])[0]);
         }

         if (LeftBoundary.size() == 0)
         {
            std::cerr << "fatal: the target quntum number is incompatible with the boundary quantum number"
               " for this unit cell.\n";
            return 1;
         }

         VectorBasis B1(HamMPO.front().GetSymmetryList());
         for (unsigned i = 0; i < LeftBoundary.size(); ++i)
         {
            B1.push_back(LeftBoundary[i], 1);
         }
         VectorBasis B2(HamMPO.front().GetSymmetryList());
         for (unsigned i = 0; i < BoundaryQ.size(); ++i)
         {
            B2.push_back(BoundaryQ[i], 1);
         }

         Psi.push_back(ConstructFromLeftBasis(FullBL[0], B1));
         for (int i = 1; i < WavefuncUnitCellSize; ++i)
         {
            Psi.push_back(ConstructFromLeftBasis(FullBL[i], Psi.get_back().Basis2()));
         }

         UR = MakeRandomMatrixOperator(Psi.Basis2(), B2);

         // adjust for periodic basis
         StateComponent x = prod(Psi.get_back(), UR);
         std::tie(R, UR);
         MatrixOperator X = Multiply(TruncateBasis2(x)); // the Basis2 is already 1-dim.  This just orthogonalizes x
         MatrixOperator U;
         SingularValueDecomposition(X, U, R, UR);
         x = prod(x, U);
         Psi.set_back(x);
      }
      else if (Create)
      {
         QShift = QuantumNumbers::QuantumNumber(HamMPO[0].GetSymmetryList(), TargetState);
         std::cout << "Target quantum number = " << QShift << '\n';

         std::cout << "Creating wavefunction.  Wavefunction unit cell size = " << WavefuncUnitCellSize << '\n';
         if (WavefuncUnitCellSize % UnitCellSize != 0)
         {
            std::cout << "mp-idmrg: fatal: the wavefunction unit cell must be a multiple of the Hamiltonian unit cell.\n";
            return 1;
         }
         std::vector<BasisList> BL = ExtractLocalBasis1(HamMPO.data());
         std::vector<BasisList> FullBL = BL;
         while (int(FullBL.size()) < WavefuncUnitCellSize)
            std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));

         QuantumNumber LBoundary, RBoundary;
         if (BoundaryState.empty())
         {
            RBoundary = QuantumNumber(HamMPO[0].GetSymmetryList());
            LBoundary = QShift;
         }
         else
         {
            RBoundary = QuantumNumber(HamMPO[0].GetSymmetryList(), BoundaryState[0]);
            std::cout << "Right boundary quantum number is " << RBoundary << '\n';
            if (BoundaryState.size() > 1)
            {
               std::cout << "WARNING: ignoring addititional boundary quantum numbers in random wavefunction\n";
            }
            QuantumNumbers::QuantumNumberList QL = transform_targets(QShift, RBoundary);
            if (QL.size() > 1)
            {
               PANIC("Don't know how to handle non-scalar non-abelian target state")(RBoundary)(QShift);
            }
            LBoundary = QL[0];
            std::cout << "Left boundary quantum number is " << LBoundary << '\n';
         }
         Psi = CreateRandomWavefunction(FullBL, LBoundary, 3, RBoundary);
         MatrixOperator X = left_orthogonalize(MatrixOperator::make_identity(Psi.Basis1()), Psi);
         MatrixOperator U;
         SingularValueDecomposition(X, U, R, UR);
         Psi.set_back(prod(Psi.get_back(), U));
      }

      WavefuncUnitCellSize = Psi.size();
      std::cout << "Wavefunction unit cell size = " << WavefuncUnitCellSize << '\n';
      if (WavefuncUnitCellSize % HamMPO.size() != 0)
      {
         std::cout << "mp-idmrg: fatal: the wavefunction unit cell must be a multiple of the Hamiltonian unit cell.\n";
         return 1;
      }

      if (vm.count("evolve"))
         std::cout << "Evolving with timestep " << EvolveDelta << '\n';

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = 0;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      StatesList MyStates(States);
      std::cout << MyStates << '\n';

      std::complex<double> InitialEnergy = 0.0;

      // replicate the HamMPO until it has the same size as the unit cell
      HamMPO = repeat(HamMPO, WavefuncUnitCellSize / HamMPO.size());
      CHECK_EQUAL(int(HamMPO.size()), WavefuncUnitCellSize);

      // Check that the local basis for the wavefunction and hamiltonian are compatible
      if (ExtractLocalBasis(Psi) != ExtractLocalBasis1(HamMPO))
      {
         std::cerr << "fatal: Hamiltonian is defined on a different local basis to the wavefunction.\n";
         return 1;
      }

      if (ExtractLocalBasis1(HamMPO) != ExtractLocalBasis2(HamMPO))
      {
         std::cerr << "fatal: Hamiltonian has different domain and co-domain.\n";
         return 1;
      }


      // Get the initial Hamiltonian matrix elements
      StateComponent BlockHamL = Initial_E(HamMPO, Psi.Basis1());
      if (StartFromFixedPoint)
      {
         std::cout << "Solving fixed-point Hamiltonian..." << std::endl;
         MatrixOperator Rho = delta_shift(scalar_prod(R, herm(R)), QShift);
         InitialEnergy = SolveHamiltonianMPO_Left(BlockHamL, Psi, QShift, HamMPO,
                                             Rho, GMRESTol, Verbose);
         std::cout << "Starting energy (left eigenvalue) = " << InitialEnergy << std::endl;

         //BlockHamL = delta_shift(BlockHamL, QShift);
      }

      StateComponent BlockHamR = Initial_F(HamMPO, Psi.Basis2());
      if (StartFromFixedPoint)
      {
         LinearWavefunction PsiR;
         RealDiagonalOperator D;
         std::tie(D, PsiR) = get_right_canonical(Wavefunction.get<InfiniteWavefunctionLeft>());

         BlockHamR = Initial_F(HamMPO, PsiR.Basis2());

         // check that we are orthogonalized
#if !defined(NDEBUG)
         MatrixOperator X = MatrixOperator::make_identity(PsiR.Basis2());
         X = inject_right(X, PsiR);
         CHECK(norm_frob(X - MatrixOperator::make_identity(PsiR.Basis1())) < 1E-12)(X);
#endif

         MatrixOperator Rho = D*D;

#if !defined(NDEBUG)
         MatrixOperator XX = Rho;
         XX = inject_left(XX, PsiR);
         CHECK(norm_frob(delta_shift(XX,QShift) - Rho) < 1E-12)(norm_frob(delta_shift(XX,QShift) - Rho) )(XX)(Rho);
#endif

         // We obtained Rho from the left side, so we need to delta shift to the right basis
         Rho = delta_shift(Rho, adjoint(QShift));

         std::complex<double> Energy = SolveHamiltonianMPO_Right(BlockHamR, PsiR, QShift, HamMPO,
                                                                 Rho, GMRESTol, Verbose);
         std::cout << "Starting energy (right eigenvalue) = " << Energy << std::endl;

      }

      // initialization complete.

      UR = delta_shift(UR, QShift);

      // Construct the iDMRG object
      iDMRG idmrg(Psi, R, UR, QShift, HamMPO, BlockHamL,
                  BlockHamR, InitialEnergy, Verbose);

      idmrg.MixingInfo.MixFactor = MixFactor;
      idmrg.MixingInfo.RandomMixFactor = RandomMixFactor;
      idmrg.SetInitialFidelity(InitialFidelity);
      idmrg.Solver().MaxTol = MaxTol;
      idmrg.Solver().MinTol = MinTol;
      idmrg.Solver().MinIter = MinIter;
      idmrg.Solver().MaxIter = NumIter;
      idmrg.Solver().FidelityScale = FidelityScale;
      idmrg.Solver().Verbose = Verbose;
      idmrg.Solver().EvolveDelta = EvolveDelta;
      if (vm.count("solver"))
      {
	 idmrg.Solver().SetSolver(vm["solver"].as<std::string>());
      }
      idmrg.Solver().SetShiftInvertEnergy(ShiftInvertEnergy);

      int ReturnCode = 0;

      try
      {
         bool First = true;
         for (int i = 0; i < MyStates.size(); ++i)
         {
            SInfo.MaxStates = MyStates[i].NumStates;

            if (i % 2 == 0)
            {
               idmrg.SweepLeft(SInfo, HMix, First);
               First = false;
            }
            else
            {
               idmrg.SweepRight(SInfo, HMix);
            }
         }
         idmrg.Finish(SInfo);

      }
      catch (ProcControl::Checkpoint& c)
      {
         ReturnCode = c.ReturnCode();
         std::cerr << "Early termination: "
                   << c.Reason() << '\n';
         EarlyTermination = true;
      }
      catch (...)
      {
         throw;      // if we got some other exception, don't even try and recover
      }

      // finished the iterations.
      std::cout << "Orthogonalizing wavefunction...\n";
      Wavefunction.Wavefunction() = InfiniteWavefunctionLeft::Construct(idmrg.Wavefunction(), idmrg.QShift, 0, Verbose+1);

      std::cerr << "Orthogonalization finished.\n";

      // any other attributes?
      Wavefunction.Attributes()["LastEnergy"] = formatting::format_complex(idmrg.Solver().LastEnergy());
      Wavefunction.SetDefaultAttributes();

      // History log
      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));

      // save wavefunction
      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Wavefunction));
      pheap::ShutdownPersistent(P);

      ProcControl::Shutdown();
      return ReturnCode;
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      if (e.Why == "File exists")
         std::cerr << "Note: use --force (-f) option to overwrite.\n";
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
