// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-dmrg.cpp
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
//
// IBC DMRG
//
// Experimental DMRG code for infinite boundary conditions

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
#include "common/prog_options.h"
#include <iostream>
#include "common/environment.h"
#include "common/statistics.h"
#include "common/prog_opt_accum.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mp-algorithms/stateslist.h"
#include "wavefunction/operator_actions.h"
#include "wavefunction/ibc.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

#if !defined(NDEBUG)
#include "mp-algorithms/triangular_mpo_solver.h"
#endif

#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "lattice/infinite-parser.h"

#include "mp-algorithms/eigensolver.h"
#include "mp-algorithms/triangular_mpo_solver.h"

namespace prog_opt = boost::program_options;

using statistics::moving_average;
using statistics::moving_exponential;

struct MixInfo
{
   double MixFactor;
   double RandomMixFactor;
};

bool EarlyTermination = false;  // we set this to true if we get a checkpoint


//#define SSC

// Apply subspace expansion / truncation on the left (C.Basis1()).
// Returns a matrix Lambda (diagonal) and a unitary
// Postcondition: U' Lambda' C' = C (up to truncation!)
std::pair<MatrixOperator, RealDiagonalOperator>
SubspaceExpandBasis1(StateComponent& C, OperatorComponent const& H, StateComponent const& RightHam,
                     MixInfo const& Mix, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& LeftHam)
{
   // truncate - FIXME: this is the s3e step
#if defined(SSC)
   MatrixOperator Lambda;
   SimpleStateComponent CX;
   std::tie(Lambda, CX) = ExpandBasis1_(C);
#else
   MatrixOperator Lambda = ExpandBasis1(C);
#endif

   MatrixOperator Rho = scalar_prod(herm(Lambda), Lambda);
   if (Mix.MixFactor > 0)
   {
#if defined(SSC)
      StateComponent RH = contract_from_right(herm(H), CX, RightHam, herm(CX));
#else
      StateComponent RH = contract_from_right(herm(H), C, RightHam, herm(C));
#endif
      MatrixOperator RhoMix;
      MatrixOperator RhoL = scalar_prod(Lambda, herm(Lambda));

      // Skip the identity and the Hamiltonian
      for (unsigned i = 1; i < RH.size()-1; ++i)
      {
         double Prefactor = trace(triple_prod(herm(LeftHam[i]), RhoL, LeftHam[i])).real();
         if (Prefactor == 0)
            Prefactor = 1;
         RhoMix += Prefactor * triple_prod(herm(RH[i]), Rho, RH[i]);
      }
      //      MatrixOperator RhoMix = operator_prod(herm(RH), Rho, RH);
      Rho += (Mix.MixFactor / trace(RhoMix)) * RhoMix;
   }
   if (Mix.RandomMixFactor > 0)
   {
      MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
      RhoMix = herm(RhoMix) * RhoMix;
      Rho += (Mix.RandomMixFactor / trace(RhoMix)) * RhoMix;
   }

   //TRACE(Rho);

   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);
   Lambda = Lambda * herm(U);

   //TRACE(Lambda);

#if defined(SSC)
   C = U*CX; //prod(U, CX);
#else
   C = prod(U, C);
#endif

   MatrixOperator Vh;
   RealDiagonalOperator D;
   SingularValueDecompositionKeepBasis2(Lambda, U, D, Vh);

   //TRACE(U)(D)(Vh);

   C = prod(Vh, C);
   return std::make_pair(U, D);
}

// Apply subspace expansion / truncation on the right (C.Basis2()).
// Returns Lambda matrix (diagonal) and a unitary matrix
// Postcondition: C' Lambda' U' = C (up to truncation!)
std::pair<RealDiagonalOperator, MatrixOperator>
SubspaceExpandBasis2(StateComponent& C, OperatorComponent const& H, StateComponent const& LeftHam,
                     MixInfo const& Mix, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& RightHam)
{
   // truncate - FIXME: this is the s3e step
   MatrixOperator Lambda = ExpandBasis2(C);

   MatrixOperator Rho = scalar_prod(Lambda, herm(Lambda));
   if (Mix.MixFactor > 0)
   {
      StateComponent LH = contract_from_left(H, herm(C), LeftHam, C);
      MatrixOperator RhoMix;

      MatrixOperator RhoR = scalar_prod(herm(Lambda), Lambda);

      for (unsigned i = 1; i < LH.size()-1; ++i)
      {
         double Prefactor = trace(triple_prod(herm(RightHam[i]), RhoR, RightHam[i])).real();
         if (Prefactor == 0)
            Prefactor = 1;
         RhoMix += Prefactor * triple_prod(LH[i], Rho, herm(LH[i]));
      }
      Rho += (Mix.MixFactor / trace(RhoMix)) * RhoMix;
   }
   if (Mix.RandomMixFactor > 0)
   {
      MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
      RhoMix = herm(RhoMix) * RhoMix;
      Rho += (Mix.RandomMixFactor / trace(RhoMix)) * RhoMix;
   }
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);

   Lambda = U * Lambda;
   C = prod(C, herm(U));

   MatrixOperator Vh;
   RealDiagonalOperator D;
   SingularValueDecompositionKeepBasis1(Lambda, U, D, Vh);

   C = prod(C, U);

   return std::make_pair(D, Vh);
}





class IBC_DMRG
{
   public:
      // Construct an iDMRG object.  It is assumed that Psi_ is in left-canonical form, with
      // LambdaR being the lambda matrix on the right edge.
      IBC_DMRG(LinearWavefunction const& Psi_, MatrixOperator const& LambdaR,
               BasicTriangularMPO const& Hamiltonian_,
               StateComponent const& LeftHam, StateComponent const& RightHam,
               int Verbose = 0);

      void SetMixInfo(MixInfo const& m);

      void SetInitialFidelity(double f)
      {
         Solver_.SetInitialFidelity(Psi.size(), f);
      }

      LocalEigensolver& Solver() { return Solver_; }

      void SweepRight(StatesInfo const& SInfo, double HMix = 0);
      void SweepLeft(StatesInfo const& SInfo, double HMix = 0);

      // returns the wavefunction, which is in center-regular form
      LinearWavefunction const& Wavefunction() const { return Psi; }

      // Above functions are implemented in terms of:

      void ConstructInitialHamiltonian();

      void Solve();

      void TruncateAndShiftLeft(StatesInfo const& States);
      void TruncateAndShiftRight(StatesInfo const& States);

      void ShowInfo(char c);

      void CheckConsistency() const;


      //   private:
      BasicTriangularMPO Hamiltonian;
      LinearWavefunction Psi;

      std::deque<StateComponent> LeftHamiltonian;
      std::deque<StateComponent> RightHamiltonian;

      LinearWavefunction::iterator C;
      BasicTriangularMPO::const_iterator H;

      int Verbose;

      // iterators pointing to the edges of the unit cell.
      // FirstSite = Psi.begin()
      // LastSite = Psi.end() - 1
      LinearWavefunction::const_iterator FirstSite, LastSite;

      LocalEigensolver Solver_;

      MixInfo    MixingInfo;
      TruncationInfo Info;

      int SweepNumber;
};



IBC_DMRG::IBC_DMRG(LinearWavefunction const& Psi_, MatrixOperator const& LambdaR,
                   BasicTriangularMPO const& Hamiltonian_,
                   StateComponent const& LeftHam, StateComponent const& RightHam,
                   int Verbose_)
   : Hamiltonian(Hamiltonian_),
     Psi(Psi_),
     Verbose(Verbose_),
     SweepNumber(1)
{
   FirstSite = Psi.begin();
   LastSite = Psi.end();
   --LastSite;

   C = Psi.end();
   --C;
   H = Hamiltonian.end();
   --H;

   *C = prod(*C, LambdaR);

   LeftHamiltonian.push_back(LeftHam);
   RightHamiltonian.push_front(RightHam);

   this->ConstructInitialHamiltonian();
}

void
IBC_DMRG::ConstructInitialHamiltonian()
{
   LinearWavefunction::const_iterator c = Psi.begin();
   BasicTriangularMPO::const_iterator h = Hamiltonian.begin();

   while (c != C)
   {
      LeftHamiltonian.push_back(contract_from_left(*h, herm(*c), LeftHamiltonian.back(), *c));
      ++c;
      ++h;
   }
   CHECK(h == H);
}

void
IBC_DMRG::TruncateAndShiftLeft(StatesInfo const& States)
{
   this->CheckConsistency();
   // Truncate right
   MatrixOperator U;
   RealDiagonalOperator Lambda;
   std::tie(U, Lambda) = SubspaceExpandBasis1(*C, *H, RightHamiltonian.front(), MixingInfo, States, Info,
                                                LeftHamiltonian.back());

   if (Verbose > 1)
   {
      std::cerr << "Truncating left basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks
   RightHamiltonian.push_front(contract_from_right(herm(*H), *C, RightHamiltonian.front(), herm(*C)));
   LeftHamiltonian.pop_back();

   // next site
   --H;
   --C;

   *C = prod(*C, U*Lambda);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   this->CheckConsistency();
}

void
IBC_DMRG::CheckConsistency() const
{
   CHECK_EQUAL(LeftHamiltonian.back().Basis2(), C->Basis1());
   CHECK_EQUAL(RightHamiltonian.front().Basis1(), C->Basis2());
   CHECK_EQUAL(LeftHamiltonian.back().LocalBasis(), H->Basis1());
   CHECK_EQUAL(RightHamiltonian.front().LocalBasis(), H->Basis2());
   CHECK_EQUAL(H->LocalBasis2(), C->LocalBasis());
}

void
IBC_DMRG::TruncateAndShiftRight(StatesInfo const& States)
{
   // Truncate right
   RealDiagonalOperator Lambda;
   MatrixOperator U;
   std::tie(Lambda, U) = SubspaceExpandBasis2(*C, *H, LeftHamiltonian.back(), MixingInfo, States, Info,
                                                RightHamiltonian.front());
   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks
   LeftHamiltonian.push_back(contract_from_left(*H, herm(*C), LeftHamiltonian.back(), *C));
   RightHamiltonian.pop_front();

   // next site
   ++H;
   ++C;

   *C = prod(Lambda*U, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   this->CheckConsistency();
}

void
IBC_DMRG::SweepRight(StatesInfo const& States, double HMix)
{
   this->Solve();
   this->ShowInfo('P');

   while (C != LastSite)
   {
      this->TruncateAndShiftRight(States);
      this->Solve();
      this->ShowInfo('R');
   }

   SweepNumber++;
}

void
IBC_DMRG::SweepLeft(StatesInfo const& States, double HMix)
{
   this->Solve();
   this->ShowInfo('Q');

   while (C != FirstSite)
   {
      this->TruncateAndShiftLeft(States);
      this->Solve();
      this->ShowInfo('L');
   }

   SweepNumber++;
}

void
IBC_DMRG::Solve()
{
   Solver_.Solve(*C, LeftHamiltonian.back(), *H, RightHamiltonian.front());
}

void
IBC_DMRG::ShowInfo(char c)
{
   std::cout << c
             << " Sweep=" << SweepNumber
             << " Energy=" << Solver_.LastEnergy()
             << " States=" << Info.KeptStates()
             << " TruncError=" << Info.TruncationError()
             << " Entropy=" << Info.KeptEntropy()
             << " Fidelity=" << Solver_.LastFidelity()
      //                         << " FidelityAv=" << FidelityAv.value()
             << " Iter=" << Solver_.LastIter()
             << " Tol=" << Solver_.LastTol()
             << '\n';
}







int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true);
   try
   {
      // flush cout if we write to cerr
      std::cout.tie(&std::cerr);

      int NumIter = 20;
      int MinIter = 4;
      int MinStates = 1;
      std::string States = "100";
      int NumSteps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;
      std::string FName;
      std::string HamStr;
      std::string CouplingFile;
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
      double ArnoldiTol = 1E-14;
      double GMRESTol = 1E-13;    // tolerance for GMRES for the initial H matrix elements.
                                  // ** 2016-01-25: 1e-14 seems too small, failed to converge with final tol 1.5e-14, increasing to 2e-14
                                  // ** 2016-02-23: increasing again to 3e-14 on an example that fails to converge
                                  // ** 2016-03-16: increased to 1e-13.  It perhaps needs to scale with the unit cell size?
      double MaxTol = 4E-4;  // never use an eigensolver tolerance larger than this
      double MinTol = 1E-16; // lower bound for the eigensolver tolerance - seems we dont really need it
      double HMix = 0;  // Hamiltonian length-scale mixing factor

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian, of the form lattice:operator")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
         ("force,f", prog_opt::bool_switch(&Force), "Allow overwriting output files")
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

      std::cout << "Starting DMRG.  Hamiltonian = " << HamStr << '\n';
      std::cout << "Wavefunction = " << FName << std::endl;

      // The main MPWavefunction object.
      MPWavefunction Wavefunction;

      // The parameters for the IBC_DMRG that we need to initialize
      IBCWavefunction Psi;

      // Initialize the filesystem

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize);
      Wavefunction = *PsiPtr;
      Psi = Wavefunction.get<IBCWavefunction>();

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
      int const UnitCellSize = Lattice.GetUnitCell().size();
      if (WavefuncUnitCellSize == 0)
         WavefuncUnitCellSize = UnitCellSize;

      optimize(HamMPO);



      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = 0;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      StatesList MyStates(States);
      if (vm.count("steps") && MyStates.size() == 1)
      {
         MyStates.Repeat(NumSteps);
      }
      std::cout << MyStates << '\n';

      BasicTriangularMPO LeftHamMPO = HamMPO;
      BasicTriangularMPO RightHamMPO = HamMPO;

      // replicate the HamMPO until it has the same size as the unit cell
      HamMPO = repeat(HamMPO, Psi.window_size() / HamMPO.size());
      CHECK_EQUAL(int(HamMPO.size()), Psi.window_size());

      // Check that the local basis for the wavefunction and hamiltonian are compatible
      // Check that the local basis for the wavefunction and hamiltonian are compatible
      if (ExtractLocalBasis(Psi.Window) != ExtractLocalBasis1(HamMPO))
      {
         std::cerr << "fatal: Hamiltonian is defined on a different local basis to the wavefunction.\n";
         return 1;
      }

      if (ExtractLocalBasis1(HamMPO) != ExtractLocalBasis2(HamMPO))
      {
         std::cerr << "fatal: Hamiltonian has different domain and co-domain.\n";
         return 1;
      }

      // Get the fixed point Hamiltonian matrix elements
      std::cout << "Solving fixed-point Hamiltonian..." << std::endl;
      StateComponent BlockHamL = Initial_E(HamMPO, Psi.Left.Basis2());
      std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, Psi.Left, LeftHamMPO,
                                                            GMRESTol, Verbose);
      std::cout << "Starting energy (left eigenvalue) = " << LeftEnergy << std::endl;

      StateComponent BlockHamR = Initial_F(HamMPO, Psi.Right.Basis1());
      std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, Psi.Right, RightHamMPO,
                                                              GMRESTol, Verbose);
      std::cout << "Starting energy (right eigenvalue) = " << RightEnergy << std::endl;

      // The LinearWavefunction representation of the Window
      LinearWavefunction PsiLinear;
      MatrixOperator Lambda;
      std::tie(PsiLinear, Lambda) = get_left_canonical(Psi.Window);

      // Construct the IBC_DMRG object
      IBC_DMRG dmrg(PsiLinear, Lambda, HamMPO, BlockHamL, BlockHamR, Verbose);

      dmrg.MixingInfo.MixFactor = MixFactor;
      dmrg.MixingInfo.RandomMixFactor = RandomMixFactor;
      dmrg.SetInitialFidelity(InitialFidelity);
      dmrg.Solver().MaxTol = MaxTol;
      dmrg.Solver().MinTol = MinTol;
      dmrg.Solver().MinIter = MinIter;
      dmrg.Solver().MaxIter = NumIter;
      dmrg.Solver().FidelityScale = FidelityScale;
      dmrg.Solver().Verbose = Verbose;
      dmrg.Solver().EvolveDelta = EvolveDelta;

      // Do the DMRG

      int ReturnCode = 0;

      try
      {
         for (int i = 0; i < MyStates.size(); ++i)
         {
            SInfo.MaxStates = MyStates[i].NumStates;

            if (i % 2 == 0)
            {
               dmrg.SweepLeft(SInfo, HMix);
            }
            else
            {
               dmrg.SweepRight(SInfo, HMix);
            }
         }
         //idmrg.Finish(SInfo);

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

      // Save the wavefunction
      std::cerr << "Orthogonalizing wavefunction...\n";

      PsiLinear = dmrg.Wavefunction();

      // We need to pull out the final Lambda matrix
      RealDiagonalOperator D;
      MatrixOperator M;
      StateComponent A = PsiLinear.get_back();
      std::tie(D, M) = OrthogonalizeBasis2(A);
      PsiLinear.set_back(A);
      M = D*M;

      Psi.Window = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiLinear), M, Verbose);

      // any other attributes?
      Wavefunction.Attributes()["LastEnergy"] = dmrg.Solver().LastEnergy();
      Wavefunction.SetDefaultAttributes();

      // History log
      Wavefunction.AppendHistory(EscapeCommandline(argc, argv));

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
