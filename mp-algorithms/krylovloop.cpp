// -*- C++ -*- $Id$

#include "krylovloop.h"
#include "linearalgebra/eigen.h"

double KrylovLoop::TimestepFromHochbruckBound(int NumKrylov, double HNorm, double EBound)
{
   double p2 = EBound * std::pow((NumKrylov+1)/2.72, NumKrylov+1) * std::sqrt(math_const::pi*2*(NumKrylov+1));
   double tn = (1.0 / HNorm) * std::pow(p2 / (4.0 * HNorm), 1.0 / NumKrylov);
   double p1 = pow(10.0, int(std::floor(std::log10(tn) - std::sqrt(0.1) + 0.5)-1));
   return std::floor(tn/p1 + 1.05) * p1;
}

double KrylovLoop::SweepRight(double MinTrunc, bool FullMixing)
{
   Solver_.ExpandLeft(FullOrtho_);
   if (TwoSite_) Solver_.ExpandRight(FullOrtho_);
   Solver_.Solve(FullOrtho_);
   double TruncPerStep = MinTrunc / Solver_.Krylov[0].size();
   double Trunc = 0;
   TruncationInfo Info = Solver_.TruncateLeft(MinStates_, MaxStates_, 
                                              TruncPerStep, MixFactor_, FullMixing, FullOrtho_);
   Trunc += Info.TruncationError();
   std::cout << "partition=(" << Solver_.LeftSize() << ',' << Solver_.RightSize() << ")"
             << " states=" << Info.KeptStates() << " trunc=" << Info.TruncationError()
             << " cumulative-trunc=" << Trunc << "\n";
   // sweep right
   while (Solver_.RightSize() > 1)
   {
      Solver_.ShiftRightAndExpand(FullOrtho_);
      if (TwoSite_) Solver_.ExpandRight(FullOrtho_);
      Solver_.Solve(FullOrtho_);
      Info = Solver_.TruncateLeft(MinStates_, MaxStates_, TruncPerStep, MixFactor_, FullMixing, 
                                  FullOrtho_);
      Trunc += Info.TruncationError();
      std::cout << "partition=(" << Solver_.LeftSize() << ',' << Solver_.RightSize() << ")"
                << " states=" << Info.KeptStates() << " trunc=" << Info.TruncationError()
                << " cumulative-trunc=" << Trunc << "\n";
   }
   return Trunc;
}

double KrylovLoop::SweepLeft(double MinTrunc, bool FullMixing)
{
   Solver_.ExpandRight(FullOrtho_);
   if (TwoSite_) Solver_.ExpandLeft(FullOrtho_);
   Solver_.Solve(FullOrtho_);
   double TruncPerStep = MinTrunc / Solver_.Krylov[0].size();
   double Trunc = 0;
   TruncationInfo Info = Solver_.TruncateRight(MinStates_, MaxStates_, TruncPerStep, 
                                               MixFactor_, FullMixing, FullOrtho_);
   Trunc += Info.TruncationError();
   std::cout << "partition=(" << Solver_.LeftSize() << ',' << Solver_.RightSize() << ")"
             << " states=" << Info.KeptStates() << " trunc=" << Info.TruncationError()
             << " cumulative-trunc=" << Trunc << "\n";
   // sweep left
   while (Solver_.LeftSize() > 1)
   {
      Solver_.ShiftLeftAndExpand(FullOrtho_);
      if (TwoSite_) Solver_.ExpandLeft(FullOrtho_);
      Solver_.Solve(FullOrtho_);
      Info = Solver_.TruncateRight(MinStates_, MaxStates_, TruncPerStep, 
                                   MixFactor_, FullMixing, FullOrtho_);
      Trunc += Info.TruncationError();
      std::cout << "partition=(" << Solver_.LeftSize() << ',' << Solver_.RightSize() << ")"
                << " states=" << Info.KeptStates() << " trunc=" << Info.TruncationError()
                << " cumulative-trunc=" << Trunc << "\n";
   }
   return Trunc;
}

void KrylovLoop::AddKrylovVector(WavefunctionType const& K)
{
   Solver_.AddKrylovVector(K);
}

void KrylovLoop::AddKrylovVector()
{
   Solver_.AddKrylovVector(Solver_.Krylov.back());
}

void KrylovLoop::EvolveNextKrylov(double RequiredVariance)
{
   int NumSweeps = 0;
   int const SweepsBailout = 2;
   double Trunc = RequiredVariance;
   double Variance = RequiredVariance + 1;  // satisfy the initial loop invariant
   while (Variance > RequiredVariance && NumSweeps < MaxSweeps_)
   {
      if (NumSweeps > 1 && NumSweeps % SweepsBailout == 0)
      {
         Trunc *= 0.25;
         std::cout << "Not converged after " << NumSweeps
                   << " sweeps, decreasing DMRG truncation to " << Trunc << '\n';
      }
      ++NumSweeps;

      this->SweepRight(Trunc, true);
      this->SweepLeft(Trunc, true);
      Variance = Solver_.Variance(FullOrtho_);
      TRACE(Variance);
   }
   std::cout << (Variance <= RequiredVariance ? "Krylov vector has converged." 
                 : "Krylov vector has NOT converged, but max-sweeps has been reached.");
   std::cout << "  Required variance = " << RequiredVariance
                << ", actual variance = " << Variance << ", number of sweeps = " << NumSweeps << std::endl;
   if (DoReductionSweep_)
   {
      std::cout << "Doing reduction sweep.\n";
      Variance = RequiredVariance + 1;  // satisfy the initial loop invariant
      while (Variance > RequiredVariance)
      {
         if (NumSweeps > 1 && NumSweeps % SweepsBailout == 0)
         {
            Trunc *= 0.25;
            std::cout << "Not converged after " << NumSweeps
                      << " sweeps, decreasing DMRG truncation to " << Trunc << '\n';
         }
         ++NumSweeps;

         this->SweepRight(Trunc, false);
         this->SweepLeft(Trunc, false);
         Variance = Solver_.Variance(FullOrtho_);
         TRACE(Variance)(RequiredVariance);
      }
   }
   Solver_.FixKrylovVector(FullOrtho_);
}

void KrylovLoop::ConstructKrylovBasis(int NumKrylov, double GuessTimestep)
{
   for (int i = 0; i < NumKrylov; ++i)
   {
      //      double MatElement = Solver_.SubMatrixElement();
      //      Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(GuessTimestep*TimeDirection_, MatElement);

      // modified
      Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(GuessTimestep*TimeDirection_);

      double Coeff = LinearAlgebra::norm_2(CoeffVec[CoeffVec.size()-1]);
      double NextError = GuessTimestep * ErrorBound_ / Coeff;
      double NextVariance = NextError * NextError / NumKrylov;
      TRACE(NumKrylov)(GuessTimestep)(CoeffVec)(NextError)(NextVariance)(Coeff);
      // TRACE(MatElement);
      this->AddKrylovVector();
      this->EvolveNextKrylov(NextVariance);
   }
   TRACE(this->SubspaceIdentity())(EigenvaluesHermitian(this->SubspaceIdentity()));
}

void KrylovLoop::ConstructKrylovBasis(std::complex<double> Timestep)
{
   //double MatElement = Solver_.SubMatrixElement();
   //TRACE(MatElement);
   //Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(Timestep, MatElement);

   // modified
   Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(Timestep);
   TRACE(CoeffVec);
   double Coeff = LinearAlgebra::norm_2_sq(CoeffVec[CoeffVec.size()-1]);

   LastNumKrylov_=1;

   double NextError = ErrorScaleFactor_ * ErrorBound_ / (Coeff * LastNumKrylov_);

   while (NextError < ErrorScaleFactor_ / KrylovCutoffFactor_)
   {
      TRACE(CoeffVec)(NextError)(Coeff);
      this->AddKrylovVector();
      this->EvolveNextKrylov(NextError);

      // MatElement = Solver_.SubMatrixElement();
      // CoeffVec = Solver_.Exponentiate(Timestep, MatElement);

      // modified
      CoeffVec = Solver_.Exponentiate(Timestep);

      Coeff = LinearAlgebra::norm_2_sq(CoeffVec[CoeffVec.size()-1]);
      
      if (OldCoeffVec_.size() >= CoeffVec.size()) 
         if (Coeff < LinearAlgebra::norm_2_sq(OldCoeffVec_[CoeffVec.size()-1]))
            Coeff = LinearAlgebra::norm_2_sq(OldCoeffVec_[CoeffVec.size()-1]);

      
      if (LastNumKrylov_ < int(Solver_.Krylov.size()+1))
         LastNumKrylov_ = Solver_.Krylov.size()+1;
      NextError = ErrorScaleFactor_ * ErrorBound_ / (Coeff * LastNumKrylov_);
   }
   LastNumKrylov_ = Solver_.Krylov.size();
   
   OldCoeffVec_ = CoeffVec;

   TRACE(CoeffVec)(NextError)(Coeff)(LastNumKrylov_);
   TRACE(EigenvaluesHermitian(this->SubspaceIdentity()));
}

double KrylovLoop::GetTimestep(double GuessTimestep) const
{
   double MinTimestep = GuessTimestep;
   double MaxTimestep = GuessTimestep;
   //double MatElement = Solver_.SubMatrixElement();
   //Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(GuessTimestep*TimeDirection_, MatElement);

   // modified
   Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(GuessTimestep*TimeDirection_);

   double Coeff = LinearAlgebra::norm_2(CoeffVec[CoeffVec.size()-1]);

   // If our error is too big, reduce the timestep until it is in range, leave MaxTimestep as the
   // last out-of-range time
   if (Coeff > GuessTimestep * ErrorBound_)
   {
      while (Coeff > GuessTimestep * ErrorBound_)
      {
         MaxTimestep = GuessTimestep;
         GuessTimestep *= 0.5;
         MinTimestep = GuessTimestep;
         // modified
         //CoeffVec = Solver_.Exponentiate(GuessTimestep*TimeDirection_, MatElement);
         CoeffVec = Solver_.Exponentiate(GuessTimestep*TimeDirection_);
         Coeff = LinearAlgebra::norm_2(CoeffVec[CoeffVec.size()-1]);
      }
   }
   else
   {
      // if our error is too small, increase the time step
      double TryCoeff = Coeff;
      double TryTimestep = GuessTimestep;
      while (TryCoeff < TryTimestep * ErrorBound_)
      {
         Coeff = TryCoeff;
         GuessTimestep = TryTimestep;

         MinTimestep = TryTimestep;
         TryTimestep *= 2;
         MaxTimestep = TryTimestep;
         //CoeffVec = Solver_.Exponentiate(TryTimestep*TimeDirection_, MatElement);
         CoeffVec = Solver_.Exponentiate(TryTimestep*TimeDirection_);
         TryCoeff = LinearAlgebra::norm_2(CoeffVec[CoeffVec.size()-1]);
      }
   }

   // Now we have fixed bounds on the possible timesteps.  Do a binary search.
   TRACE(MinTimestep)(MaxTimestep);

   while (Coeff < 0.9 * (GuessTimestep * ErrorBound_ ))
   {
      double TryTimestep = 0.5 * (MinTimestep + MaxTimestep);
      //CoeffVec = Solver_.Exponentiate(TryTimestep*TimeDirection_, MatElement);
      CoeffVec = Solver_.Exponentiate(TryTimestep*TimeDirection_);
      double TryCoeff = LinearAlgebra::norm_2(CoeffVec[CoeffVec.size()-1]);
      if (TryCoeff < 0.9 * (TryTimestep * ErrorBound_))
      {
         MinTimestep = TryTimestep;
         GuessTimestep = TryTimestep;
         Coeff = TryCoeff;
      }
      else if (TryCoeff > TryTimestep * ErrorBound_)
      {
         MaxTimestep = TryTimestep;
      }
      else
      {
         GuessTimestep = TryTimestep;
         Coeff = TryCoeff;
      }
   }

   std::cout << "Using timestep " << GuessTimestep << " for an error term of "
             << (Coeff / GuessTimestep) << " (per unit time).  Remainder term is " << Coeff << "\n";

   return GuessTimestep;
}

void KrylovLoop::Evolve(double Timestep, StatesInfo const& SInfo)
{
   Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(Timestep*TimeDirection_);
   TRACE(CoeffVec);
   TRACE("evolving")(Timestep*TimeDirection_);
   CenterWavefunction Psi = Solver_.ConstructExpansion(CoeffVec, SInfo);
   PsiNorm_ *= norm_frob(Psi.Center());
   Psi.normalize();
   Solver_ = KrylovSolver(Solver_.H, Solver_.H2, Psi);
   RealTime_ -= Timestep * TimeDirection_.imag();
   Beta_ -= Timestep * TimeDirection_.real();
}

void KrylovLoop::Evolve(std::complex<double> Timestep, StatesInfo const& SInfo)
{
   Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(Timestep);
   TRACE(CoeffVec);
   TRACE("evolving")(Timestep);
   CenterWavefunction Psi = Solver_.ConstructExpansion(CoeffVec, SInfo);
   PsiNorm_ *= norm_frob(Psi.Center());
   Psi.normalize();
   Solver_ = KrylovSolver(Solver_.H, Solver_.H2, Psi);
   RealTime_ -= Timestep.imag();
   Beta_ -= Timestep.real();
}

void KrylovLoop::Evolve(double Timestep, StatesInfo const& SInfo, SplitOperator const& NewH)
{
   Vector<std::complex<double> > CoeffVec = Solver_.Exponentiate(Timestep*TimeDirection_);
   CenterWavefunction Psi = Solver_.ConstructExpansion(CoeffVec, SInfo);
   PsiNorm_ *= norm_frob(Psi.Center());
   Psi.normalize();

   SplitOperator NewH2 = prod(NewH, NewH, NewH.TransformsAs());
   Solver_ = KrylovSolver(NewH, NewH2, Psi);

   RealTime_ -= Timestep * TimeDirection_.imag();
   Beta_ -= Timestep * TimeDirection_.real();
}

KrylovLoop::WavefunctionType KrylovLoop::Wavefunction() const
{
   WavefunctionType Ret = Solver_.Krylov[0] * PsiNorm_;
   Ret.Attributes()["Time"] = RealTime_;
   Ret.Attributes()["Beta"] = Beta_;
   return Ret;
}
