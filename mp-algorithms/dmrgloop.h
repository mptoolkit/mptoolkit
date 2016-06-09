// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/dmrgloop.h
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

#include "common/statistics.h"
#include "common/messagelogger.h"
#include "pstream/pstream.h"
#include "common/conflist.h"
#include "stateslist.h"
#include "common/proccontrol.h"
#include "matrixproduct/mpwavefunction.h"

template <typename SolverType>
class DMRGLoop;

template <typename SolverType>
PStream::opstream& operator<<(PStream::opstream& out, DMRGLoop<SolverType> const& s);

template <typename SolverType>
PStream::ipstream& operator>>(PStream::ipstream& in, DMRGLoop<SolverType>& s);



template <typename SolverType>
class DMRGLoop
{
   public:
      typedef typename SolverType::WavefunctionType WavefunctionType;

      DMRGLoop() {}

      DMRGLoop(SolverType const& Solver, StatesList const& States);

      void Run(std::string const& BasePathFile, ConfList const& Conf);

      bool Finished() const;

      void SetPreviousCPUTime(double n) { PreviousCPUTime_ = n; }
      void SetPreviousElapsedTime(double n) { PreviousElapsedTime_ = n; }

      double PreviousCPUTime() const { return PreviousCPUTime_; }
      double PreviousElapsedTime() const { return PreviousElapsedTime_; }

   private:
      void DoIteration();

      void SaveWavefunction() const;

   private:
      StatesList States_;
      SolverType Solver_;
      int SweepRecNum_;
      int SweepNum_;
      int Direction_;  // +1 for right-moving, -1 for left moving
      bool DoneShift_;  // true if we have done the shift for the next iteration, false otherwise
      double PreviousCPUTime_;
      double PreviousElapsedTime_;

      // The following data is not persistent, but loaded from the configuration 
      std::string BinPath_;
      std::string BasePath_;
      std::string FileName_;
      int NumIter_;
      int MinStates_;
      int MaxStates_;
      bool ShowOutput_;
      double CFactor_;
      bool TwoSite_;
      StatesInfo DefaultStates_;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, DMRGLoop<SolverType> const& s);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, DMRGLoop<SolverType>& s);
};

template <typename SolverType>
PStream::opstream& operator<<(PStream::opstream& out, DMRGLoop<SolverType> const& s)
{
   return out << s.States_ 
              << s.Solver_ 
              << s.SweepRecNum_
              << s.SweepNum_
              << s.Direction_
              << s.DoneShift_
              << s.PreviousCPUTime_
              << s.PreviousElapsedTime_;
}

template <typename SolverType>
PStream::ipstream& operator>>(PStream::ipstream& in, DMRGLoop<SolverType>& s)
{
   return in >> s.States_ 
             >> s.Solver_ 
             >> s.SweepRecNum_
             >> s.SweepNum_
             >> s.Direction_
             >> s.DoneShift_
             >> s.PreviousCPUTime_
             >> s.PreviousElapsedTime_;
}

template <typename SolverType>
DMRGLoop<SolverType>::DMRGLoop(SolverType const& Solver, StatesList const& States)
   : States_(States), 
     Solver_(Solver), 
     SweepRecNum_(0),
     SweepNum_(0),
     Direction_(1),
     DoneShift_(false),
     PreviousCPUTime_(0),
     PreviousElapsedTime_(0)
{
   if (Solver_.LeftSize() > Solver_.RightSize()) Direction_ = -1;

   // start the first sweep
   if (!this->Finished()) 
   {
      double NewMix = States_.MixFactor(SweepRecNum_);
      if (NewMix != -1)
         CFactor_ = NewMix;
      Solver_.StartSweep(true, States_.Broadening(SweepRecNum_));
   }
}

template <typename SolverType>
void DMRGLoop<SolverType>::Run(std::string const& BasePathFile, ConfList const& Conf)
{
   NumIter_ = Conf.Get("NumIterations", 10);
   msg_log(2, "InfoLog") << "NumIterations = " << NumIter_ << '\n';

   ShowOutput_ = Conf.Get("ShowOutput", true);

   std::tie(BasePath_, FileName_) = pheap::SplitPathFile(BasePathFile);
   if (BasePath_.empty()) BasePath_ = "./";
   BinPath_ = Conf.Get("BinPath", std::string("."));
   if (BinPath_[BinPath_.size()-1] != '/') BinPath_ += '/';

   CFactor_ = Conf.Get("MixFactor", 0.0);
   TwoSite_ = Conf.Get("TwoSite", false);
   MinStates_ = Conf.Get("MinStates", 0);
   MaxStates_ = Conf.Get("MaxStates", DefaultMaxStates);

   DefaultStates_.MinStates = MinStates_;
   DefaultStates_.MaxStates = MaxStates_;

   Solver_.RestoreLogFiles(BasePathFile, Conf);

   while (!this->Finished())
      this->DoIteration();
}

template <typename SolverType>
bool DMRGLoop<SolverType>::Finished() const
{
   return SweepRecNum_ == States_.NumSweeps();
}

template <typename SolverType>
void DMRGLoop<SolverType>::SaveWavefunction() const
{
   std::string FileName = BasePath_ + FileName_
      + ".psi."+boost::lexical_cast<std::string>(SweepRecNum_+1);
   pvalue_ptr<MPWavefunction> OutPsi(new LinearWavefunction(Solver_.Wavefunction().AsLinearWavefunction()));
   pheap::ExportHeap(FileName, OutPsi);
}

template <typename SolverType>
void DMRGLoop<SolverType>::DoIteration()
{
   PRECONDITION(!this->Finished());

   std::cout.precision(12);

   if (Direction_ == 1)  // moving to the right?
   {
      if (!DoneShift_)
      {
         if (Solver_.LeftSize() == 1 && States_.TestConverge(SweepRecNum_))
            Solver_.PrepareConvergenceTest();

         Solver_.StartIteration();

         Solver_.ShiftRightAndExpand();

         // special case for the boundary, make sure we expand the boundary sites.
         if (TwoSite_ || Solver_.RightSize() == 1) Solver_.ExpandRight();
         DoneShift_ = true;
      }

      ProcControl::TestAsyncCheckpoint();

      double E = Solver_.Solve(NumIter_);

      TruncationInfo States = Solver_.TruncateLeft(States_.GetStatesInfo(DefaultStates_, SweepRecNum_), CFactor_);

      Solver_.EndIteration();

      if (ShowOutput_)
      {
         std::cout << '(' << Solver_.LeftSize() << ',' << Solver_.RightSize()
                   << ") " << E << ' ' << States.KeptStates() << ' ' << States.TruncationError() << '\n';
      }

      // Are we at the end of the sweep?
      if (Solver_.RightSize() == 1) 
      {
         Solver_.EndSweep();

         if (States_.SaveState(SweepRecNum_) 
             && (!States_.WaitConverge(SweepRecNum_) || Solver_.IsConverged()))
         {
            this->SaveWavefunction();
         }
         
	 bool DoneSweep = !States_.WaitConverge(SweepRecNum_) || Solver_.IsConverged();
	 if (DoneSweep) ++SweepRecNum_;

         if (!this->Finished())
         {
            double NewMix = States_.MixFactor(SweepRecNum_);
            if (NewMix != -1)
               CFactor_ = NewMix;
            Solver_.StartSweep(DoneSweep, States_.Broadening(SweepRecNum_));
         }

         Direction_ = -1;
      }
   }
   else
   {
      // Direction_ == -1, moving to the left
      if (!DoneShift_)
      {
         if (Solver_.RightSize() == 1 && States_.TestConverge(SweepRecNum_))
            Solver_.PrepareConvergenceTest();

         Solver_.StartIteration();

         Solver_.ShiftLeftAndExpand();

         // special case for the boundary, make sure we expand the boundary sites.
         if (TwoSite_ || Solver_.LeftSize() == 1) Solver_.ExpandLeft();
         DoneShift_ = true;
      }

      ProcControl::TestAsyncCheckpoint();

      double E = Solver_.Solve(NumIter_);

      TruncationInfo States = Solver_.TruncateRight(States_.GetStatesInfo(DefaultStates_, SweepRecNum_), CFactor_);

      Solver_.EndIteration();

      if (ShowOutput_)
      {
         std::cout << '(' << Solver_.LeftSize() << ',' << Solver_.RightSize()
                   << ") " << E << ' ' << States.KeptStates() << ' ' << States.TruncationError() << '\n';
      }

      // Are we at the end of the sweep?
      if (Solver_.LeftSize() == 1) 
      {
         Solver_.EndSweep();

         if (States_.SaveState(SweepRecNum_) 
             && (!States_.WaitConverge(SweepRecNum_) || Solver_.IsConverged()))
         {
            this->SaveWavefunction();
         }
         
         if (!States_.WaitConverge(SweepRecNum_) || Solver_.IsConverged())
            ++SweepRecNum_;

         if (!this->Finished())
         {
            double NewMix = States_.MixFactor(SweepRecNum_);
            if (NewMix != -1)
               CFactor_ = NewMix;
            Solver_.StartSweep(!States_.WaitConverge(SweepRecNum_) || Solver_.IsConverged(),
			       States_.Broadening(SweepRecNum_));
         }

         Direction_ = 1;
      }
   }

   DoneShift_ = false;
}
