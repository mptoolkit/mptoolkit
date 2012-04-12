// -*- C++ -*- $Id$

#include "stateslist.h"
#include <sstream>
#include <iterator>
#include <iostream>

// StatesList

StatesList::StatesList(char const* str)
{
   Append(str);
}

StatesList::StatesList(std::string const& str)
{
   Append(str.c_str());
}

StatesList::~StatesList()
{ 
}

int const DefaultNumStates = 50000;

void StatesList::Append(char const* str)
{
   std::istringstream Stream(str);
   std::string S;
   while (Stream >> S)
   {
      InfoType I;
      I.Physics = false;
      I.Corrs = false;
      I.Wait = false;
      I.Save = false;
      I.Test = false;
      I.Broad = 0.0;
      I.Trunc = -1.0;
      I.MixFactor = -1.0;
      int Count = 1;
      int Increment = 0;
      // the first part is a number
      std::istringstream Token(S);
      Token >> I.NumStates;
      if (I.NumStates == 0)
	 I.NumStates = DefaultNumStates;
      char c;
      while (Token >> c)
      {
	 switch (c)
	 {
 	    case 'c' : I.Corrs = true;               break;
 	    case 'p' : I.Physics = true;             break;
 	    case 's' : I.Save = true;                break;
 	    case 'w' : I.Wait = true; I.Test = true; break;
 	    case 't' : I.Test = true;                break;
 	    case '*' : I.Wait = true;                break;
	    case 'b' : Token >> I.Broad;             break;
	    case 'r' : Token >> I.Trunc;             break;
            case 'f' : Token >> I.MixFactor;         break;
	    case 'x' : Token >> Count;               break;
	    case '+' : Token >> Increment;           break;
	 }
      }
      for (int i = 0; i < Count; ++i)
      {
         Info.push_back(I);
	 I.NumStates += Increment;
      }
   }
}

void StatesList::OverrideNoPhysics()
{
   for (std::size_t i = 0; i < Info.size(); ++i)
   {
      Info[i].Physics = false;
   }
}

void StatesList::OverrideNoCorrelations()
{
   for (std::size_t i = 0; i < Info.size(); ++i)
   {
      Info[i].Corrs = false;
   }
}

void StatesList::ShowOptions(std::ostream& out) const
{
   out << "Number of sweeps: " << NumSweeps() << '\n'
       << "Number of states at each sweep: ";
   for (int i = 0; i < NumSweeps(); ++i)
   {
      out << NumStates(i) << ' ';
   }

   bool AnyPhysics = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (CalculatePhysics(i)) AnyPhysics = true;
   if (AnyPhysics)
   {
      out << '\n'
          << "Calculating physics on sweeps ";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (CalculatePhysics(i)) 
            out << (i+1) << ' ';
      }
   }

   bool AnyCorrs = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (CalculateCorrelations(i)) AnyCorrs = true;
   if (AnyCorrs)
   {
      out << '\n'
          << "Calculating correlations on sweeps ";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (CalculateCorrelations(i)) 
            out << (i+1) << ' ';
      }
   }

   bool AnyTest = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (TestConverge(i)) AnyTest = true;
   if (AnyTest)
   {
      out << '\n'
          << "Test convergence on sweeps ";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (TestConverge(i)) 
            out << (i+1) << ' ';
      }
   }

   bool AnyWait = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (WaitConverge(i)) AnyWait = true;
   if (AnyWait)
   {
      out << '\n'
          << "Wait for convergence on sweeps ";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (WaitConverge(i)) 
            out << (i+1) << ' ';
      }
   }

   bool AnySave = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (SaveState(i)) AnySave = true;
   if (AnySave)
   {
      out << '\n'
          << "Save wavefunction on sweeps ";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (SaveState(i)) 
            out << (i+1) << ' ';
      }
   }

   bool AnyTrunc = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (TruncationError(i) != -1) AnyTrunc = true;
   if (AnyTrunc)
   {
      out << '\n'
          << "Changing the truncation error as\n";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (TruncationError(i) != -1.0) 
            out << "    Sweep " << (i+1) << " truncation error = " << TruncationError(i) << '\n';
      }
   }

   bool AnyBroad = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (Broadening(i) != 0) AnyBroad = true;
   if (AnyBroad)
   {
      out << '\n'
          << "Changing the broadening as\n";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (Broadening(i) != 0) 
            out << "    Sweep " << (i+1) << " broadening = " << Broadening(i) << '\n';
      }
   }
   
   bool AnyMixFactor = false;
   for (int i = 0; i < NumSweeps(); ++i)
      if (MixFactor(i) != -1) AnyMixFactor = true;
   if (AnyMixFactor)
   {
      out << '\n'
          << "Changing the mix factor as\n";
      for (int i = 0; i < NumSweeps(); ++i)
      {
         if (MixFactor(i) != -1) 
            out << "    Sweep " << (i+1) << " mix factor = " << MixFactor(i) << '\n';
      }
   }

   out << '\n';
}

bool StatesList::AnyPhysics() const 
{
   for (int i = 0; i < NumSweeps(); ++i)
   {
      if (CalculatePhysics(i)) return true;
   }
   return false;
}

bool StatesList::AnyCorrelations() const 
{
   for (int i = 0; i < NumSweeps(); ++i)
   {
      if (CalculateCorrelations(i)) return true;
   }
   return false;
}

PStream::opstream& operator<<(PStream::opstream& out, StatesList::InfoType const& i)
{
   return out << i.NumStates << i.Physics << i.Corrs << i.Wait << i.Test << i.Save 
              << i.Trunc << i.Broad << i.MixFactor;
}

PStream::ipstream& operator>>(PStream::ipstream& in, StatesList::InfoType& i)
{
   return in >> i.NumStates >> i.Physics >> i.Corrs >> i.Wait >> i.Test >> i.Save 
             >> i.Trunc >> i.Broad >> i.MixFactor;
}

PStream::opstream& operator<<(PStream::opstream& out, StatesList const& s)
{
   return out << s.Info;
}

PStream::ipstream& operator>>(PStream::ipstream& in, StatesList& s)
{
   return in >> s.Info;
}

// WeightsList

WeightsList:: WeightsList(char const* str)
{
   std::istringstream Stream(str);
   double W;
   while (Stream >> W)
   {
      Weights.push_back(W);
   }
}

void WeightsList::ShowOptions(std::ostream& out) const
{
   out << "Number of wavefunctions to calculate: " << this->NumWeights() 
       << "\nWeights: ";
   std::copy(this->Weights.begin(), this->Weights.end(), std::ostream_iterator<double>(out, " "));
   out << '\n';
}

PStream::opstream& operator<<(PStream::opstream& out, WeightsList const& s)
{
   return out << s.Weights;
}

PStream::ipstream& operator>>(PStream::ipstream& in, WeightsList& s)
{
   return in >> s.Weights;
}

StatesInfo StatesList::GetStatesInfo(StatesInfo const& Default, int s) const
{
   StatesInfo Result = Default;
   Result.MaxStates = this->NumStates(s);
   Result.TruncationCutoff = this->TruncationError(s);
   return Result;
}
