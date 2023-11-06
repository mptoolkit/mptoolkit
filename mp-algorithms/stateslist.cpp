// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/stateslist.cpp
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

#include "stateslist.h"
#include <sstream>
#include <iterator>
#include <iostream>
#include <cmath>

// StateParams

StateParams::StateParams()
   : NumStates(0),
     Wait(false),
     Test(false),
     Save(false),
     Variance(false),
     ZeroEnv(false)
{
}

StateParams
StateParams::ContinueFrom(StateParams const& Other)
{
   StateParams Result;
   Result.NumStates = Other.NumStates;
   return Result;
}

PStream::opstream& operator<<(PStream::opstream& out, StateParams const& Info)
{
   return out << Info.NumStates << Info.Wait << Info.Test << Info.Save << Info.Variance;
}

PStream::ipstream& operator>>(PStream::ipstream& in, StateParams& Info)
{
   return in >> Info.NumStates >> Info.Wait >> Info.Test >> Info.Save >> Info.Variance;
}

std::ostream& operator<<(std::ostream& out, StateParams const& Info)
{
   out << Info.NumStates;
   if (Info.Wait)
      out << 'w';
   if (Info.Test)
      out << 't';
   if (Info.Save)
      out << 's';
   if (Info.Variance)
      out << 'v';
   if (Info.ZeroEnv)
      out << 'z';
   return out;
}

// StatesList

StatesList::StatesList(char const* str)
{
   Append(str);
}

StatesList::StatesList(std::string const& str)
{
   Append(str.c_str());
}

StatesList::StatesList(int NumSweeps, int NumStates)
{
   StateParams i;
   i.NumStates = NumStates;
   Info = std::vector<StateParams>(NumSweeps, i);
}

// integer sign, returns +1, 0, or -1
int sign(int i)
{
   return int(i>0) - int(i<0);
}

void StatesList::AppendToken(char const* s)
{
   // Allowed notation:
   // 10   -> 10
   // 10x2 -> 10,10
   // 10+5x2 -> 15,20
   // 10..20x2 -> 15,20
   // 10+5..30
   // 10+5..30x10  // short-hand form for 10+5..30,30x9
   // 10*2..100x10 // short-hand form for 10*2..100,100x9
   //
   // Multiplicative increase:
   // 10^2 -> 10,10  (just for completeness)
   // 10*1.5^2 -> 15,23 (the half-open interval makes the arithmetic work here too!)
   // 10..30^2 -> 17,30 (multiplicatively increase from >10 to 30 in 2 sweeps; i.e. sweep n is 10*(30/10)^{n/2} )
   // 10*1.5..30 -> 15,23,30 (multiply by 1.5 each sweep, stopping at 30)
   // NOTE: since 2023-06-23, the initial number N is NOT included; the range is half-open

   // If there is only one sweep, we can append some flags, w,t,s,v although these arae currently unused
   // Even if there are multiple sweeps, we can add some flags,
   // z
   //
   // With the 10+5..30x10 and 10*2..100x10, adding a 'z' flag at the end only applies the 'z' to the final n sweeps 10, in the examples, so the number of sweeps with 100 states is actually 11).

   StateParams I;
   if (!Info.empty())
      I = StateParams::ContinueFrom(Info.back());

   int InitialSweepCount = Info.size();  // keep track of the number of sweeps so far, in case we need to go back and add flags

   int InitialStates = 0;
   int NumSweeps = 0;
   int FinalStates = 0;
   int Increment = 0;
   double Factor = 0.0;

   char* p = const_cast<char*>(s);
   // parse the leading number
   InitialStates = std::strtol(s, &p, 10);
   if (p == s)
   {
      PANIC("StatesList format error: number expected")(s);
   }
   s = p;

   if (s[0] == '.' && s[1] == '.')
   {
      // we have N..MxA notation
      s += 2;
      // parse the final number
      FinalStates = std::strtol(s, &p, 10);
      if (p == s)
      {
         PANIC("StatesList format error: N.. must be followed by another number")(s);
      }
      // and now we must get the "x"S part
      s = p;
      if (s[0] == 'x')
      {
         ++s;
         NumSweeps = std::strtol(s, &p, 10);
         if (p == s)
         {
            PANIC("StatesList format error: N..Mx must be followed by a number")(s);
         }
         s = p;
         // We have the form InitialStates..FinalStatesxNumSweeps
         for (int i = 0; i < NumSweeps; ++i)
         {
            I.NumStates = InitialStates + int((FinalStates-InitialStates) * (double(i+1)/NumSweeps) + 0.5);
            Info.push_back(I);
         }
      }
      else if (s[0] == '^')
      {
         ++s;
         NumSweeps = std::strtol(s, &p, 10);
         if (p == s)
         {
            PANIC("StatesList format error: N..M^ must be followed by a number")(s);
         }
         s = p;
         // We have the form InitialStates..FinalStates^NumSweeps
         for (int i = 0; i < NumSweeps; ++i)
         {
            I.NumStates = InitialStates + int(std::pow(double(FinalStates)/InitialStates, double(i+1)/NumSweeps) + 0.5);
            Info.push_back(I);
         }
      }
      else
      {
         PANIC("StatesList format error: N..M must be followed by x<number> or ^<number>")(s);
      }
   }
   else if (s[0] == '+')
   {
      // we have N+IxS or N+I..M notation
      ++s; // skip over the '+'
      // parse the increment
      Increment = std::strtol(s, &p, 10);
      if (p == s)
      {
         PANIC("StatesList format error: N+ must be followed by another number")(s);
      }
      s = p;
      if (s[0] == '.' && s[1] == '.')
      {
         // We have notation N+I..M
         s += 2;
         FinalStates = std::strtol(s, &p, 10);
         if (p == s)
         {
            PANIC("StatesList format error: N+I.. must be followed by another number")(s);
         }
         s = p;
         // we have the form InitialStates+Increment..FinalStates
         // check that Increment and (FinalStates-InitialStates) have the same sign
         int ss = sign(Increment);
         if (sign(FinalStates-InitialStates) != ss)
         {
            PANIC("StatesList format error: In N+I..M, the sign of the increment must match (M-N)")(s);
         }
         I.NumStates = InitialStates + Increment;
         // Take the sign into comparison in the loop.  If Increment is positive, then we want to loop while
         // states < FinalStates.  If Increment is negative, then we want to loop while states > FinalStates.
         while (I.NumStates*ss < FinalStates*ss)
         {
            Info.push_back(I);
            I.NumStates += Increment;
         }
         I.NumStates = FinalStates;
         Info.push_back(I);
         // Extension notation: we allow xS at the end, to keep iterating with the same number of states.
         if (s[0] == 'x')
         {
            // reset the sweep count, so flags only apply to subsequent sweeps
            InitialSweepCount = Info.size();
            ++s;
            NumSweeps = std::strtol(s, &p, 10);
            if (p == s)
            {
               PANIC("StatesList format error: N+I..Mx, expecting a number of sweeps")(s);
            }
            if (NumSweeps < 1)
            {
               PANIC("StatesList format error: number of sweeps must be >= 1");
            }
            for (int i = 0; i < NumSweeps; ++i)
            {
               Info.push_back(I);
            }
         }
      }
      else if (s[0] == 'x')
      {
         // we have notation N+IxS
         ++s; // skip over the 'x'
         NumSweeps = std::strtol(s, &p, 10);
         if (p == s)
         {
            PANIC("StatesList format error: N+Ix must be followed by a number")(s);
         }
         s = p;
         I.NumStates = InitialStates;
         for (int i = 0; i < NumSweeps; ++i)
         {
            I.NumStates += Increment;
            Info.push_back(I);
         }
      }
      else
      {
         PANIC("StatesList format error: N+I must be followed by ..M or xS")(s);
      }
   }
   else if (s[0] == '*')
   {
      // we have N*F^S or N*F..M notation
      ++s; // skip over the '*'
      // parse the factor
      Factor = std::strtod(s, &p);
      if (p == s)
      {
         PANIC("StatesList format error: N* must be followed by another number")(s);
      }
      s = p;
      if (Factor < 0.0)
      {
         PANIC("StatesList format error: negative scale factor is not allowed in N*F")(s);
      }
      // if the factor is an integer, then the first . will get eaten by the floating point number
      if (s[-1] == '.')
         --s;
      if (s[0] == '.' && s[1] == '.')
      {
         // We have notation N*I..M
         s += 2;
         FinalStates = std::strtol(s, &p, 10);
         if (p == s)
         {
            PANIC("StatesList format error: N*F.. must be followed by another number")(s);
         }
         s = p;
         // we have the form InitialStates*Factor..FinalStates
         I.NumStates = int(InitialStates*Factor + 0.5);
         int n = 1;
         if (FinalStates > InitialStates)
         {
            if (Factor <= 1.0)
            {
               PANIC("StatesList format error: N*F..M has M>N, but F is not > 1")(s);
            }
            while (I.NumStates < FinalStates)
            {
               Info.push_back(I);
               ++n;
               I.NumStates = int(InitialStates*std::pow(Factor, n) + 0.5);
            }
         }
         else if (FinalStates < InitialStates)
         {
            if (Factor >= 1.0)
            {
               PANIC("StatesList format error: N*F..M has M<N, but F is not < 1")(s);
            }
            while (I.NumStates > FinalStates)
            {
               Info.push_back(I);
               I.NumStates = int(InitialStates*std::pow(Factor, n) + 0.5);
               ++n;
            }
         }
         // The final case is FinalStates == InitialStates, in which case we just do 1 sweep with FinalStates
         // (which we do anyway)
         I.NumStates = FinalStates;
         Info.push_back(I);
         // Extension notation: we allow xS at the end, to keep iterating with the same number of states.
         if (s[0] == 'x')
         {
            // reset the sweep count, so flags only apply to subsequent sweeps
            InitialSweepCount = Info.size();
            ++s;
            NumSweeps = std::strtol(s, &p, 10);
            if (p == s)
            {
               PANIC("StatesList format error: N*F..Mx, expecting a number of sweeps")(s);
            }
            s = p;
            if (NumSweeps < 1)
            {
               PANIC("StatesList format error: number of sweeps must be >= 1");
            }
            for (int i = 0; i < NumSweeps; ++i)
            {
               Info.push_back(I);
            }
         }
      }
      else if (s[0] == '^')
      {
         // we have notation N*F^S
         ++s; // skip over the '^'
         NumSweeps = std::strtol(s, &p, 10);
         if (p == s)
         {
            PANIC("StatesList format error: N*F^S must be followed by a number")(s);
         }
         s = p;
         I.NumStates = InitialStates;
         for (int i = 0; i < NumSweeps; ++i)
         {
            I.NumStates = int(InitialStates*std::pow(Factor,i+1) + 0.5);
            Info.push_back(I);
         }
      }
      else
      {
         PANIC("StatesList format error: N*I must be followed by ..M or ^S")(s);
      }

   }
   else if (s[0] == 'x' || s[0] == '^')
   {
      // we can have MxS or M^S on its own
      ++s; // skip the character
      NumSweeps = std::strtol(s, &p, 10);
      if (p == s)
      {
         PANIC("StatesList format error: expected a number of states")(s);
      }
      s = p;
      I.NumStates = InitialStates;
      for (int i = 0; i < NumSweeps; ++i)
      {
         Info.push_back(I);
      }
   }
   else
   {
      // Parse any flags - these are only allowed if we have a single sweep
      while (s[0] != '\0')
      {
         switch (s[0])
         {
         case 'w' : I.Wait = true; I.Test = true; break;
         case 't' : I.Test = true;                break;
         case 'v' : I.Variance = true;            break;
         case 's' : I.Save = true;                break;
         case 'z' : I.ZeroEnv = true;             break;
         default  : PANIC("Unknown flag in StatesList")(s);
         }
         ++s;
      }
      I.NumStates = InitialStates;
      Info.push_back(I);
   }

   // Flags that can apply to multiple sweeps
   if (s[0] == 'z')
   {
      ++s;
      // zero the environment; this applies to all sweeps in this block
      for (int i = InitialSweepCount; i < Info.size(); ++i)
      {
         Info[i].ZeroEnv = true;
      }
   }
   if (s[0] != '\0')
   {
      PANIC("StatesList format error: unexpected extra characters")(s);
   }
}

void StatesList::Append(char const* s)
{
   char const* End = s;
   while (End[0] != '\0')
   {
      while (End[0] != '\0' && !isspace(*End) && End[0] != ',')
         ++End;

      // did we find a token?
      if (End != s)
         this->AppendToken(std::string(s,End).c_str());

      // skip the trailing whitespace
      while (isspace(*End) || End[0] == ',')
         ++End;

      // next token
      s = End;
   }
}

void
StatesList::Repeat(int n)
{
   while (n > 1)
   {
      Info.push_back(Info.back());
      --n;
   }
}

std::ostream& operator<<(std::ostream& out, StatesList const& States)
{
   out << "Number of sweeps: " << States.size() << '\n'
       << "Number of states at each sweep: ";
   for (int i = 0; i < States.size(); ++i)
   {
      out << States[i] << ' ';
   }
   out << '\n';
   return out;
}

PStream::opstream& operator<<(PStream::opstream& out, StatesList const& s)
{
   return out << s.Info;
}

PStream::ipstream& operator>>(PStream::ipstream& in, StatesList& s)
{
   return in >> s.Info;
}
