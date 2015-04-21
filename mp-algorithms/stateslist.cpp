// -*- C++ -*- $Id$

#include "stateslist.h"
#include <sstream>
#include <iterator>
#include <iostream>

// StateParams

StateParams::StateParams()
   : NumStates(0),
     Wait(false),
     Test(false),
     Save(false),
     Variance(false)
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

void StatesList::AppendToken(char const* s)
{
   // notation is (symbols in quotes are literals, capitals are numbers)
   // N[[".."M | "+"A] "x"A]["w"]["t"]

   StateParams I;
   if (!Info.empty())
      I = StateParams::ContinueFrom(Info.back());

   int NumSweeps = 1;
   int FinalStates = 0;
   int Increment = 0;

   char* p = const_cast<char*>(s);
   // parse the leading number
   I.NumStates = std::strtol(s, &p, 10);
   if (p == s)
   {
      PANIC("Did not find a number in the StatesList")(s);
   }
   s = p;

   // see what we have next
   if (s[0] == '\0')
   {
      // we have the complete info
      Info.push_back(I);
      return;
   }
   else if (s[0] == '.' && s[1] == '.')
   {
      // we have N..MxA notation
      s += 2;
      // parse the final number
      FinalStates =  std::strtol(s, &p, 10);
      if (p == s)
      {
	 PANIC("Did not find a final number in N..MxS notation in the StatesList")(s);
      }
      // and now we must get the "x"S part
      s = p;
      if (s[0] != 'x')
      {
	 PANIC("N..M must be followed by x<number> in StatesList");
      }
      ++s;
      NumSweeps = std::strtol(s, &p, 10);
      if (p == s)
      {
	 PANIC("N..M must be followed by x<number> in StatesList");
      }
      s = p;
   }
   else if (s[0] == '+')
   {
      // we have N+MxA notation
      ++s; // skip over the '+'
      // parse the increment
      Increment = std::strtol(s, &p, 10);
      if (p == s)
      {
	 PANIC("Did not find an increment in N+M notation in the StatesList")(s);
      }
      // and now we must get the "x"S part
      s = p;
      if (s[0] != 'x')
      {
	 PANIC("N+M must be followed by x<number> in StatesList");
      }
      ++s;
      NumSweeps = std::strtol(s, &p, 10);
      if (p == s)
      {
	 PANIC("N+M must be followed by x<number> in StatesList");
      }
      s = p;
   }
   else if (s[0] == 'x')
   {
      // we can have MxA on its own
      ++s;
      NumSweeps = std::strtol(s, &p, 10);
      if (p == s)
      {
	 PANIC("'x' must be followed by <number> in StatesList");
      }
      s = p;
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
	 default  : PANIC("Unknown flag in StatesList")(s);
	 }
	 ++s;
      }
   }

   if (s[0] != '\0')
   {
      PANIC("Extra characters at end of NumStates line")(s);
   }

   int InitialStates = I.NumStates;
   // Now put it together if we have repeats
   for (int i = 0; i < NumSweeps; ++i)
   {
      if (FinalStates != 0)
      {
	 I.NumStates = InitialStates + i*double((FinalStates-InitialStates)/double(NumSweeps-1));
      }
      else if (Increment != 0)
      {
	 I.NumStates = InitialStates + i*Increment;
      }
      Info.push_back(I);
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
